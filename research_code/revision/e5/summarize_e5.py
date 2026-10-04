"""Validated E5 summaries with staged rendering and rollback publication.

See README_SUMMARY_VALIDATION.md for the input contract and transaction limits.
No model fitting, input/protocol edits, or new reference draws occur here.
"""
from pathlib import Path
import argparse
import hashlib
import io
import json
import math
import os
import shutil
import tempfile

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PHENOTYPES = ('1', '2', '4', '5')
ROLES = ('calibration', 'evaluation')
LEVELS = (0, 400, 1000, 2500)
REPLICATES = tuple(range(100))
FOLDS, ARMS, SAMPLE_N = 5, 2, 599
DIFFERENCE_ATOL, DIFFERENCE_RTOL = 1e-12, 1e-10
INTERVAL_TOL = 1e-12
OUTPUT_NAMES = ('operating_rates.csv', 'bootstrap_draws.npz',
                'Figure_S1_authentic_wheat.pdf', 'Figure_S1_authentic_wheat.svg',
                'Figure_S1_authentic_wheat.png', 'summary.json')
LOCK_NAME = '.e5-summary.lock'
JOURNAL_NAME = '.e5-summary-transaction.json'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json_object(raw, label):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f'{label}: duplicate JSON key {key!r}')
            result[key] = value
        return result
    def invalid_constant(value):
        raise ValueError(f'{label}: non-finite JSON constant {value}')
    try:
        result = json.loads(raw, object_pairs_hook=pairs, parse_constant=invalid_constant)
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise ValueError(f'{label}: invalid JSON') from exc
    if not isinstance(result, dict):
        raise ValueError(f'{label}: expected a JSON object')
    return result


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _recorded_integer(manifest, key, expected):
    value = manifest.get(key)
    _require(type(value) is int and value == expected,
             f'run_manifest {key}: expected integer {expected}, got {value!r}')


def _seed(role, replicate, stream):
    value = f'patterns-v14-e5|20260908|{role}|{replicate}|{stream}'
    return int.from_bytes(hashlib.sha256(value.encode()).digest()[:8], 'little') % 2**32


def validate_inputs(root):
    """Read/validate one byte snapshot, before statistics or output writes.

    Model-fit counts are cross-checked against the completed score grid and the
    recorded manifest (two arms times five fits per OOF score). This does not
    independently reconstruct the underlying model fits from raw data.
    """
    root = Path(root)
    try:
        manifest_bytes = (root / 'run_manifest.json').read_bytes()
    except OSError as exc:
        raise ValueError('run_manifest.json is missing or unreadable') from exc
    manifest = _json_object(manifest_bytes, 'run_manifest')
    _require(manifest.get('status') == 'COMPLETE', 'run_manifest status must be COMPLETE')
    _require(isinstance(manifest.get('failed'), list) and not manifest['failed'],
             'run_manifest failed must be an empty list')

    # These files are the unchanged, recorded reader exposition and its identity.
    protocol_dir = Path(__file__).resolve().parents[3] / 'protocol/e5'
    try:
        freeze = _json_object((protocol_dir / 'freeze.json').read_bytes(), 'protocol freeze')
        protocol_hash = sha(protocol_dir / 'scientific_protocol.md')
    except OSError as exc:
        raise ValueError('Bundled E5 protocol/freeze is missing or unreadable') from exc
    _require(freeze.get('scientific_parameters_changed') is False,
             'E5 reader freeze must record unchanged scientific parameters')
    _require(protocol_hash == freeze.get('protocol_sha256'), 'Bundled E5 protocol checksum mismatch')
    _require(manifest.get('protocol_sha256') in
             (freeze.get('source_protocol_sha256'), protocol_hash),
             'run_manifest protocol_sha256 is not a recorded E5 protocol identity')
    _require(manifest.get('raw_sha256') == freeze.get('raw_sha256'),
             'run_manifest raw_sha256 differs from the recorded E5 input')
    _require(manifest.get('code_sha256') == sha(Path(__file__).with_name('run_e5.py')),
             'run_manifest fitting-code checksum differs from unchanged run_e5.py')
    try:
        score_bytes = (root / 'scores.csv').read_bytes()
    except OSError as exc:
        raise ValueError('scores.csv is missing or unreadable') from exc
    score_hash = hashlib.sha256(score_bytes).hexdigest()
    _require(score_hash == manifest.get('scores_sha256'), 'scores.csv checksum mismatch')
    try:
        d = pd.read_csv(io.BytesIO(score_bytes), float_precision='round_trip',
                        dtype={'phenotype': str, 'reference_id': str, 'role': str, 'status': str})
    except Exception as exc:
        raise ValueError('scores.csv cannot be parsed') from exc
    required = {'reference_id', 'role', 'replicate', 'phenotype', 'L', 'R2_B', 'R2_C',
                'T_scope', 'fold_seed', 'candidate_seed', 'n', 'model', 'status'}
    _require(required <= set(d.columns), f'scores.csv missing columns: {sorted(required-set(d.columns))}')
    _require(d['status'].eq('COMPLETE').all(), 'Every score row status must be COMPLETE')
    _require(set(d.phenotype) == set(PHENOTYPES), 'Phenotype set must be exactly 1,2,4,5')
    _require(set(d.role) == set(ROLES), 'Role set must be calibration and evaluation')
    for col in ['replicate', 'L', 'n', 'fold_seed', 'candidate_seed']:
        try:
            values = pd.to_numeric(d[col], errors='raise').to_numpy(dtype=float)
        except (ValueError, TypeError) as exc:
            raise ValueError(f'{col} must contain finite integers') from exc
        _require(np.isfinite(values).all() and (values == np.floor(values)).all(),
                 f'{col} must contain finite integers')
        _require((np.abs(values) <= 2**53).all(), f'{col} exceeds the exact integer range')
        d[col] = values.astype(np.int64)
    _require(d['n'].eq(SAMPLE_N).all(), f'Every score row must report n={SAMPLE_N}')
    _require(d['model'].eq('ridge_alpha1').all(), 'Model must be ridge_alpha1')
    _require(d.replicate.isin(REPLICATES).all(), 'Replicate IDs must be integers 0..99')
    expected_ids = d.role + '_' + d.replicate.map(lambda i: f'{i:03d}')
    _require(d.reference_id.eq(expected_ids).all(), 'reference_id must match role and replicate')
    cal_ids = set(d.loc[d.role == 'calibration', 'reference_id'])
    eval_ids = set(d.loc[d.role == 'evaluation', 'reference_id'])
    _require(cal_ids.isdisjoint(eval_ids), 'Calibration and evaluation reference IDs overlap')
    _require(not d.duplicated(['reference_id', 'phenotype', 'L']).any(),
             'Duplicate reference_id/phenotype/L score')
    expected = {(role, i, ph, level)
                for role in ROLES for i in REPLICATES for ph in PHENOTYPES
                for level in ((0,) if role == 'calibration' else LEVELS)}
    actual = set(d[['role', 'replicate', 'phenotype', 'L']].itertuples(index=False, name=None))
    _require(actual == expected and len(d) == len(expected),
             f'Incomplete or invalid phenotype/role/L/reference grid: '
             f'{len(expected-actual)} missing, {len(actual-expected)} unexpected')
    for row in d[['role', 'replicate', 'fold_seed', 'candidate_seed']].drop_duplicates().itertuples(index=False):
        _require(row.fold_seed == _seed(row.role, row.replicate, 'folds') and
                 row.candidate_seed == _seed(row.role, row.replicate, 'candidates'),
                 'Recorded seeds disagree with the unchanged E5 seed rule')
    for col in ['R2_B', 'R2_C', 'T_scope']:
        try:
            d[col] = pd.to_numeric(d[col], errors='raise')
        except (ValueError, TypeError) as exc:
            raise ValueError(f'{col} must be finite numeric values') from exc
        _require(np.isfinite(d[col].to_numpy(dtype=float)).all(), f'{col} contains non-finite values')
    _require(np.allclose(d.T_scope, d.R2_C - d.R2_B, rtol=DIFFERENCE_RTOL, atol=DIFFERENCE_ATOL),
             'T_scope differs from R2_C - R2_B beyond atol=1e-12, rtol=1e-10')
    score_n = len(d)
    reference_n = d.reference_id.nunique()
    fit_n = score_n * FOLDS * ARMS
    for key, value in {'reference_n': reference_n, 'phenotype_n': d.phenotype.nunique(),
                       'planned_score_rows': len(expected), 'successful_score_rows': score_n,
                       'planned_model_fits': len(expected) * FOLDS * ARMS,
                       'completed_model_fits': fit_n}.items():
        _recorded_integer(manifest, key, int(value))
    counts = {'reference_n': int(reference_n), 'score_n': score_n, 'fit_n': fit_n,
              'failed_n': int((d.status != 'COMPLETE').sum()) + len(manifest['failed']),
              'excluded_n': len(expected-actual),
              'manifest_sha256': hashlib.sha256(manifest_bytes).hexdigest(), 'scores_sha256': score_hash}
    return d, counts


def wilson(k, n):
    _require(type(k) is int and type(n) is int and n > 0 and 0 <= k <= n,
             'Wilson counts must be integers with 0 <= k <= n and n > 0')
    z=1.959963984540054; p=k/n; den=1+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return center-half, center+half


def validate_interval(lo, hi, lower_bound, upper_bound, label):
    lo, hi = np.asarray(lo, dtype=float), np.asarray(hi, dtype=float)
    _require(np.isfinite(lo).all() and np.isfinite(hi).all(), f'{label}: non-finite interval')
    _require((lo <= hi).all(), f'{label}: reversed interval')
    _require((lo >= lower_bound-INTERVAL_TOL).all() and (hi <= upper_bound+INTERVAL_TOL).all(),
             f'{label}: interval outside permitted bounds')


def interval_yerr(y, lo, hi):
    """Only normalize roundoff <=1e-12 after validating genuine rate intervals."""
    y, lo, hi = (np.asarray(a, dtype=float) for a in (y, lo, hi))
    _require(np.isfinite(y).all() and ((0 <= y) & (y <= 1)).all(), 'Non-finite or invalid rate')
    validate_interval(lo, hi, 0, 1, 'Wilson rate')
    err = np.array([y-lo, hi-y])
    _require((err >= -INTERVAL_TOL).all(), 'Wilson interval does not contain its observed rate')
    return np.where(err < 0, 0.0, err)


def compute_statistics(d, counts):
    # Identical draw order, rank rule, arithmetic and CSV column order to v14.
    n_cal = d.query('role=="calibration" and phenotype=="1"').reference_id.nunique()
    n_eval = d.query('role=="evaluation" and phenotype=="1" and L==0').reference_id.nunique()
    rng=np.random.default_rng(2026090814)
    ci_ids=rng.integers(0,n_cal,size=(2000,n_cal))
    ei_ids=rng.integers(0,n_eval,size=(2000,n_eval))
    rank=math.ceil((n_cal+1)*.95)
    rows=[]; boot={}
    for ph in PHENOTYPES:
        p=d[d.phenotype==ph]
        cal=p.query('role=="calibration"').sort_values('replicate').T_scope.to_numpy()
        ev=p.query('role=="evaluation"').pivot(index='replicate',columns='L',values='T_scope').sort_index()
        tau=np.sort(cal)[rank-1]; btau=np.sort(cal[ci_ids],axis=1)[:,rank-1]
        a=ev.to_numpy(); bev=a[ei_ids]
        br=(bev>btau[:,None,None]).mean(axis=1); bf=(bev>.05).mean(axis=1)
        boot['ph'+ph+'_threshold']=btau; boot['ph'+ph+'_rate']=br; boot['ph'+ph+'_difference']=br-bf
        for j,L in enumerate(ev.columns):
            vals=a[:,j]; n=len(vals); k=int((vals>tau).sum()); f=int((vals>.05).sum())
            lo,hi=wilson(k,n); flo,fhi=wilson(f,n)
            ql,qh=np.quantile(br[:,j],[.025,.975],method='linear')
            dl,dh=np.quantile((br-bf)[:,j],[.025,.975],method='linear')
            interval_yerr([k/n],[lo],[hi]); interval_yerr([f/n],[flo],[fhi])
            validate_interval(ql,qh,0,1,'Calibration-bootstrap rate')
            validate_interval(dl,dh,-1,1,'Paired-bootstrap difference')
            group=p.query('role=="evaluation" and L==@L')
            valid_n=int(group.status.eq('COMPLETE').sum())
            failed_n=int(group.status.ne('COMPLETE').sum())
            rows.append(dict(phenotype=ph,L=int(L),threshold=tau,calibration_n=len(cal),calibration_rank=rank,attempted_n=len(group),valid_n=valid_n,failed_n=failed_n,excluded_n=counts['excluded_n'],exceedance_k=k,rate=k/n,ci_lower=lo,ci_upper=hi,fixed_threshold=.05,fixed_k=f,fixed_rate=f/n,fixed_ci_lower=flo,fixed_ci_upper=fhi,calibration_bootstrap_rate_lower=ql,calibration_bootstrap_rate_upper=qh,rate_difference=(k-f)/n,rate_difference_bootstrap_lower=dl,rate_difference_bootstrap_upper=dh,mean_R2_B=float(group.R2_B.mean()),mean_R2_C=float(group.R2_C.mean()),additional_action_over_correct_provenance=False))
    table=pd.DataFrame(rows)
    _require(np.isfinite(table.select_dtypes(include=np.number).to_numpy()).all(), 'Non-finite computed table value')
    _require(all(np.isfinite(a).all() for a in boot.values()), 'Non-finite bootstrap value')
    summary={'status':'COMPLETE','phenotypes_reported':list(PHENOTYPES),'raw_outcomes_authentic':True,'scope':'conditional split/candidate randomization on fixed X,y; not biological-null calibration','reference_realizations':counts['reference_n'],'score_rows':counts['score_n'],'model_fits':counts['fit_n'],'failures':counts['failed_n'],'additional_practitioner_utility_established':False,'reason':'Correct provenance inspection requires repairing C at every L, irrespective of either numeric screen. Injected candidates and randomized folds supply no new untouched deployment decision.'}
    return table, boot, summary


def render_outputs(stage, table, boot, summary):
    table.to_csv(stage/'operating_rates.csv',index=False,float_format='%.17g')
    np.savez_compressed(stage/'bootstrap_draws.npz',**boot)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(2,2,figsize=(7.1,6.1),sharex=True,sharey=True)
    try:
        for ax,ph in zip(axes.flat,PHENOTYPES):
            q=table[table.phenotype==ph]; x=np.arange(len(LEVELS))
            for name,color,label in [('','darkcyan','Conditional reference cutoff'),('fixed_','#AE6040','Analyst-defined 0.05')]:
                y=q[name+'rate'].to_numpy(); lo=q[name+'ci_lower'].to_numpy(); hi=q[name+'ci_upper'].to_numpy()
                ax.errorbar(x,y,yerr=interval_yerr(y,lo,hi),color=color,marker='o',capsize=3,label=label)
            ax.set_title(f'Yield environment {ph}; cutoff {q.threshold.iloc[0]:.4f}')
            ax.set(xticks=x,xticklabels=list(LEVELS),ylim=(-.03,1.05)); ax.spines[['top','right']].set_visible(False)
        axes[1,0].set_xlabel('Additional candidates L'); axes[1,1].set_xlabel('Additional candidates L')
        axes[0,0].set_ylabel('Conditional exceedance'); axes[1,0].set_ylabel('Conditional exceedance')
        fig.legend(*axes[0,0].get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False)
        fig.tight_layout(rect=(0,.05,1,1))
        for ext in ['pdf','svg','png']:
            fig.savefig(stage/f'Figure_S1_authentic_wheat.{ext}',dpi=600,bbox_inches='tight')
    finally:
        plt.close(fig)
    (stage/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


def _write_json(path, value):
    with path.open('w') as f:
        json.dump(value, f, indent=2)
        f.write('\n'); f.flush(); os.fsync(f.fileno())


def publish_outputs(root, stage):
    """Per-file atomic replacements with backups and catchable-error rollback.

    NOT a group-atomic snapshot for concurrent readers or abrupt termination.
    A journal/lock makes an interrupted transaction visible and fails closed.
    """
    journal=root/JOURNAL_NAME
    backup=stage/'backup'; backup.mkdir()
    existed={}
    for name in OUTPUT_NAMES:
        target=root/name
        _require(not target.is_symlink() and (not target.exists() or target.is_file()),
                 f'Output destination must be a regular file or absent: {name}')
        existed[name]=target.exists()
        if target.exists(): shutil.copy2(target,backup/name)
        _require((stage/name).is_file() and (stage/name).stat().st_size > 0,
                 f'Staged output is absent or empty: {name}')
    record={'status':'PREPARED','stage_directory':stage.name,'output_order':list(OUTPUT_NAMES),
            'previous_files':{n:{'existed':existed[n],'sha256':sha(backup/n) if existed[n] else None} for n in OUTPUT_NAMES},
            'staged_sha256':{n:sha(stage/n) for n in OUTPUT_NAMES},
            'semantics':'Per-file atomic replacement; summary last; catchable errors roll back. Not group atomic or power-loss safe.'}
    _write_json(journal,record)
    try:
        for name in OUTPUT_NAMES:
            os.replace(stage/name,root/name)
        journal.unlink()
    except BaseException as commit_error:
        errors=[]
        for name in reversed(OUTPUT_NAMES):
            try:
                if existed[name]:
                    restore=stage/('restore-'+name)
                    shutil.copy2(backup/name,restore)
                    os.replace(restore,root/name)
                elif (root/name).exists():
                    (root/name).unlink()
            except BaseException as exc:
                errors.append(f'{name}: {type(exc).__name__}: {exc}')
        if errors:
            record.update(status='ROLLBACK_FAILED',rollback_errors=errors)
            try: _write_json(journal,record)
            except OSError: pass  # The prepared journal still identifies backups.
            raise RuntimeError('E5 commit/rollback failed; retain journal and stage backups for manual recovery') from commit_error
        journal.unlink()
        raise RuntimeError('E5 output commit failed; previous output set restored') from commit_error


def summarize(root):
    root=Path(root)
    d,counts=validate_inputs(root)
    _require(not (root/JOURNAL_NAME).exists(), 'Unresolved E5 transaction journal; recover previous outputs before rerunning')
    lock=root/LOCK_NAME
    try:
        fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600)
    except FileExistsError as exc:
        raise ValueError('E5 summary writer lock exists; concurrent writer or interrupted run requires inspection') from exc
    stage=None
    try:
        with os.fdopen(fd,'w') as f: f.write(f'pid={os.getpid()}\n')
        table,boot,summary=compute_statistics(d,counts)
        stage=Path(tempfile.mkdtemp(prefix='.e5-summary-stage-',dir=root))
        render_outputs(stage,table,boot,summary)
        publish_outputs(root,stage)
    finally:
        if not (root/JOURNAL_NAME).exists():
            if stage is not None: shutil.rmtree(stage)
            lock.unlink(missing_ok=True)
    print(table[['phenotype','L','threshold','exceedance_k','fixed_k']].to_string(index=False))
    return table


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',required=True)
    args=parser.parse_args()
    summarize(args.root)
