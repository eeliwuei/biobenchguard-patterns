"""Regression and fault-injection checks; no biological model refitting.

Default fixtures use the bundled canonical derived score/operating tables plus
recorded v14 manifest/bootstrap fixtures. E5_REFERENCE_ROOT can instead point to
an original completed E5 directory; tests always copy inputs to temporary paths.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / 'summarize_e5.py'
PACKAGE = SCRIPT.parents[3]
FIXTURES = Path(__file__).parent / 'fixtures'
REFERENCE = Path(os.environ.get('E5_REFERENCE_ROOT', PACKAGE / 'results_summary/e5'))
spec = importlib.util.spec_from_file_location('e5_summary_under_test', SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.fixture
def run_root(tmp_path):
    root = tmp_path / 'run'
    root.mkdir()
    shutil.copyfile(REFERENCE / 'scores.csv', root / 'scores.csv')
    manifest = REFERENCE / 'run_manifest.json'
    if not manifest.exists():
        manifest = FIXTURES / 'run_manifest_v14.json'
    shutil.copyfile(manifest, root / 'run_manifest.json')
    return root


def read_scores(root):
    return pd.read_csv(root/'scores.csv', float_precision='round_trip', dtype={'phenotype': str})


def write_scores(root, data, update_checksum=True):
    data.to_csv(root/'scores.csv', index=False, float_format='%.17g')
    if update_checksum:
        manifest = json.loads((root/'run_manifest.json').read_text())
        manifest['scores_sha256'] = digest(root/'scores.csv')
        (root/'run_manifest.json').write_text(json.dumps(manifest))


def edit_manifest(root, **values):
    data = json.loads((root/'run_manifest.json').read_text())
    data.update(values)
    (root/'run_manifest.json').write_text(json.dumps(data))


def prior_outputs(root):
    # Real old outputs when present; fallback fixtures retain valid numeric files
    # and minimal stand-ins only for figures absent from the public source tables.
    for name in m.OUTPUT_NAMES:
        source = REFERENCE / name
        if source.exists():
            shutil.copyfile(source, root/name)
        elif name == 'bootstrap_draws.npz':
            shutil.copyfile(FIXTURES/'bootstrap_draws_v14.npz', root/name)
        else:
            (root/name).write_bytes(b'prior figure fixture: '+name.encode())
    return {n: (root/n).read_bytes() for n in m.OUTPUT_NAMES}


def unchanged(root, before):
    assert {n: (root/n).read_bytes() for n in m.OUTPUT_NAMES if (root/n).exists()} == before
    assert not (root/m.LOCK_NAME).exists()
    assert not (root/m.JOURNAL_NAME).exists()
    assert not list(root.glob('.e5-summary-stage-*'))


def test_normal_matches_recorded_v14_table_bootstraps_and_summary(run_root):
    d, counts = m.validate_inputs(run_root)
    table, boot, summary = m.compute_statistics(d, counts)
    path = run_root/'computed.csv'
    table.to_csv(path, index=False, float_format='%.17g')
    assert path.read_bytes() == (REFERENCE/'operating_rates.csv').read_bytes()
    baseline = REFERENCE/'bootstrap_draws.npz'
    if not baseline.exists(): baseline = FIXTURES/'bootstrap_draws_v14.npz'
    with np.load(baseline) as old:
        assert set(boot) == set(old.files)
        for key, values in boot.items(): np.testing.assert_array_equal(values, old[key])
    old = json.loads((REFERENCE/'summary.json').read_text())
    for key in ['reference_realizations','score_rows','model_fits','failures','phenotypes_reported']:
        assert summary[key] == old[key]
    assert summary['status'] == 'COMPLETE'
    assert counts['score_n'] == len(d)
    assert summary['failures'] == len(json.loads((run_root/'run_manifest.json').read_text())['failed'])


@pytest.mark.parametrize('previous', [False, True])
@pytest.mark.parametrize('fault', [
    'row_failed','manifest_failed','manifest_failed_list','cal_nan','cal_inf','r2_inf',
    'missing_group','duplicate_id','wrong_grid','wrong_id','overlap_id','wrong_replicate',
    'wrong_seed','difference_mismatch','checksum','wrong_count','wrong_protocol',
    'missing_manifest','corrupt_manifest','missing_scores'])
def test_bad_input_rejected_before_statistics_or_outputs(run_root, monkeypatch, previous, fault):
    before = prior_outputs(run_root) if previous else {}
    if fault.startswith('manifest_'):
        edit_manifest(run_root, **({'status':'FAILED'} if fault=='manifest_failed' else {'failed':['injected']}))
    elif fault == 'wrong_count': edit_manifest(run_root, completed_model_fits=19999)
    elif fault == 'wrong_protocol': edit_manifest(run_root, protocol_sha256='0'*64)
    elif fault == 'missing_manifest': (run_root/'run_manifest.json').unlink()
    elif fault == 'corrupt_manifest': (run_root/'run_manifest.json').write_text('{broken')
    elif fault == 'missing_scores': (run_root/'scores.csv').unlink()
    else:
        d = read_scores(run_root)
        if fault == 'row_failed': d.loc[0,'status']='FAILED'
        elif fault in ['cal_nan','cal_inf']: d.loc[0,'T_scope']=np.nan if fault=='cal_nan' else np.inf
        elif fault == 'r2_inf': d.loc[0,'R2_B']=-np.inf
        elif fault == 'missing_group': d=d[~((d.role=='evaluation') & (d.phenotype=='4') & (d.L==400))]
        elif fault == 'duplicate_id': d=pd.concat([d,d.iloc[[0]]],ignore_index=True)
        elif fault == 'wrong_grid': d.loc[d.index[d.L==400][0],'L']=401
        elif fault == 'wrong_id': d.loc[0,'reference_id']='calibration_999'
        elif fault == 'overlap_id': d.loc[d.index[d.role=='evaluation'][0],'reference_id']='calibration_000'
        elif fault == 'wrong_replicate': d.loc[0,'replicate']=100
        elif fault == 'wrong_seed': d.loc[0,'fold_seed']+=1
        elif fault == 'difference_mismatch': d.loc[0,'T_scope']+=.01
        elif fault == 'checksum': d.loc[0,'status']='FAILED'
        write_scores(run_root,d,update_checksum=fault!='checksum')
    def should_not_compute(*args):
        raise RuntimeError('Statistics were incorrectly reached with invalid input')
    monkeypatch.setattr(m,'compute_statistics',should_not_compute)
    with pytest.raises(ValueError): m.summarize(run_root)
    unchanged(run_root,before)


def test_manifest_and_score_identity_from_single_byte_snapshot(run_root):
    d, counts=m.validate_inputs(run_root)
    assert counts['scores_sha256']==digest(run_root/'scores.csv')
    assert counts['manifest_sha256']==digest(run_root/'run_manifest.json')
    assert d.reference_id.nunique()==counts['reference_n']


@pytest.mark.parametrize('previous',[False,True])
def test_render_failure_preserves_output_set(run_root,monkeypatch,previous):
    before=prior_outputs(run_root) if previous else {}
    def broken(stage,table,boot,summary):
        (stage/'operating_rates.csv').write_text('new partial staged file')
        raise RuntimeError('injected render failure')
    monkeypatch.setattr(m,'render_outputs',broken)
    with pytest.raises(RuntimeError,match='injected render'):m.summarize(run_root)
    unchanged(run_root,before)


@pytest.mark.parametrize('previous',[False,True])
def test_actual_savefig_failure_preserves_output_set(run_root,monkeypatch,previous):
    before=prior_outputs(run_root) if previous else {}
    savefig=m.plt.Figure.savefig
    calls=[]
    def interrupted(figure,path,*args,**kwargs):
        calls.append(str(path))
        if len(calls)==2:raise RuntimeError('injected second-format render failure')
        return savefig(figure,path,*args,**kwargs)
    monkeypatch.setattr(m.plt.Figure,'savefig',interrupted)
    with pytest.raises(RuntimeError,match='second-format render'):m.summarize(run_root)
    assert len(calls)==2
    unchanged(run_root,before)


def staged_fixture(stage,table,boot,summary):
    # Commit fault tests exercise publication independently of Matplotlib rendering.
    for name in m.OUTPUT_NAMES:(stage/name).write_bytes(b'new staged artifact '+name.encode())


@pytest.mark.parametrize('previous',[False,True])
@pytest.mark.parametrize('fail_at',[1,3,6])
def test_commit_failure_rolls_back_entire_previous_set(run_root,monkeypatch,previous,fail_at):
    before=prior_outputs(run_root) if previous else {}
    monkeypatch.setattr(m,'render_outputs',staged_fixture)
    replace=m.os.replace;calls=[]
    def interrupted(source,dest):
        calls.append((Path(source).name,Path(dest).name))
        if len(calls)==fail_at:raise OSError('injected commit failure')
        return replace(source,dest)
    monkeypatch.setattr(m.os,'replace',interrupted)
    with pytest.raises(RuntimeError,match='previous output set restored'):m.summarize(run_root)
    unchanged(run_root,before)
    assert calls[0][1]==m.OUTPUT_NAMES[0]


def test_success_commits_summary_last_after_rendering(run_root,monkeypatch):
    monkeypatch.setattr(m,'render_outputs',staged_fixture)
    replace=m.os.replace;order=[]
    def observed(source,dest):
        assert (run_root/m.JOURNAL_NAME).exists()
        order.append(Path(dest).name)
        return replace(source,dest)
    monkeypatch.setattr(m.os,'replace',observed)
    m.summarize(run_root)
    assert order==list(m.OUTPUT_NAMES)
    assert order[-1]=='summary.json'
    assert not (run_root/m.JOURNAL_NAME).exists()


def test_rollback_failure_retains_recovery_evidence_and_fails_closed(run_root,monkeypatch):
    prior_outputs(run_root)
    monkeypatch.setattr(m,'render_outputs',staged_fixture)
    def denied(*args):raise PermissionError('injected persistent filesystem failure')
    monkeypatch.setattr(m.os,'replace',denied)
    with pytest.raises(RuntimeError,match='manual recovery'):m.summarize(run_root)
    journal=json.loads((run_root/m.JOURNAL_NAME).read_text())
    assert journal['status']=='ROLLBACK_FAILED'
    assert all((run_root/journal['stage_directory']/'backup'/n).exists() for n in m.OUTPUT_NAMES)
    with pytest.raises(ValueError,match='Unresolved E5 transaction'):m.summarize(run_root)


def test_concurrent_or_stale_writer_lock_rejected(run_root):
    (run_root/m.LOCK_NAME).write_text('existing writer')
    with pytest.raises(ValueError,match='writer lock'):m.summarize(run_root)
    assert not (run_root/'summary.json').exists()


@pytest.mark.parametrize('bad',[(np.nan,.2,.4),(.3,np.nan,.4),(.3,.4,.2),(.3,.31,.4),(.3,-.01,.4),(.3,.2,1.01)])
def test_invalid_intervals_are_never_hidden(bad):
    with pytest.raises(ValueError):m.interval_yerr(*([v] for v in bad))


def test_tiny_roundoff_in_valid_interval_only_is_normalized():
    np.testing.assert_array_equal(m.interval_yerr([0.0],[1e-16],[.05]),[[0.0],[.05]])
    with pytest.raises(ValueError):m.interval_yerr([0.0],[1e-9],[.05])


@pytest.mark.parametrize('optimization',[[],['-O']])
def test_cli_rejects_failed_row_even_under_optimization(run_root,optimization):
    d=read_scores(run_root);d.loc[0,'status']='FAILED';write_scores(run_root,d)
    result=subprocess.run([sys.executable,*optimization,str(SCRIPT),'--root',str(run_root)],capture_output=True,text=True)
    assert result.returncode!=0
    assert 'Every score row status must be COMPLETE' in result.stderr
    unchanged(run_root,{})
