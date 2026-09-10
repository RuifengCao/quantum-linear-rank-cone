#!/usr/bin/env python3
"""run_session: ONE server session that runs every pending tier-X job in order,
stages each result folder, and writes a summary -- so a rented server is used
once, end to end, without round trips.

  python3 scripts/run_session.py --hours 8            # full session (default plan below)
  python3 scripts/run_session.py --smoke               # sandbox rehearsal (minutes)

Plan (full mode):
  1. selftest --full                     (~1 min; aborts the session if not green)
  2. a1-campaign                          (budget = hours - 1.5 h; --workers = cpus-1; --max-orbits 1.5M)
  3. s7-pools                             (~1 h on 25 cores: 50 M GF(2) samples + full F_3 layer)
  4. SESSION_SUMMARY.md staged under results/<date>_session/
Each step is executed through run_job.py, so every step is staged on its own
(also on failure / interruption) and the session continues past a failed step
where that makes sense.  Commit results/ once at the end.
"""
import argparse, datetime, json, os, re, subprocess, sys, time

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
os.chdir(ROOT)
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--hours', type=float, default=8.0, help='total wall-clock budget of the session')
ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 1))
ap.add_argument('--max-orbits', type=int, default=1500000)
ap.add_argument('--skip', nargs='*', default=[], help='job ids to skip, e.g. --skip s7-pools')
ap.add_argument('--smoke', action='store_true')
a = ap.parse_args()

def run(job, extra=None, smoke=False):
    cmd = [sys.executable, 'scripts/run_job.py', 'run', job] + (['--smoke'] if smoke else [])
    if extra and not smoke:
        cmd += ['--'] + extra
    t0 = time.time()
    print(f'\n##### session step: {" ".join(cmd[2:])}', flush=True)
    rc = subprocess.call(cmd)
    return rc, time.time() - t0

def last_match(path, pattern):
    try:
        hits = re.findall(pattern, open(path, errors='replace').read())
        return hits[-1] if hits else ''
    except Exception:
        return ''

steps = []
t_session = time.time()
# 1. self-test (fast in smoke mode, full otherwise)
rc, el = run('selftest', smoke=False) if not a.smoke else (subprocess.call([sys.executable, 'selftest.py']), 0)
steps.append(('selftest', rc, el, 'ALL PASS' if rc == 0 else 'FAILED'))
if rc != 0:
    print('self-test failed -- session aborted (fix the environment first)')
    sys.exit(rc)
# 2. campaign
if 'a1-campaign' not in a.skip:
    budget = int(max(600, (a.hours - 1.5) * 3600))
    rc, el = run('a1-campaign', ['--budget', str(budget), '--workers', str(a.workers),
                                 '--max-orbits', str(a.max_orbits)], smoke=a.smoke)
    log = 'a1-campaign-smoke.log' if a.smoke else 'a1-campaign.log'
    steps.append(('a1-campaign', rc, el, last_match(log, r'STATE: .*')))
# 3. pools
if 's7-pools' not in a.skip:
    rc, el = run('s7-pools', smoke=a.smoke)
    log = 's7-pools-smoke.log' if a.smoke else 's7-pools.log'
    steps.append(('s7-pools', rc, el, ' | '.join(re.findall(r'\[[AB]\][^\n]*', open(log, errors='replace').read())[-2:]) if os.path.exists(log) else ''))
# 4. summary
day = datetime.date.today().isoformat()
lines = [f'# Session summary ({day}, {"smoke" if a.smoke else "full"} mode)', '',
         f'host: {os.uname().nodename}, cpus: {os.cpu_count()}, workers: {a.workers}, hours: {a.hours}', '',
         '| step | exit | elapsed | key line |', '| --- | :-: | --- | --- |']
for name, rc, el, key in steps:
    lines.append(f'| {name} | {rc} | {el/60:.1f} min | {key} |')
lines += ['', f'total elapsed: {(time.time()-t_session)/3600:.2f} h', '',
          'Next: commit `results/` (all folders created by this session) and push.']
open('SESSION_SUMMARY.md', 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
if not a.smoke:
    subprocess.call([sys.executable, 'scripts/s13_results_commit.py', '--tag', 'session',
                     '--note', f'{len(steps)} steps, {a.hours} h budget', 'SESSION_SUMMARY.md'])
    print(f'\nsuggested commit message:  results: server session {day}')
sys.exit(max((rc for _, rc, _, _ in steps), default=0))
