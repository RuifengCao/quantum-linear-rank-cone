#!/usr/bin/env python3
"""run_job: uniform runner for the jobs registered in jobs.json.

  python3 scripts/run_job.py --list
  python3 scripts/run_job.py run <job-id> [--smoke] [--tag TAG] [--no-stage] [--no-resume] [-- EXTRA ARGS]
  e.g. python3 scripts/run_job.py run a1-campaign -- --budget 7200 --workers 12

What it does (the hand-off protocol between the analysis sandbox and a server):
  1. preflight: python deps, gcc / lrs / mplrs as declared by the job
  2. resume: if the job declares a resume file and the latest results/*/ folder
     holds one, it is copied to the repository root (never overwrites a newer local file)
  3. run: the command (or its --smoke variant) with the log captured to
     <job>.log; the exit code is captured directly (no `tee` pitfall)
  4. JOB_STATUS.json: host, nproc, start/end, elapsed, exit code, mode
  5. stage: on success the declared outputs + log + status are staged into
     results/<date>_<tag>/ via s13_results_commit.py; on failure only the log
     and status are staged under <tag>-failed so the failure is committable too
Tier S jobs are meant for the sandbox, tier X jobs for a server; smoke
variants of X jobs must pass in the sandbox before hand-off.
"""
import argparse, datetime, glob, json, os, shutil, socket, subprocess, sys, time

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
os.chdir(ROOT)
REG = json.load(open('jobs.json'))['jobs']

def have(tool):
    return shutil.which(tool) is not None

def latest_results_file(name):
    cands = sorted(glob.glob(os.path.join('results', '*', name)))
    return cands[-1] if cands else None

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--list', action='store_true')
ap.add_argument('action', nargs='?', choices=['run'])
ap.add_argument('job', nargs='?')
ap.add_argument('--smoke', action='store_true')
ap.add_argument('--tag', default=None)
ap.add_argument('--no-stage', action='store_true')
ap.add_argument('--no-resume', action='store_true')
a, unknown = ap.parse_known_args()   # unknown options (e.g. -- --budget 7200) are appended to the command
extra = ' '.join(x for x in unknown if x != '--')

if a.list or not a.job:
    for k, j in REG.items():
        print(f"{k:22s} tier {j['tier']}  {j.get('runtime',''):40s} smoke={'yes' if 'smoke' in j else 'no'}")
    sys.exit(0)

j = REG[a.job]
mode = 'smoke' if a.smoke else 'full'
if a.smoke and 'smoke' not in j:
    raise SystemExit(f'{a.job} has no smoke variant')
cmd = j['smoke'] if a.smoke else (j['cmd'] + (' ' + extra if extra else ''))

# 1. preflight
missing = []
try:
    import numpy, scipy  # noqa: F401
except Exception:
    missing.append('python: pip install -r requirements.txt')
for tool in j.get('needs', []):
    if not have(tool):
        missing.append(tool)
if missing and not a.smoke:
    raise SystemExit('preflight failed, missing: ' + ', '.join(missing))
if missing:
    print('preflight (smoke): missing ' + ', '.join(missing) + ' -- continuing where possible')

# 2. resume
if j.get('resume') and not a.smoke and not a.no_resume:
    src = latest_results_file(j['resume'])
    dst = j['resume']
    if src and (not os.path.exists(dst) or os.path.getmtime(src) > os.path.getmtime(dst)):
        shutil.copy2(src, dst); print(f'resume: copied {src} -> {dst}')
    elif os.path.exists(dst):
        print(f'resume: using local {dst}')

# 3. run
log = f'{a.job}{"-smoke" if a.smoke else ""}.log'
t0 = time.time(); start = datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')
print(f'== run_job {a.job} [{mode}] on {socket.gethostname()} ({os.cpu_count()} cpus) ==\n$ {cmd}', flush=True)
interrupted = False
with open(log, 'w') as lf:
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    try:
        for line in p.stdout:
            sys.stdout.write(line); lf.write(line)
    except KeyboardInterrupt:
        interrupted = True
        lf.write('\n[run_job] interrupted by user; waiting for the child to save its state\n')
    rc = p.wait()
    if interrupted and rc == 0:
        rc = 130
elapsed = time.time() - t0

# 4. status
status = {'job': a.job, 'mode': mode, 'host': socket.gethostname(), 'nproc': os.cpu_count(),
          'start_utc': start, 'end_utc': datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
          'elapsed_s': round(elapsed, 1), 'exit_code': rc, 'cmd': cmd, 'checks': j.get('checks', '')}
json.dump(status, open('JOB_STATUS.json', 'w'), indent=1)
print(f'== exit {rc} after {elapsed:.0f}s; checks to eyeball: {j.get("checks","")} ==')

# 5. stage
if not a.no_stage and not a.smoke:
    tag = a.tag or j.get('tag', a.job)
    files = [f for f in j.get('outputs', []) if os.path.isfile(f)]   # partial outputs are worth keeping
    files += [log, 'JOB_STATUS.json']
    if interrupted:
        tag += '-interrupted'
    elif rc != 0:
        tag += '-failed'
    subprocess.run([sys.executable, 'scripts/s13_results_commit.py', '--tag', tag,
                    '--note', f'run_job {a.job} exit {rc} on {socket.gethostname()}'] + files, check=False)
sys.exit(rc)
