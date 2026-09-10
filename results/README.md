# results/

Committed outputs of server sessions and analysis rounds, one folder per
session: `results/<YYYY-MM-DD>_<tag>/` (e.g. `2026-09-10_a1-batch2`).

Each folder is created by `scripts/s13_results_commit.py`, which copies the
files, writes `manifest.json` (sha256 of every file; canonical row-family sha
for 2-D integer `.npy` files) and prints the suggested commit message.
Files larger than 50 MB (e.g. the s7 pools) do not go here: attach them to a
GitHub Release instead and record the release tag in the manifest note.

Round-trip protocol with GitHub Desktop:
1. Run the job on the server (or in the container).
2. `python3 scripts/s13_results_commit.py --tag <tag> <files...>` → creates
   the folder; copy it to the local clone if the job ran elsewhere.
3. GitHub Desktop shows the new folder as changes → commit with the printed
   message → Push.
4. The analysis side pulls `main` and reads `results/` directly.
