#!/bin/sh
# One-command server session (R10): s7 -> s8(#19 certificate) -> s10 -> s11.
# Run from the kit/ directory.  ~25 minutes total on 25 vCPU.
# s8 certificate policy: try the UNCONDITIONAL pool first (GF(2) 12-vertex);
# only if that fails, merge in the F_3 layer -- which is CONDITIONAL (see
# README), so a certificate from the merged pool carries that premise and
# is labelled accordingly in the log.
set -e
python3 -c "import numpy, scipy" 2>/dev/null || {
  echo "Missing dependencies: run  pip install -r requirements.txt  first"; exit 1; }
command -v gcc >/dev/null || { echo "Missing gcc: apt install build-essential"; exit 1; }
step_ok() { [ -e "$1" ] || { echo "Step failed: expected output $1 not found (see the corresponding .log)"; exit 1; }; }
H=QLR_H_pure28.npy
[ -f "$H" ] || python3 scripts/s1_build_qlr.py --variant pure28 --workers "$(nproc)"
echo "== s7 expand (full) =="
python3 scripts/s7_expand_realizable.py --H "$H" 2>&1 | tee s7.log
step_ok s7_out/realizable_gf2_12v.npy
echo "== s8 decompose #19 (unconditional GF(2) pool) =="
python3 scripts/s8_decompose.py --H "$H" \
    --pool s7_out/realizable_gf2_12v.npy 2>&1 | tee s8.log
if grep -q "no certificate" s8.log; then
  echo "== s8 retry: merged pool incl. CONDITIONAL F_3 layer =="
  python3 scripts/s8_decompose.py --H "$H" \
      --pool s7_out/realizable_gf2_12v.npy s7_out/realizable_f3_conditional.npy \
      --out s8_certs_conditional.npz 2>&1 | tee -a s8.log
  echo "NOTE: any certificate above rests on the F_3 conditional premise (README)." | tee -a s8.log
fi
echo "== s10 facet classes (84 LPs) =="
python3 scripts/s10_facets.py --H "$H" 2>&1 | tee s10.log
step_ok pure_facet_classes.npz
echo "== s11 collect =="
python3 scripts/s11_collect.py "$H" 's7_out/*.npy' 's8_certs*.npz' \
    'pure_facet_classes.npz' 's7.log' 's8.log' 's10.log'
echo "DONE -- upload the s11 tar.gz back to the conversation."
