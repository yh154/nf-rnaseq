#!/usr/bin/env bash
# Every bin/ script (and tests/make_synthetic.py): --help => exit 0 + "Usage"; --version => 0;
# bad/missing arguments => exit 2 + "Usage". main.nf --help is tested in run_tests.sh.
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fail=0
check() {
    local script="$1" interp="$2" out rc
    out=$($interp "$script" --help 2>&1); rc=$?
    if [[ $rc -ne 0 || "$out" != *Usage* ]]; then echo "FAIL --help: $script (rc=$rc)"; fail=1; return; fi
    $interp "$script" --version >/dev/null 2>&1 || { echo "FAIL --version: $script"; fail=1; }
    out=$($interp "$script" --definitely-bad-arg 2>&1); rc=$?
    if [[ $rc -ne 2 || "$out" != *Usage* ]]; then echo "FAIL bad args: $script (rc=$rc)"; fail=1; return; fi
    echo "ok   $(basename "$script")"
}
for f in "$ROOT"/bin/*; do
    case "$f" in
        *.py) check "$f" "python3 -I" ;;
        *.R)  command -v Rscript >/dev/null && check "$f" Rscript || echo "skip $(basename "$f") (Rscript missing)" ;;
        *.sh) check "$f" bash ;;
    esac
done
check "$ROOT/tests/make_synthetic.py" "python3 -I"
exit $fail
