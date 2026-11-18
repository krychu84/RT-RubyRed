#!/usr/bin/env bash
# RT-RubyRed installer
#   1. conda environment "qiime2" (created as in the RubyRed README if it does not exist)
#   2. the original RubyRed pipeline from https://github.com/reillyeo/RubyRed, including its
#      reference database (classifier, sequences, taxonomy), which is stored with Git LFS
#   3. torchlex and the Python packages RT-RubyRed needs
#   4. RT-RubyRed on PATH
#
# Usage:  ./install.sh                  install / check
#         ./install.sh --update-rubyred also pull the latest RubyRed (and its database)
# Environment: RUBYRED_HOME (default ~/my_scripts/RubyRed), QIIME_ENV (default qiime2),
#              CONDA_BASE (default: ~/miniforge3, ~/miniconda3 or ~/anaconda3)
set -euo pipefail

RT_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
RUBYRED_HOME="${RUBYRED_HOME:-$HOME/my_scripts/RubyRed}"
ENV_NAME="${QIIME_ENV:-qiime2}"
RUBYRED_URL="https://github.com/reillyeo/RubyRed"
QIIME_YML="https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2024.5/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml"
UPDATE=0
[ "${1:-}" = "--update-rubyred" ] && UPDATE=1

say()  { printf '\033[1;31m▸\033[0m %s\n' "$*"; }
ok()   { printf '  \033[32mok\033[0m  %s\n' "$*"; }
die()  { printf '\033[31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }

# ---- conda
if [ -z "${CONDA_BASE:-}" ]; then
    for d in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3"; do
        [ -f "$d/etc/profile.d/conda.sh" ] && CONDA_BASE="$d" && break
    done
fi
[ -n "${CONDA_BASE:-}" ] || die "conda not found; install Miniforge or set CONDA_BASE"
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"

say "conda environment '$ENV_NAME'"
if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "  not found: creating it as described in the RubyRed README (QIIME 2 2024.5 amplicon; takes a while)"
    conda env create -n "$ENV_NAME" --file "$QIIME_YML"
    conda install -y -n "$ENV_NAME" -c bioconda -c conda-forge chopper seqkit parallel
fi
conda activate "$ENV_NAME"
ok "$ENV_NAME ($(qiime --version 2>/dev/null | head -1 || echo 'qiime not found'))"

# ---- original RubyRed (+ database via Git LFS)
say "RubyRed ($RUBYRED_URL) in $RUBYRED_HOME"
if ! command -v git-lfs >/dev/null 2>&1; then
    echo "  installing git-lfs into $ENV_NAME (the RubyRed database is stored with Git LFS)"
    conda install -y -n "$ENV_NAME" --freeze-installed -c conda-forge git-lfs
fi
if [ ! -d "$RUBYRED_HOME/.git" ]; then
    mkdir -p "$(dirname "$RUBYRED_HOME")"
    git lfs install --skip-repo >/dev/null
    git clone "$RUBYRED_URL" "$RUBYRED_HOME"
elif [ "$UPDATE" = 1 ]; then
    git -C "$RUBYRED_HOME" pull --ff-only
else
    echo "  already installed ($(git -C "$RUBYRED_HOME" log -1 --format='%h %cs')); use --update-rubyred to pull the latest version"
fi
git -C "$RUBYRED_HOME" lfs install --local >/dev/null
git -C "$RUBYRED_HOME" lfs pull
chmod +x "$RUBYRED_HOME/RubyRed" 2>/dev/null || true
ok "RubyRed $(git -C "$RUBYRED_HOME" log -1 --format='%h (%cs)')"

# ---- torchlex and Python packages
say "torchlex and Python packages"
if ! command -v torchlex >/dev/null 2>&1; then
    pip install -e "$RUBYRED_HOME/scripts/torchlex"
fi
python -c "import pandas, psutil" 2>/dev/null || conda install -y -n "$ENV_NAME" --freeze-installed -c conda-forge pandas psutil
for c in qiime torchlex cutadapt chopper vsearch seqkit parallel biom git-lfs; do
    command -v "$c" >/dev/null && ok "$c" || die "$c missing from the $ENV_NAME environment"
done

# ---- database check: the .qza files must be the real files, not LFS pointers
say "RubyRed database"
RUBYRED_HOME="$RUBYRED_HOME" python - "$RT_DIR" <<'PY'
import os, sys
sys.path.insert(0, sys.argv[1])
import rt_rubyred as r
rr = os.path.join(os.environ["RUBYRED_HOME"], "RubyRed")
db = r.rubyred_database(rr)
bad = r.lfs_problems(db.values()) if db else ["no default database found in " + rr]
if bad:
    sys.exit("ERROR: database not usable:\n  " + "\n  ".join(bad))
for k, p in db.items():
    print(f"  \033[32mok\033[0m  {k}: {os.path.relpath(p, os.environ['RUBYRED_HOME'])} ({os.path.getsize(p)/1e6:.1f} MB)")
PY

# ---- RT-RubyRed on PATH
say "RT-RubyRed"
chmod +x "$RT_DIR/RT-RubyRed"
if ! grep -qs "$RT_DIR" "$HOME/.bashrc"; then
    echo "export PATH=\"$RT_DIR:\$PATH\"" >> "$HOME/.bashrc"
    echo "  added $RT_DIR to PATH in ~/.bashrc"
fi
if [ "$RUBYRED_HOME" != "$HOME/my_scripts/RubyRed" ] && ! grep -qs "RUBYRED_HOME=" "$HOME/.bashrc"; then
    echo "export RUBYRED_HOME=\"$RUBYRED_HOME\"" >> "$HOME/.bashrc"
    echo "  added RUBYRED_HOME to ~/.bashrc"
fi
ok "$("$RT_DIR/RT-RubyRed" -v)"
echo
echo "Done. Open a new terminal (or: source ~/.bashrc), then:  RT-RubyRed -w /path/to/fastq_pass/barcodeNN"
