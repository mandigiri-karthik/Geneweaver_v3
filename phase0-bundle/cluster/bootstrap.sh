#!/usr/bin/env bash
# Environment setup for the IDA GPU pod (macavaney/idabox, OpenShift scc=anyuid).
#
# Creates a DEDICATED conda env for the Phase 0 pipeline, separate from the Bakta env.
# That separation is deliberate: Bakta depends on diamond, pyrodigal, pyhmmer and biopython
# itself, and its solver pins versions of its own choosing. Phase 0's numbers are tied to
# DIAMOND 2.2.8 / pyrodigal 3.7.1 / pyhmmer 0.12.3 / biopython 1.88, so sharing one env would
# silently change the tool versions the baseline was measured with.
#
# Everything lives on /mnt/primary: the other PVCs are 95-100% full.
#
#   bash cluster/bootstrap.sh
set -euo pipefail

PRIMARY="${PRIMARY:-/mnt/primary}"
ENV_PREFIX="${ENV_PREFIX:-$PRIMARY/envs/geneweaver}"
DIAMOND_VERSION="${DIAMOND_VERSION:-2.2.8}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Keep conda's package cache off the root filesystem.
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-$PRIMARY/conda-pkgs}"
mkdir -p "$CONDA_PKGS_DIRS" "$PRIMARY/envs"

echo "== disk =="
df -h "$PRIMARY" | tail -1

echo
echo "== conda env: $ENV_PREFIX =="
if [ -d "$ENV_PREFIX" ]; then
    echo "   exists - reusing (delete it to rebuild: conda env remove -p $ENV_PREFIX)"
else
    # Python 3.12: the pod image is py310, and the project requires >=3.12.
    conda create -y -p "$ENV_PREFIX" -c conda-forge python=3.12 pip
fi

echo
echo "== pinned Python packages =="
conda run -p "$ENV_PREFIX" python -m pip install --quiet --upgrade pip
conda run -p "$ENV_PREFIX" python -m pip install --quiet -r "$HERE/../requirements.lock.txt"

echo
echo "== DIAMOND $DIAMOND_VERSION (static binary, into the env's bin) =="
DBIN="$ENV_PREFIX/bin/diamond"
if [ -x "$DBIN" ] && "$DBIN" --version 2>/dev/null | grep -q "$DIAMOND_VERSION"; then
    echo "   already present: $("$DBIN" --version)"
else
    tmp="$(mktemp -d -p "$PRIMARY")"
    curl -fSL --retry 3 -o "$tmp/diamond.tar.gz" \
        "https://github.com/bbuchfink/diamond/releases/download/v${DIAMOND_VERSION}/diamond-linux64.tar.gz"
    tar -xzf "$tmp/diamond.tar.gz" -C "$tmp" diamond
    install -m 0755 "$tmp/diamond" "$DBIN"
    rm -rf "$tmp"
    echo "   installed: $("$DBIN" --version)"
fi

echo
echo "== versions in $ENV_PREFIX =="
conda run -p "$ENV_PREFIX" python -c \
    "import pyrodigal,pyhmmer,Bio;print('pyrodigal',pyrodigal.__version__,'| pyhmmer',pyhmmer.__version__,'| biopython',Bio.__version__)"
"$DBIN" --version

echo
echo "== register the env as a Jupyter kernel =="
conda run -p "$ENV_PREFIX" python -m pip install --quiet ipykernel
# --user puts it in ~/.local/share/jupyter/kernels, which every Jupyter server
# scans. Installing with --prefix into another env's tree gives a kernelspec the
# running server may never find.
conda run -p "$ENV_PREFIX" python -m ipykernel install --user \
    --name geneweaver --display-name "Python 3.12 (geneweaver)" >/dev/null
echo "   kernel installed; servers visible to this shell now list:"
if command -v jupyter >/dev/null; then
    jupyter kernelspec list 2>/dev/null | sed 's/^/     /'
    jupyter kernelspec list 2>/dev/null | grep -q geneweaver \
        && echo "   OK - 'geneweaver' is registered" \
        || echo "   WARNING - 'geneweaver' not listed; see the paths above"
else
    echo "     (jupyter not on PATH in this shell - check in JupyterLab)"
fi

cat <<EOF

Done.

  Activate in a terminal:   conda activate $ENV_PREFIX
  In JupyterLab:            Kernel > Change Kernel > "Python 3.12 (geneweaver)"
  Bakta env (separate):     conda activate $PRIMARY/envs/bakta

Next: edit cluster/paths.json, then run cluster/check_env.py (or 00_setup.ipynb).
EOF
