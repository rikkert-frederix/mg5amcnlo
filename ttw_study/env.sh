# Source this file from the checkout root before study commands.
export TTW_ROOT=/export/slow1/rikkert/mg5amcnlo
export PATH="$TTW_ROOT/LHAPDF/bin:$TTW_ROOT/ttw_study/local/fastjet/bin:$PATH"
export LD_LIBRARY_PATH="$TTW_ROOT/LHAPDF/lib:$TTW_ROOT/ttw_study/local/fastjet/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PYTHONPATH="$TTW_ROOT/LHAPDF/local/lib/python3.12/dist-packages:$TTW_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export LHAPDF_DATA_PATH="$TTW_ROOT/LHAPDF/share/LHAPDF"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
