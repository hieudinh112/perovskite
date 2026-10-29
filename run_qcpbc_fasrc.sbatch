#!/bin/bash
#SBATCH -p test
#SBATCH -c 32
#SBATCH --mem=900G
#SBATCH -t 12:00:00
#SBATCH -J pvk_hse
#SBATCH -o /n/netscratch/joonholee_lab/Lab/hieudinh/perovskite/slurm-%j.out

# One qcpbc input on 32 threads, with the run dir, scratch and Slurm log on netscratch (copied from
# ~/Fe_Ni/pbe_soc/run_qcpbc_fasrc.sbatch). The binary is the isolated HEAD build a3459caa2 (its RPATH points at
# its own tree, so rebuilds of ~/qcpbc/build cannot touch a running job).
#
#   sbatch run_qcpbc_fasrc.sbatch P-1_100K_hse06_thcaok_svp_k111.in
#
# The input is a file name in ~/perovskite or an absolute path. GUESS_FROM=<scratch folder of an earlier run of
# the same cell and mesh> copies its saved density (dm, nocc*) into the new scratch folder, for scf_guess = read.
# QCPBC_EXE=<qcpbc.exe> overrides the binary. RUN_LABEL=<label> appends _<label> to <name> below.
# For input <name>.in and job <id>:
#   scratch    /n/netscratch/joonholee_lab/Lab/hieudinh/perovskite/scr/<name>.<id>
#   run dir    /n/netscratch/joonholee_lab/Lab/hieudinh/perovskite/<name>.<id>    input copy, output, memory logs
#   results    ~/perovskite/outputs/<name>.<id>.out    copy of the output at the end
#   points     .../perovskite/thc/isdf_pts.<name>.<id>.txt    ISDF points of a save_isdf_points run (they depend
#              on the cell, basis and THC grid only, so later runs can load_isdf_points them)

INP=${1:?usage: run_qcpbc_fasrc.sbatch <input.in>}
SRC=/n/home02/hieudinh/perovskite
case "$INP" in /*) INPFILE=$INP ;; *) INPFILE=$SRC/$INP ;; esac
OUTDIR=$SRC/outputs
NAME=$(basename "$INPFILE" .in)${RUN_LABEL:+_$RUN_LABEL}
ID=${SLURM_JOB_ID:?run this through sbatch or srun}
BASE=/n/netscratch/joonholee_lab/Lab/hieudinh/perovskite
EXE=${QCPBC_EXE:-/n/netscratch/joonholee_lab/Everyone/hieudinh/qcpbc_build_head_a3459caa2/libpbc/qcpbc/qcpbc.exe}

source /etc/profile > /dev/null 2>&1
module purge
module load gmp/6.3.0-fasrc01 mpc/1.3.1-fasrc01 mpfr/4.2.1-fasrc01 intel/24.2.1-fasrc01 intelmpi/2021.13-fasrc01
export LD_LIBRARY_PATH=/n/sw/helmod-rocky8/apps/Core/gcc/14.2.0-fasrc01/lib64:/n/sw/helmod-rocky8/apps/Core/mpfr/4.2.1-fasrc01/lib64:$LD_LIBRARY_PATH
export QCAUX=/n/home02/hieudinh/qcpbc/qcaux
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-32} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-32} OMP_STACKSIZE=512M

SCR=$BASE/scr/$NAME.$ID
RUN=$BASE/$NAME.$ID
mkdir -p "$SCR" "$RUN" "$OUTDIR" "$BASE/thc" || exit 1
cp "$INPFILE" "$RUN/$NAME.in" || exit 1
if [ -n "$GUESS_FROM" ]; then
    cp -p "$GUESS_FROM"/* "$SCR"/ 2> /dev/null     # files only: the autosad/frgmt_* folders are not needed
    ls "$SCR"/dm* > /dev/null || { echo "no density matrix in $GUESS_FROM"; exit 1; }
fi
cd "$RUN" || exit 1

echo "node=$(hostname) start=$(date) threads=$OMP_NUM_THREADS input=$INPFILE"
echo "run dir=$RUN  scratch=$SCR${GUESS_FROM:+  guess from=$GUESS_FROM}"
echo "exe=$EXE md5=$(md5sum "$EXE" | cut -c1-12)"

# Peak memory: /usr/bin/time -v (exact, at exit) and a 10 s sampler (survives a time-limit kill).
( while true; do
    rss=$(ps -u "$USER" -o rss=,args= | awk -v n="$NAME.in" '$0 ~ /qcpbc.exe/ && index($0, n) {s+=$1} END{print s+0}')
    echo "$(date +%s) rss_kb=$rss"
    sleep 10
  done ) > mem_samples.txt &
SAMPLER=$!

/usr/bin/time -v "$EXE" "$NAME.in" "$SCR" > "$NAME.out" 2> time_v.txt
rc=$?
kill $SAMPLER

{
  echo "qcpbc exit code: $rc   end: $(date)"
  grep -E "Maximum resident set size|Elapsed \(wall clock\)|Percent of CPU" time_v.txt
  echo "sampled peak qcpbc.exe RSS (KB): $(sed -E 's/.*rss_kb=([0-9]+).*/\1/' mem_samples.txt | sort -n | tail -1)"
  echo "scratch used: $(du -sh "$SCR" | cut -f1)  ($SCR)"
} > mem_summary.txt
cat mem_summary.txt

[ -f "$SCR/isdf_pts.txt" ] && cp "$SCR/isdf_pts.txt" "$BASE/thc/isdf_pts.$NAME.$ID.txt"
cp "$NAME.out" "$OUTDIR/$NAME.$ID.out"
exit $rc
