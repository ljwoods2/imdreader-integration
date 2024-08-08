#!/bin/bash
#SBATCH -J GMX_IMDREADER
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --gres=gpu:1

echo "hostname: $(hostname)"
echo "starting: $(date)"
echo "SLURM env:"
env | grep ^SLURM | sort
env | grep CUDA

module load gromacs-2023.1-gcc-11.2.0 
module load mamba/latest
source activate imdreader-test

output_dir="output_${SLURM_JOB_ID}"
mkdir -p $output_dir
OUTPUT_FILE="${output_dir}/slurm.out"

await_gmx_imd() {
    grep -q "IMD: Will wait until I have a connection and IMD_GO orders." "$OUTPUT_FILE"
}

echo "Using Conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"
gmx mdrun -s md.tpr -deffnm "${output_dir}/test" -imdwait -maxh 0.03 &> "$OUTPUT_FILE" &

while ! await_gmx_imd; do
    echo "Waiting for GROMACS IMD readiness in $OUTPUT_FILE..."
    sleep 5 
done

echo "GROMACS is ready. Running IMDReader"

python client.py

source deactivate

echo "Finished: $(date)"