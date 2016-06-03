  sbatch \
    --time 00:01:00  --mem=12g --cpus-per-task=1 \
    -D `pwd` -J "example[57,67-69]" \
    -a 57,67-69 \
    -o example.%A_%a.out \
    example.sh
