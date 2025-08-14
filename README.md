# likely need correct version of Java:
module load Java

rm -r results/; rm -r work/; rm -r .singularity/; rm -r .nextflow*

nextflow run main.nf -profile slurm,singularity -params-file config.yaml

