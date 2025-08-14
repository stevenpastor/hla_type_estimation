# HLA Type Estimation from RNA or DNA

## BACKGROUND

* Estimation of HLA haplotypes from RNA or DNA sequencing fastq files

* Identifies Class I, II, and non-canonical types using 3 algorithms: Optitype, PHLAT, and ArcasHLA

* NOTE: ArcasHLA limited to RNA-seq

## INSTALLATION

* Only three requirements: Nextflow and a Linux environment with Singularity

* Install Nextflow, go here: https://www.nextflow.io/docs/latest/install.html

* Make sure you pay attention to the Java instructions in the above link

* Taken from their website: "Nextflow requires Bash 3.2 (or later) and Java 17 (or later, up to 24) to be installed. To see which version of Java you have, run the following command:"

```
java --version
```

## SETUP

* Check the config.yaml file for an example setup

* Pipeline is compatible with paired-end (PE) or single-end (SE) DNA or RNA fastq files

* You need to have a PHLAT index (CHOP users ask me)

* Current config.yaml cpus and memory likely fine for most use cases but increase it if needed (e.g., 100+ samples)

## RUN

* After editing config.yaml:

```
nextflow run main.nf -profile slurm,singularity -params-file config.yaml
```

* Outputs will be in results, with one dir per algorithm

* Optitype: *_result.tsv file

* PHLAT: *.sum file

* ArcasHLA: *.genotype.json file
