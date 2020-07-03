# SV pipeline

Snakemake pipeline for Delly and CNVnator.
Written to be run on uppmax using the slurm scheduler.

Prior to running the pipeline, specify the reference-files, directory with bam-files,
and singularity containers to be used, must be specified in `Snakefile`.
Slurm parameters must be specified in `cluster.yaml`

Pipeline can be run with

```bash
snakemake -p -n --cluster-config cluster.yaml \
 --cluster-status ./slurm-status.py \
 --use-singularity \
 --cluster "sbatch --parsable -A {cluster.account} -p {cluster.partition} -n {cluster.n} -t {cluster.time} --job-name {cluster.name} --error {cluster.error} --output {cluster.output}" \
 all
```
