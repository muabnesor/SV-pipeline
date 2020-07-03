from pathlib import Path

__version__ == "0.1.0"

bam_dir_path = "/path/to/bams"
slurm_log_dir = "slurm_log"
scripts_dir = "scripts"

reference_fasta = "/path/to/reference.fa"
gnomad_sv = "/path/to/gnomad_SV.fa"
exclude_tsv = "/path/to/delly.human.hg19.excl.tsv"
vep_cache_dir = "/path/to/vep_cache"

def find_bams(bam_dir):
	bam_dict = dict()
	for sub_file in Path(bam_dir).iterdir():
		if not sub_file.is_file:
			continue
		if not sub_file.suffix == ".bam":
			continue

		sample_name = sub_file.name.split(".")[0]
		bam_dict[sample_name] = sub_file

	return bam_dict

bam_dict = find_bams(bam_dir=bam_dir_path)

def find_bam(wildcards):
	return f"{bam_dict[wildcards.sample]}"

delly_dir_path = "delly_out"
cnvnator_dir_path = "cnvnator_out"

rule delly_call:
	input:
		bam_file = find_bam,
		reference_fasta = f"{reference_fasta}",
		exclude_tsv = f"{exclude_tsv}"

	output:
		bcf_file = f"{delly_dir_path}/{{sample}}.bcf"

	params:
		slurm_log_dir = f"{slurm_log_dir}"

	threads: 1

	singularity:
		f"/path/to/delly.sif"

	shell:
		"delly call -g {input.reference_fasta} -o {output.bcf_file} -x {input.exclude_tsv} {input.bam_file}"


rule delly_merge:
	input:
		bcf_files = expand(f"{delly_dir_path}/{{sample}}.bcf", sample=bam_dict.keys())

	output:
		merged_bcf = f"{delly_dir_path}/sites.bcf"

	params:
		slurm_log_dir = f"{slurm_log_dir}"

	threads: 1

	singularity:
		f"/path/to/delly.sif"

	shell:
		"delly merge -o {output.merged_bcf} -m 50 -p {input.bcf_files}"


rule delly_genotype:
	input:
		reference_fasta = f"{reference_fasta}",
		exclude_tsv = f"{exclude_tsv}",
		sites_bcf = f"{delly_dir_path}/sites.bcf",
		bam_file = find_bam

	output:
		genotyped_bcf = f"{delly_dir_path}/{{sample}}_genotyped.bcf"

	params:
		slurm_log_dir = f"{slurm_log_dir}"

	threads: 1

    singularity:
		f"/path/to/delly.sif"

	shell:
		"delly call -g {input.reference_fasta} -v {input.sites_bcf} -o {output.genotyped_bcf} -x {input.exclude_tsv} {input.bam_file}"


rule bcftools_merge:
	input:
		bcf_files = expand(f"{delly_dir_path}/{{sample}}_genotyped.bcf", sample=bam_dict.keys())

	output:
		merged_bcf = f"{delly_dir_path}/merged.bcf"

	params:
		slurm_log_dir = f"{slurm_log_dir}"

	threads: 1

    singularity:
        f"/path/to/bcftools.sif"

	shell:
		"bcftools merge -m id -O b -o {output.merged_bcf} {input.bcf_files}"

rule bcftools_index:
	input:
		merged_bcf = f"{delly_dir_path}/merged.bcf"

	output:
		bcf_index = f"{delly_dir_path}/merged.bcf.csi"

	params:
		slurm_log_dir = f"{slurm_log_dir}"

	threads: 1

    singularity:
        f"/path/to/bcftools.sif"

	shell:
		"bcftools index {input.merged_bcf}"

rule delly_filter:
	input:
		merged_bcf = f"{delly_dir_path}/merged.bcf",
		bcf_index = f"{delly_dir_path}/merged.bcf.csi"

	output:
		filtered_bcf = f"{delly_dir_path}/filtered.bcf"

	params:
		slurm_log_dir = f"{slurm_log_dir}"

	threads: 1

	singularity:
		f"/path/to/delly.sif"

	shell:
		"delly filter -f germline -o {output.filtered_bcf} {input.merged_bcf}"


rule bcf_to_vcf:
	input:
		filtered_bcf = f"{delly_dir_path}/filtered.bcf"

	output:
		filtered_vcf = f"{delly_dir_path}/filtered.vcf"

	params:
		slurm_log_dir = f"{slurm_log_dir}"

	threads: 1

    singularity:
        f"/path/to/bcftools.sif"

	shell:
		"bcftools view {input.filtered_bcf} > {output.filtered_vcf}"

rule vep_annotate:
	input:
		filtered_vcf = f"{delly_dir_path}/filtered.vcf",
        vep_cache_dir = f"{vep_cache_dir}",
		reference_fasta = f"{reference_fasta}"

	output:
		annotated_vcf = f"{delly_dir_path}/filtered_annotated.vcf"

	params:
		slurm_log_dir = f"{slurm_log_dir}"

	threads: 4

    singularity:
        f"/path/to/vep.sif"

	shell:
		"vep --cache "
		"--offline "
		"--output_file {output.annotated_vcf} "
		"--vcf "
		"--input_file {input.filtered_vcf} "
		"--dir_cache {input.vep_cache_dir}"
		"--species homo_sapiens "
		"--assembly GRCh37 "
		"--everything "
		"--force_overwrite "
		"--per_gene "
		"--fasta {input.reference_fasta} "
		"--fork {threads}"


rule gnomad_delly:
	input:
		in_vcf = f"{delly_dir_path}/filtered_annotated.vcf",
		gnomad_sv_vcf = f"{gnomad_sv}",
		vep_cache_dir = f"{vep_cache_dir}"

	output:
		out_vcf = f"{delly_dir_path}/gnomad_filtered_annotated.vcf"

	params:
		slurm_log_dir = f"{slurm_log_dir}"

	threads: 4

	singularity:
        f"/path/to/vep.sif"

	shell:
		"vep "
		"--offline "
		"--output_file {output.out_vcf} "
		"--vcf "
		"--input_file {input.in_vcf} "
		"--dir_cache {input.vep_cache_dir}/ "
		"--species homo_sapiens "
		"--assembly GRCh37 "
		"--plugin StructuralVariantOverlap,file={input.gnomad_sv_vcf},cols=AF,same_type=1,distance=3500,label=gnomad_sv "
		"--force_overwrite "
		"--fork {threads}"


rule delly_final_filter:
	input:
		in_vcf = f"{delly_dir_path}/gnomad_filtered_annotated.vcf"

	output:
		out_vcf = f"{delly_dir_path}/final_gnomad_filtered_annotated.vcf"

	params:
		slurm_log_dir = f"{slurm_log_dir}",
		script = f"{scripts_dir}/SV-filter.py",
		caller = "delly",
		max_gnomad = 0.01

	threads: 1

	shell:
		"module load bioinfo-tools cyvcf2/0.8.8 && "
		"python {params.script} "
		"{input.in_vcf} "
		"{params.caller} "
		"{params.max_gnomad} "
		"> {output.out_vcf}"


### CNVNATOR

rule cnvnator:
	input:
		bam_file = find_bam,
		reference_fasta = f"{reference_fasta}"

	output:
		root_file = f"{cnvnator_dir_path}/{{sample}}.root",
		tsv_file = f"{cnvnator_dir_path}/{{sample}}.tsv",
		vcf_file = f"{cnvnator_dir_path}/{{sample}}.vcf"

	params:
		bin_size = 1000,
		slurm_log_dir = f"{slurm_log_dir}"

	threads: 2

	singularity:
        f"/path/to/cnvnator.sif"

	shell:
		"cnvnator -root {output.root_file} -tree {input.bam_file} && "
		"cnvnator -root {output.root_file} -his {params.bin_size} -fasta {input.reference_fasta} && "
		"cnvnator -root {output.root_file} -stat {params.bin_size} && "
		"cnvnator -root {output.root_file} -partition {params.bin_size} && "
		"cnvnator -root {output.root_file} -call {params.bin_size} > {output.tsv_file} && "
		"cnvnator2VCF.pl {output.tsv_file} > {output.vcf_file}"


rule filter_cnvnator_calls:
	input:
		vcf_file = f"{cnvnator_dir_path}/{{sample}}.vcf"

	output:
		filtered_vcf = f"{cnvnator_dir_path}/{{sample}}_filtered.vcf"

	params:
		slurm_log_dir = f"{slurm_log_dir}"

	threads: 1

    singularity:
        "/path/to/bcftools.sif"

	shell:
		"bcftools view -i 'INFO/natorQ0 < 0.5 & INFO/natorP1 < 0.5 & INFO/natorP3 < 0.5' {input.vcf_file} > {output.filtered_vcf}"

rule survivor_cnvnator:
	input:
		vcf_files = expand(f"{cnvnator_dir_path}/{{sample}}_filtered.vcf", sample=bam_dict.keys())

	output:
		vcf_files_list = f"cnvnator_vcf_files.txt",
		merged_vcf = f"{cnvnator_dir_path}/cnvnator_merged.vcf"

	params:
		slurm_log_dir = f"{slurm_log_dir}",
		max_dist = 3500,

	threads: 1

    singularity:
        "/path/to/suvivor.sif"

	shell:
		"for vcf_file in {input.vcf_files} ; do echo $vcf_file >> {output.vcf_files_list} ; done && "
		"SURVIVOR merge {output.vcf_files_list} {params.max_dist} 1 1 0 1 1000 {output.merged_vcf} && "
        "rm {output.vcf_files_list}"


rule vep_cnvnator:
	input:
		vcf_file = f"{cnvnator_dir_path}/cnvnator_merged.vcf",
		reference_fasta = f"{reference_fasta}",
		gnomad_sv_vcf = f"{gnomad_sv}",
		vep_cache_dir = f"{vep_cache_dir}"

	output:
		annotated_vcf = f"{cnvnator_dir_path}/annotated_merged_cnvnator.vcf"

	params:
		slurm_log_dir = f"{slurm_log_dir}"

	threads: 4

	singularity:
		f"/Path/to/vep.sif"

	shell:
		"vep --cache "
		"--offline "
		"--output_file {output.annotated_vcf} "
		"--vcf "
		"--no_check_variants_order "
		"--input_file {input.vcf_file} "
		"--dir_cache {input.vep_cache_dir}/ "
		"--species homo_sapiens "
		"--assembly GRCh37 "
		"--everything "
		"--force_overwrite "
		"--per_gene "
		"--fasta {input.reference_fasta} "
		"--fork {threads}"


rule gnomad_cnvnator:
	input:
		in_vcf = f"{cnvnator_dir_path}/annotated_merged_cnvnator.vcf",
		gnomad_sv_vcf = f"{gnomad_sv}",
		vep_cache_dir = f"{vep_cache_dir}"

	output:
		out_vcf = f"{cnvnator_dir_path}/gnomad_annotated_merged_cnvnator.vcf"

	params:
		slurm_log_dir = f"{slurm_log_dir}"

	threads: 4

	singularity:
		f"/path/to/vep.sif"

	shell:
		"vep "
		"--offline "
		"--output_file {output.out_vcf} "
		"--vcf "
		"--no_check_variants_order "
		"--input_file {input.in_vcf} "
		"--dir_cache {input.vep_cache_dir}/ "
		"--species homo_sapiens "
		"--assembly GRCh37 "
		"--plugin StructuralVariantOverlap,file={input.gnomad_sv_vcf},cols=AF,same_type=1,distance=3500,label=gnomad_sv "
		"--force_overwrite "
		"--fork {threads}"


rule cnvnator_final_filter:
	input:
		in_vcf = f"{cnvnator_dir_path}/gnomad_annotated_merged_cnvnator.vcf"

	output:
		out_vcf = f"{cnvnator_dir_path}/final_gnomad_annotated_merged_cnvnator.vcf"

	params:
		slurm_log_dir = f"{slurm_log_dir}",
		script = f"{scripts_dir}/SV-filter.py",
		caller = "cnvnator",
		max_gnomad = 0.01

	threads: 1

	shell:
		"module load bioinfo-tools cyvcf2/0.8.8 && "
		"python {params.script} "
		"{input.in_vcf} "
		"{params.caller} "
		"{params.max_gnomad} "
		"> {output.out_vcf}"

rule all:
	input:
		delly_final = f"{delly_dir_path}/final_gnomad_filtered_annotated.vcf",
		cnvnator_final = f"{cnvnator_dir_path}/final_gnomad_annotated_merged_cnvnator.vcf"
