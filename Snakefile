# Snakefile for Sparrow Chromosome Z Project
# Focus: 50kb window size analysis
# Usage: snakemake --use-conda --cores 4

configfile: "config.yaml"

# Target output files
rule all:
    input:
        "results/figures/Figure1_Genome_Dxy.pdf",
        "results/figures/Figure2_ChrZ_Zoom_50kb.pdf"

# Pre-processing the VCF: removing outgroup and keeping only biallelic SNPs
rule subset_biallelic:
    input:
        config["raw_vcf"]
    output:
        temp("intermediate/ProjTaxa.noOutgroup.biallelic.vcf.gz")
    params:
        # Exclude outgroup and keep only biallelic SNPs
        extra = f"-s ^{config['outgroup']} -m2 -M2 -v snps"
    wrapper:
        "v9.0.0/bio/bcftools/view"

# Apply hard filters (missing data, depth, quality, etc)
rule filter_vcf:
    input:
        "intermediate/ProjTaxa.noOutgroup.biallelic.vcf.gz"
    output:
        temp("results/filtered_vcf/ProjTaxa.filtered.final.vcf")
    params:
        extra = lambda wildcards: (
            f"--max-missing {config['max_missing']} "
            f"--min-meanDP {config['min_mean_dp']} "
            f"--max-meanDP {config['max_mean_dp']} "
            f"--minDP {config['min_dp']} "
            f"--maxDP {config['max_dp']} "
            f"--minGQ {config['min_gq']} "
            f"--maf {config['min_maf']} "
            "--recode-INFO-all"
        )
    wrapper:
        "v9.0.0/bio/vcftools/filter"

# Accession requires bgzipped VCF for Pixy
rule compress_vcf:
    input:
        "results/filtered_vcf/ProjTaxa.filtered.final.vcf"
    output:
        "results/filtered_vcf/ProjTaxa.filtered.final.vcf.gz"
    params:
        extra=""
    wrapper:
        "v3.13.4/bio/bgzip"

# Indexing VCF for random access
rule index_vcf:
    input:
        "results/filtered_vcf/ProjTaxa.filtered.final.vcf.gz"
    output:
        "results/filtered_vcf/ProjTaxa.filtered.final.vcf.gz.tbi"
    wrapper:
        "v8.1.1/bio/tabix/index"

# Calculate Dxy (Absolute Divergence) in 50kb sliding windows
rule pixy_dxy:
    input:
        vcf = "results/filtered_vcf/ProjTaxa.filtered.final.vcf.gz",
        index = "results/filtered_vcf/ProjTaxa.filtered.final.vcf.gz.tbi"
    output:
        "results/pixy/ProjTaxa_dxy.txt"
    conda:
        "envs/pixy.yaml"
    params:
        pop_file = config["pixy_pop_file"],
        window_size = config["pixy_window_size"],
        out_dir = "results/pixy",
        stats = "dxy"
    threads: 4
    shell:
        """
        pixy --stats {params.stats} \
             --vcf {input.vcf} \
             --populations {params.pop_file} \
             --window_size {params.window_size} \
             --output_folder {params.out_dir} \
             --output_prefix ProjTaxa \
             --n_cores {threads} \
             --bypass_invariant_check
        """

# Generate final plots and extract candidate genes
rule final_report:
    input:
        dxy = "results/pixy/ProjTaxa_dxy.txt",
        script = "final_analysis.R"
    output:
        "results/figures/Figure1_Genome_Dxy.pdf",
        "results/figures/Figure2_ChrZ_Zoom_50kb.pdf"
    shell:
        "Rscript {input.script}"
