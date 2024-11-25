# imports
import matplotlib.pyplot as plt
import pandas as pd
import re

# global variables
sample_names, = glob_wildcards("/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz")
db = "/lustre1/project/stg_00079/teaching/hg38_21/chr21.fa"
genome_version = "hg38"
snpeff_jar = "/lustre1/project/stg_00079/teaching/I0U19a_conda_2024/share/snpeff-5.2-0/snpEff.jar"

# searches directory for files that end in 7 and saves file names to list
pattern = re.compile(r'^HG\d+7\.GRCh38DH\.exome\.chr21$')
filtered_list = [sample for sample in sample_names if pattern.match(sample)]
print(filtered_list)

# inputs needed to run full workflow
rule all:
    input: 
        fq = expand("fastqc_output/{sample}_fastqc.zip", sample=filtered_list),
        zip = expand("fastqc_output/{sample}_fastqc.zip", sample=filtered_list),
        bam = expand("bwa_output/{sample}.bam", sample=filtered_list),
        bam_bai = expand("bwa_output/{sample}.bam.bai", sample=filtered_list),
        vcf = "snpcall_output/raw_snps.vcf",
        clean = "snpcall_output/clean_snps.vcf",
        annotated = "snpcall_output/annotated_snps.vcf",
        extracted_snps = "genes.vcf",
        rep1 = expand("fastqc_output/{sample}_fastqc/Images/per_base_quality.png", sample = filtered_list),
        stats="barplot_stats/stats.txt",
        snp_counts="barplot_stats/snp_counts.txt",
        rep4 = report("barplot.png", category='SNPs', subcategory="SNPs per individual")

# copy zipped fastq files that end with 7 and save in fastq_files directory
rule copy_files:
    input: 
        fq_gz = expand("/staging/leuven/stg_00079/teaching/1000genomes/{sample}.fq.gz", sample=filtered_list)
    output: 
        dir = directory("fastq_files"),
        cp_fq = expand("fastq_files/{sample}.fq.gz", sample=filtered_list),
    params:
        stdout = "logs/stdout.txt",
        stderr = "logs/stderr.txt"

    shell: 
        """
        if [ ! -d "logs" ]; then
            mkdir logs
        fi
        
        for f in {input.fq_gz}; do
            cp "$f" "{output.dir}"
            if [[ $(sha256sum "$f" | awk '{{print $1}}') != $(sha256sum "{output.dir}/$(basename "$f")" | awk '{{print $1}}') ]]; then
                echo "$f not copied" >> {params.stderr}
            else
                echo "$f copied succesfully" >> {params.stdout}
            fi
        done
        """


# unzip fastq files and save in fastq_unzip directory
rule unzip:
    input: 
        fq_gz = expand("fastq_files/{sample}.fq.gz", sample=filtered_list)
    output: 
        fq_unzip = expand("fastq_unzip/{sample}.fastq", sample=filtered_list)
    params:
        stdout = "logs/stdout.txt",
        stderr = "logs/stderr.txt"    

    shell: 
        """
        for f in {input.fq_gz}; do
            gunzip -c $f > fastq_unzip/$(basename $f .fq.gz).fastq;

            if [[ $(basename "$f" .fq.gz) != $(basename "{output}/$(basename "$f" .fq.gz)") ]]; then
                echo "$f not unzipped" >> {params.stderr}
            else
                echo "$f unzipped succesfully" >> {params.stdout}
            fi
        done
        """

# run fastqc on fastq files and save output in fastqc_output directory
rule fastqc:
    input: 
        fq = "fastq_unzip/{sample}.fastq"
    output: 
        html = "fastqc_output/{sample}_fastqc.html",
        zip = "fastqc_output/{sample}_fastqc.zip",
        summarydata = "fastqc_output/{sample}_fastqc/fastqc_data.txt",
        rep1 = report("fastqc_output/{sample}_fastqc/Images/per_base_quality.png", category="Fastqc", subcategory="Per base quality"),
        rep2 = report("fastqc_output/{sample}_fastqc/Images/per_base_sequence_content.png", category="Fastqc", subcategory="Per base sequence content"),
        rep3 = report("fastqc_output/{sample}_fastqc/summary.txt", category="Fastqc", subcategory="Summary text")
    params:
        stdout = "logs/stdout.txt",
        stderr = "logs/stderr.txt"
    shell: 
        """
        fastqc -o fastqc_output {input.fq} --extract

        if grep FAIL {output.rep3}; then
            echo 'fastqc FAIL for {input}' >> {params.stderr}
        else
            echo 'fastqc PASS for {input}' >> {params.stdout}
        fi
        """

# map sequences onto hg38 chromosome 21 reference genome for chromosome 21 and save output in bwa_output directory
rule bwa:
    input: fq = "fastq_unzip/{sample}.fastq"
    output: 
        bam = "bwa_output/{sample}.bam",
        bam_bai = "bwa_output/{sample}.bam.bai"
    params:
        stdout = "logs/stdout.txt",
        stderr = "logs/stderr.txt"

    shell: 
        """
        bwa mem {db} {input.fq} | samtools sort -o {output.bam} && samtools index {output.bam}

        for f in {output.bam}; do
            mapped_percentage=$(samtools flagstat "$f" | grep -oP '\d+\.\d+%' | head -n1)
                
            if [[ "$mapped_percentage" < 99.0 ]]; then
                echo "WARNING: mapped_percentage for $f is $mapped_percentage" >> {params.stderr}
            else
                echo "mapped_percentage for $f is $mapped_percentage" >> {params.stdout}
            fi
        done
        """

# call variants with bcftools and save raw_snps.vcf in snpcall_output directory
rule snpcall:
    input:
        bam = expand("bwa_output/{sample}.bam", sample = filtered_list)

    output: vcf = "snpcall_output/raw_snps.vcf"

    params:
        stdout = "logs/stdout.txt",
        stderr = "logs/stderr.txt"

    shell: 
        """
        bcftools mpileup -Ou -f {db} {input.bam} | bcftools call -mv -Ov -o {output.vcf}

        header_length=$(grep -c '^#' {output.vcf})

        if [ $header_length -lt 28 ]; then
            echo "WARNING: raw_snps.vcf looks too short" >> {params.stderr}
        else
            echo "raw_snps.vcf created successfully" >> {params.stdout}
        fi
        """

# normalize and filter variant calls and save clean_snps.vcf in snpcall_output directory 
rule cleanup:
    input: vcf = "snpcall_output/raw_snps.vcf",

    output: clean = "snpcall_output/clean_snps.vcf",

    params:
        stdout = "logs/stdout.txt",
        stderr = "logs/stderr.txt"

    shell: 
        """
        cat {input.vcf} | vt decompose - | vt normalize -n -r {db} - | vt uniq - | vt view -f \"QUAL>20\" -h - > {output.clean}

        header_length=$(grep -c '^#' {output.clean})

        if [ $header_length -lt 30 ]; then
            echo "WARNING: clean_snps.vcf looks too short" >> {params.stderr}
        else
            echo "clean_snps.vcf created successfully" >> {params.stdout}
        fi
        """

# annotate snps with snpeff and save annotated_snps.vcf in snpcall_output directory
rule snpeff:
    input: clean = "snpcall_output/clean_snps.vcf"

    output:
        annotated = "snpcall_output/annotated_snps.vcf",
    
    params:
        stdout = "logs/stdout.txt",
        stderr = "logs/stderr.txt"

    shell:
        """
        java -Xmx3400m -jar {snpeff_jar} eff {genome_version} -dataDir /staging/leuven/stg_00079/teaching/snpeff_db {input.clean} > {output.annotated}

        header_length=$(grep -c '^#' {output.annotated})

        if [ $header_length -lt 35 ]; then
            echo "WARNING: annotated_snps.vcf looks too short" >> {params.stderr}
        else
            echo "annotated_snps.vcf created successfully" >> {params.stdout}
        fi
        """

# extract SNPs for APP, SOD1, and DYRK1A genes and save genes.vcf in main directory
rule extract_snps:
    input: 
        annotated = "snpcall_output/annotated_snps.vcf"

    output:
        extracted_snps = "genes.vcf"

    params:
        regions = "chr21:25880535-26171128,chr21:31659666-31668931,chr21:37365573-37526358",
        stdout = "logs/stdout.txt",
        stderr = "logs/stderr.txt"

    shell:
        """ 
        bgzip -c {input.annotated} > {input.annotated}.gz
        tabix -p vcf {input.annotated}.gz
        bcftools view -r {params.regions} {input.annotated}.gz > {output.extracted_snps}

        header_length=$(grep -c '^#' {output.extracted_snps})

        if [ $header_length -lt 37 ]; then
            echo "WARNING: genes.vcf looks too short" >> {params.stderr}
        else
            echo "genes.vcf created successfully" >> {params.stdout}
        fi
        """

rule barplot:
    input: "genes.vcf"
    output: 
        stats="barplot_stats/stats.txt",
        snp_counts="barplot_stats/snp_counts.txt",
        rep4 = report("barplot.png", category='SNPs', subcategory="SNPs per individual")
    shell:
        """
        bcftools stats -s- {input} > {output.stats}
        grep '^PSC' {output.stats} | awk '{{split($3, arr, "/"); split(arr[2], sample, "."); print sample[1], $5+$6}}' > {output.snp_counts}
        
        python -c "

df = pd.read_csv('barplot_stats/snp_counts.txt', sep=' ', header=None, names=['Sample_ID', 'SNP_Counts'])

plt.figure(figsize=(10, 6))
plt.bar(df['Sample_ID'], df['SNP_Counts'], color='skyblue')
plt.xlabel('Sample ID')
plt.ylabel('SNP Counts')
plt.title('SNP Counts per Individual')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

plt.savefig('barplot.png')
plt.close()
"
        """
