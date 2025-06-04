import glob
import re

files = glob.glob("*_R1.fastq.gz")
SAMPLES = [re.match(r"(.*)_R1.fastq.gz", f).group(1) for f in files]
SAMPLES = sorted(set(SAMPLES))

rule all:
    input: expand("Paired/{sample}_paired_R1.fastq.gz", sample=SAMPLES),
           expand("Unpaired/{sample}_unpaired_R1.fastq.gz", sample=SAMPLES),
           expand("Paired/{sample}_paired_R2.fastq.gz", sample=SAMPLES),
           expand("Unpaired/{sample}_unpaired_R2.fastq.gz", sample=SAMPLES),
           expand("TrimQC/{sample}_paired_R1_fastqc.html", sample=SAMPLES),
           expand("TrimQC/{sample}_paired_R1_fastqc.zip", sample=SAMPLES),
           expand("TrimQC/{sample}_paired_R2_fastqc.html", sample=SAMPLES),
           expand("TrimQC/{sample}_paired_R2_fastqc.zip", sample=SAMPLES),
           "hiv1cref.fasta.amb",
            "hiv1cref.fasta.ann",
            "hiv1cref.fasta.bwt",
            "hiv1cref.fasta.pac",
            "hiv1cref.fasta.sa",
            expand("Alignment/{sample}_aligned.bam", sample=SAMPLES),
            expand("Alignment/{sample}_aligned_sorted.bam", sample=SAMPLES),
            expand("Alignment/{sample}_alignment_stats.txt", sample=SAMPLES),
            expand("Variants/{sample}_variants.tsv", sample=SAMPLES)
          
rule trimming:
    input: R1="{sample}_R1.fastq.gz",
           R2="{sample}_R2.fastq.gz"
    output: R1p="Paired/{sample}_paired_R1.fastq.gz",
            R1u="Unpaired/{sample}_unpaired_R1.fastq.gz",
            R2p="Paired/{sample}_paired_R2.fastq.gz",
            R2u="Unpaired/{sample}_unpaired_R2.fastq.gz"
    shell:
        "java -jar /home/sj/Installations/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 {input.R1} {input.R2} {output.R1p} {output.R1u} {output.R2p} {output.R2u} ILLUMINACLIP:/home/sj/Installations/Trimmomatic-0.39/adapters/NEBnext-PE.fa:2:30:10 SLIDINGWINDOW:4:30 MINLEN:100"



rule Trim_qc_R1:
    input: "Paired/{sample}_paired_R1.fastq.gz",
    output: "TrimQC/{sample}_paired_R1_fastqc.html",
            "TrimQC/{sample}_paired_R1_fastqc.zip"
    threads: 2
    shell:
        "fastqc -t 2 --outdir=TrimQC {input}"

rule trim_qc_R2:
    input:  "Paired/{sample}_paired_R2.fastq.gz"
    output: "TrimQC/{sample}_paired_R2_fastqc.html",
            "TrimQC/{sample}_paired_R2_fastqc.zip"
    threads: 2
    shell:
        "fastqc -t 2 --outdir=TrimQC {input}"

rule index:
    input: "hiv1cref.fasta"
    output: "hiv1cref.fasta.amb",
            "hiv1cref.fasta.ann",
            "hiv1cref.fasta.bwt",
            "hiv1cref.fasta.pac",
            "hiv1cref.fasta.sa"
    shell:
        "bwa index {input}"

rule bwa_align:
    input: ref="hiv1cref.fasta",
           s1="Paired/{sample}_paired_R1.fastq.gz",
           s2="Paired/{sample}_paired_R2.fastq.gz"
    output:
        "Alignment/{sample}_aligned.bam"
    threads: 4
    shell:
        "bwa mem -t 4 {input.ref} {input.s1} {input.s2} | samtools view -F 4 -o {output}"

rule sorting:
    input: "Alignment/{sample}_aligned.bam"
    output: "Alignment/{sample}_aligned_sorted.bam"
    shell:
        "samtools sort {input} -o {output}"
            
rule alignment_statistics:
    input: "Alignment/{sample}_aligned_sorted.bam"
    output: "Alignment/{sample}_alignment_stats.txt"
    shell:
        "samtools flagstat {input} >> {output}"

rule variant_calling:
    input: s="Alignment/{sample}_aligned_sorted.bam",
            ref="hiv1cref.fasta",
            gf="hiv1cref.gff3"
    output: "Variants/{sample}_variants.tsv"
    shell:
        "samtools mpileup -aa -A -d 0 -B -Q 0 {input.s} | ivar variants -p {output} -q 30 -t 0.001 -r {input.ref} -g {input.gf}"
        
