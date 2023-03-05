import os

configfile: "config.yaml"
workdir: ".."

samples = config["SAMPLES"]
snp_vcf = config["REF"]["snp_vcf_dir"]
java_opt = "-Xmx10G"

rule all:
    input:
        expand("result/fastqc/{sample}", sample=samples),
        expand("result/mapping/{sample}/ApplyBqsr.bam", sample=samples),

rule cutadapt:
    input:
        fq1 = "data/{sample}_R1.fq.gz",
        fq2 = "data/{sample}_R2.fq.gz",
    output:
        fq1 = "result/cutadapt/{sample}_cutadapt_R1.fq",
        fq2 = "result/cutadapt/{sample}_cutadapt_R2.fq",
        st1 = "result/cutadapt/{sample}_tooshort_R1.fq.gz",
        st2 = "result/cutadapt/{sample}_tooshort_R2.fq.gz",
    params:
        adp_seq1 = "AGATCGGAAGAGCACACGTCTGAAC",
        adp_seq2 = "AGATCGGAAGAGCGTCGTGTAGGGA",
    log: "result/cutadapt/{sample}_cutadapt.log",
    threads: 5
    shell:"""
        set +u; source ~/miniconda3/bin/activate seq; set -u
        cutadapt -q 20 --max-n 0.1 -m 75 -j {threads} \
            -a {params.adp_seq1} -A {params.adp_seq2} \
            -o {output.fq1} -p {output.fq2} \
            {input.fq1} {input.fq2} \
            --too-short-output {output.st1} \
            --too-short-paired-output {output.st2} \
            > {log}
        set +u; conda deactivate; set -u
        """

rule fastqc:
    input:
        rules.cutadapt.input,
        rules.cutadapt.output.fq1,
        rules.cutadapt.output.fq2,
    output:
        directory("result/fastqc/{sample}"),
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        mkdir -p {output}
        fastqc -o {output} {input}
        set +u; conda deactivate; set -u
        """

rule bwa_mapping:
    input:
        fq1 = rules.cutadapt.output.fq1,
        fq2 = rules.cutadapt.output.fq2,
    output:
        bam = "result/mapping/{sample}/Aligned.sort.BWA.bam",
        stat = "result/mapping/{sample}/Aligned.sort.flagstat.log",
    params:
        bwt_pfx = config["REF"]["bwa_pfx"],
        header = "@RG\\tID:{sample}\\tLB:{sample}_WES_L4\\tPL:ILLUMINA\\tSM:{sample}",
    threads: 5
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        bwa mem -t {threads} -M -R '{params.header}' \
            {params.bwt_pfx} {input.fq1} {input.fq2} |
        samtools sort -@ {threads} -o {output.bam}
        samtools index -@ {threads} {output.bam}
        samtools flagstat -@ {threads} {output.bam} > {output.stat}
        set +u; conda deactivate; set -u
        """

rule build_dict:
    input: config["REF"]["fa"],
    output: config["REF"]["gatk_dict"],
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        samtools faidx {input}
        gatk CreateSequenceDictionary -R {input} -O {output}
        set +u; conda deactivate; set -u
        """

rule mark_dup:
    input: rules.bwa_mapping.output.bam,
    output:
        bam = "result/mapping/{sample}/Markdup.bam",
        mtx = "result/mapping/{sample}/Markdup.metrics.txt",
    params:
        java_opt = java_opt,
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        gatk --java-options '{params.java_opt}' MarkDuplicates \
            -I {input} -O {output.bam} -M {output.mtx}
        set +u; conda deactivate; set -u
        """

rule fix_mate:
    input: rules.mark_dup.output.bam,
    output: "result/mapping/{sample}/FixMate.bam",
    params:
        java_opt = java_opt,
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        gatk --java-options '{params.java_opt}' FixMateInformation \
            -I {input} -O {output} -SO coordinate
        samtools index {output}
        set +u; conda deactivate; set -u
        """

rule base_recall:
    input:
        bam = rules.fix_mate.output,
        fa = rules.build_dict.input,
        gatk_dict = rules.build_dict.output,
        snp_db = snp_vcf + "/dbsnp_146.hg38.vcf.gz",
        indel_db = snp_vcf + "/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
    output: 
        bqsr = "result/mapping/{sample}/BaseRecalibrator.txt",
    params:
        java_opt = java_opt,
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        gatk --java-options '{params.java_opt}' BaseRecalibrator \
            -R {input.fa} -I {input.bam} -O {output} \
            --known-sites {input.snp_db} \
            --known-sites {input.indel_db} 
        set +u; conda deactivate; set -u
        """

rule apply_bqsr:
    input:
        bam = rules.fix_mate.output,
        bqsr = rules.base_recall.output,
        fa = rules.build_dict.input,
    output:
        bam = "result/mapping/{sample}/ApplyBqsr.bam",
    params:
        java_opt = java_opt,
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        gatk --java-options '{params.java_opt}' ApplyBQSR \
            -R {input.fa} -I {input.bam} -bqsr {input.bqsr} -O {output.bam}
        set +u; conda deactivate; set -u
        """
    
