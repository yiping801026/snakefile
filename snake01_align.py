import os

configfile: "config.yaml"
workdir: ".."

samples = config["SAMPLES"]
reads = ["R1", "R2"]
genomes = ["rtRNA", "hg38"]

rule all:
    input:
        "result/reads_QC/fastqc",
        expand("result/mapping/{genome}/{sample}.unmaped.fq.gz", genome=genomes, sample=samples),
        expand("result/reads_QC/check_ss/{sample}_strand_specificity.log", sample=samples),

rule srr2fq:
    input:
        soft = "/beegfs/zhoulab/yangjiayi/software/sratoolkit/bin/fasterq-dump",
        srr = "data/{sample}",
    output: 
        fq = temp("result/reads_QC/fasterq_dump/{sample}.fastq"),
    threads: 14
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        {input.soft} {input} -e {threads} --split-spot --stdout > {output.fq}
        set +u; conda deactivate; set -u
        """

rule fastqc:
    input:
        expand(rules.srr2fq.output, sample=samples),
    output:
        directory("result/reads_QC/fastqc"),
    threads: 8
    shell:"""
        set +u; source ~/miniconda3/bin/activate seq; set -u
        mkdir -p {output}
        fastqc -o {output} --threads {threads} {input}
        multiqc -o {output} {output}
        set +u; conda deactivate; set -u
        """

rule rtRNA_mapping:
    input:
        fq = rules.srr2fq.output.fq,
    output:
        un = temp("result/mapping/rtRNA/{sample}.unmaped.fq"),
    params:
        bt_index = "/beegfs/zhoulab/yangjiayi/reference/hg38/bt_index/rtRNA/rtRNA",
    log: "result/mapping/rtRNA/{sample}.bowtie.log",
    threads: 14
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        bowtie {params.bt_index} {input.fq} \
            -p {threads} -v 2 -k 1 --best \
            --chunkmbs 250 --un {output.un} \
            > /dev/null 2> {log}
        set +u; conda deactivate; set -u
        """

rule hg38_mapping:
    input:
        fq = rules.rtRNA_mapping.output.un,
    output:
        bam = "result/mapping/hg38/{sample}.sort.bam",
        un = temp("result/mapping/hg38/{sample}.unmaped.fq"),
        stat = "result/mapping/hg38/{sample}.flagstat.log",
    params:
        bt_index = "/beegfs/zhoulab/yangjiayi/reference/hg38/bt_index/hg38/hg38",
    log: "result/mapping/hg38/{sample}.bowtie.log",
    threads: 14
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        bowtie {params.bt_index} {input.fq} \
            -p {threads} -v 2 -k 20 --best --strata \
            -S --chunkmbs 250 --un {output.un} \
            2> {log} |
        samtools view -hbS -F 256 -F 4 |
        samtools sort -@ {threads} -o {output.bam}
        samtools index -@ {threads} {output.bam}
        samtools flagstat -@ {threads} {output.bam} > {output.stat}
        set +u; conda deactivate; set -u
        """

rule zip_unmap:
    input: "result/mapping/{genome}/{sample}.unmaped.fq",
    output: "result/mapping/{genome}/{sample}.unmaped.fq.gz",
    threads: 14
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        pigz -p {threads} -c {input} > {output}
        set +u; conda deactivate; set -u
        """

#ss:"1+-,1-+,2++,2--"
rule check_ss:
    input:
        bed = config["REF"]["bed"],
        bam = rules.hg38_mapping.output.bam,
    output:
        res = "result/reads_QC/check_ss/{sample}_strand_specificity.log",
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        echo "sample name:{wildcards.sample}" > {output.res}
        infer_experiment.py -r {input.bed} -i {input.bam} \
            -s 1000000 >> {output.res} 2>&1
        set +u; conda deactivate; set -u
        """
