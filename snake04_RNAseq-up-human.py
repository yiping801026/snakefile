import os

configfile: "config.yaml"

samples = config["SAMPLES"]

workdir: "."

RNAseq = "/Sshare/home/liuchang/miniconda3/envs/RNAseq/"

rule all:
    input:
        expand("../results/mapping/{sample}/Aligned.out.bam", sample=samples),
        expand("../results/mapping/{sample}/Unmapped_R{read}.fq.gz", sample=samples, read=[1,2]),
        expand("../results/mapping/{sample}/samtools_flagstat.log", sample=samples),
        expand("../results/mapping/{sample}/primary_sort_rmdup.bam", sample=samples),
        expand("../results/featurecounts/{sample}_featureCounts.txt",sample=samples),
        "../results/featurecounts/featureCounts_all.txt",


rule star_index: #创建参考基因组索引
    input:
        gtf = config["REF"]["gtf"],
        fa = config["REF"]["fa"],
    output:
        index = directory("../results/index"),
        log = "../results/index/log.txt"
    threads: 10
    shell:"""
        set +u; source ~/miniconda3/bin/activate {RNAseq}; set -u
        echo "start build index `date "+%Y-%m-%d %H:%M:%S"`" > {output.log}
        STAR --runMode genomeGenerate \
            --runThreadN  {threads}\
            --genomeDir {output.index} \
            --genomeFastaFiles {input.fa} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang  149
        echo "end build index `date "+%Y-%m-%d %H:%M:%S"`" >> {output.log}
        set +u; conda deactivate; set -u
    """

rule reads_unzip: 
     input: "/fshare2/linlab/SunXin/N2011852_WX_80-576948958_eukRNASEQ/result/fastqc/data/00data/{sample}/{sample}_combined_{read}.fastq.gz",
     output: temp("../data/raw-data/{sample}_{read}.fq"),
     threads: 4
     shell: """
         set +u; source ~/miniconda3/bin/activate {RNAseq}; set -u
         pigz -p {threads} -dc {input} > {output}
         set +u; conda deactivate; set -u
         """

rule mapping: #序列比对
    input:
        index = config["REF"]["STARindex"]["mm"],
        gtf = config["REF"]["gtf"],
        fq1 = "../data/raw-data/{sample}_R1.fq",
        fq2 = "../data/raw-data/{sample}_R2.fq",
        log = rules.star_index.output.log
    output:
        bam = "../results/mapping/{sample}/Aligned.out.bam",
        un1 = temp("../results/mapping/{sample}/Unmapped.out.mate1"), #mate是指read1对应的read2
        un2 = temp("../results/mapping/{sample}/Unmapped.out.mate2"),
    params:
        outdir = "../results/mapping/{sample}",
    threads: 8
    shell: """
        set +u; source ~/miniconda3/bin/activate {RNAseq}; set -u
        if [[ -e {params.outdir} ]]; then
            rm -r {params.outdir}
        fi
        mkdir -p {params.outdir}
        cp {input.log} {params.outdir}
        STAR --runThreadN {threads} \
            --genomeDir {input.index} \
            --sjdbGTFfile {input.gtf} \
            --sjdbScore 1 \
            --readFilesIn {input.fq1} {input.fq2} \
            --outFileNamePrefix {params.outdir}/ \
            --outSAMattributes NH HI AS NM MD \
            --outSAMunmapped None \
            --outFilterType BySJout \
            --outFilterMultimapNmax 20 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outReadsUnmapped Fastx \
            --outSAMtype BAM Unsorted 
        set +u; conda deactivate; set -u
        """

rule zip_unmap: #对未匹配的reads进行打包压缩
    input: "../results/mapping/{sample}/Unmapped.out.mate{read}",
    output: "../results/mapping/{sample}/Unmapped_R{read}.fq.gz",
    threads: 4
    shell: """
        set +u; source ~/miniconda3/bin/activate {RNAseq}; set -u
        pigz -p {threads} -c {input} > {output}
        set +u; conda deactivate; set -u
        """

  
rule sort_bam: #对bam文件进行排序
    input:
        bam = rules.mapping.output.bam,
    output:
        bam = temp("../results/mapping/{sample}/uniq_sort.bam"),
    threads: 4
    shell:"""
        set +u; source ~/miniconda3/bin/activate {RNAseq}; set -u
        samtools view -@ {threads} -hb -q 255 {input.bam} |
        samtools sort -@ {threads} -o {output.bam}
        set +u; conda deactivate; set -u
        """

# Remove duplication reads
rule picard_rmdup: #移除重复的reads
    input:
        bam = rules.sort_bam.output.bam,
    output:
        bam = temp("../results/mapping/{sample}/rmdup.bam"),
        mark = "../results/mapping/{sample}/rmdup_markdup.txt",
    log: "../results/mapping/{sample}/rmdup.log",
    shell:""" 
        set +u; source ~/miniconda3/bin/activate {RNAseq}; set -u
        picard MarkDuplicates \
            INPUT={input.bam} \
            OUTPUT={output.bam} \
            METRICS_FILE={output.mark} \
            REMOVE_DUPLICATES=true \
            ASSUME_SORTED=true 2> {log}
        set +u; conda deactivate; set -u
        """

rule index_bam:
    input:
        bam = rules.picard_rmdup.output.bam,
    output:
        bam = "../results/mapping/{sample}/primary_sort_rmdup.bam",
        bai = "../results/mapping/{sample}/primary_sort_rmdup.bam.bai",
    threads: 4
    shell:"""
        set +u; source ~/miniconda3/bin/activate {RNAseq}; set -u
        samtools sort -@ {threads} -o {output.bam} {input.bam}
        samtools index -@ {threads} {output.bam} #这里为什么不用输出output.bai？
        set +u; conda deactivate; set -u
        """

rule flagstat: #结果总结
    input:
        bam = rules.index_bam.output.bam,
    output:
        log = "../results/mapping/{sample}/samtools_flagstat.log",
    shell:"""
        set +u; source ~/miniconda3/bin/activate {RNAseq}; set -u
        samtools flagstat {input.bam} > {output.log}
        set +u; conda deactivate; set -u
        """


rule featurecounts: 
    input:
        bam = "../results/mapping/{sample}/primary_sort_rmdup.bam",
        gtf = config["REF"]["gtf"],
    output:
        txt = "../results/featurecounts/{sample}_featureCounts.txt",
    log: "../results/featurecounts/{sample}_featureCounts.log",
    threads: 4
    shell:"""
        set +u; source ~/miniconda3/bin/activate {RNAseq}; set -u
        featureCounts -a {input.gtf} -o {output.txt} {input.bam} \
            -g gene_id -s 0 -T {threads} -p > {log} 2>&1
        set +u; conda deactivate; set -u
        """

rule featurecountsAll: #对gene-level上计数得到的文件进行整合
    input: 
        bam = expand("../results/mapping/{sample}/primary_sort_rmdup.bam",sample = samples),
        bai = expand("../results/mapping/{sample}/primary_sort_rmdup.bam.bai",sample = samples),
        gtf = config["REF"]["gtf"],
    output:
        txt = "../results/featurecounts/featureCounts_all.txt",
    log: "../results/featurecounts/featureCounts_all.log",
    threads: 12
    shell: """
        set +u; source ~/miniconda3/bin/activate {RNAseq}; set -u
        featureCounts -a {input.gtf} -o {output.txt} {input.bam} \
            -t gene -g gene_id -s 0 -T {threads} -p > {log} 2>&1
        set +u; conda deactivate; set -u
        """