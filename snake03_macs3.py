import os

configfile: "config.yaml"
workdir: ".."

samples = config["SAMPLES"]
reads = ["R1", "R2"]
genomes = ["rtRNA", "hg38"]

rule all:
    input:
        expand("result/macs3/predictd/{sample}", sample=samples)

# rule predictd:
#     input:
#         bam = "result/mapping/hg38/{sample}.sort.bam",
#     output:
#         directory("result/macs3/predictd/{sample}"),
#     log: "result/macs3/predictd/{sample}.log",
#     threads: 1
#     shell: """
#         set +u; source ~/miniconda3/bin/activate macs3; set -u
#         macs3 predictd -i {input.bam} -g hs --outdir {output} \
#             > {log} 2>&1
#         """

rule callpeak:
    input:
        t1 = "result/mapping/hg38/{sample}1_input.sort.bam",
        t2 = "result/mapping/hg38/{sample}2_input.sort.bam",
        c1 = "result/mapping/hg38/{sample}1_IP.sort.bam",
        c2 = "result/mapping/hg38/{sample}2_IP.sort.bam",
    output:
        directory("result/callpeak/macs3/{sample}"),
    log: "result/callpeak/macs3/{sample}.log",
    threads: 5
    shell: """
        set +u; source ~/miniconda3/bin/activate macs3; set -u
        macs3 callpeak --format BAM \
            -t {input.t1} {input.t2} \
            -t {input.c1} {input.c2} \
            --outdir {output} \
            -n {wildcards.sample} \
            -B --SPMR \
            --pvalue 0.01 --gsize hs \
            --seed 89 \
            > {log} 2>&1
        set +u; conda deactivate; set -u
        """

rule plotProfile_ARS:
    input:
        bw_t = "result/bigwig/{treat}.bw",
        bw_c = "result/bigwig/Input.bw",
        bed = rules.get_bed_all.output.ars_bed,
    output:
        mat = "result/profile_plot/{treat}/{treat}_allgene.mat.gz",
        pdf = "result/profile_plot/{treat}/{treat}_allgene.pdf",
        png = "result/profile_plot/{treat}/{treat}_allgene.png",
    threads: 1
    shell:"""
        set +u; source ~/miniconda3/bin/activate bowtie2; set -u
        computeMatrix reference-point -p {threads} \
            --referencePoint TSS -b 2000 -a 2000 \
            -S {input.bw_t} {input.bw_c} -R {input.bed} \
            --skipZeros -out {output.mat} &&
        plotProfile --dpi 300 -m {output.mat} -out {output.pdf} \
            --plotFileFormat pdf --perGroup \
            --plotTitle "{wildcards.treat} ChIP profile plot" \
            --yAxisLabel "PolII singnal density" --legendLocation best &&
        plotHeatmap -m {output.mat} -out {output.png} \
            --colorMap Blues --whatToShow 'heatmap and colorbar' \
            --zMin -4 -4 0 --zMax 0 0 5 --dpi 300 --boxAroundHeatmaps no
        set +u; conda deactivate; set -u
        """