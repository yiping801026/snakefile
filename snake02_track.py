import os

configfile: "config.yaml"
workdir: ".."

samples = config["SAMPLES"]
reads_nums = config["READS"]
colors = config["COLORS"]

rule all:
    input:
        expand("result/track/bw/{sample}.bw", sample=samples),
        expand("result/track/norm2mean/{sample}_norm.bedGraph.gz", sample=samples),

rule bam2bw:
    input:
        bam = "result/mapping/hg38/{sample}.sort.bam",
        chr_info = config["REF"]["chrsize"]["hg38"],
    output:
        wig = temp("result/track/bw/{sample}.wig"),
        bw = "result/track/bw/{sample}.bw",
    log: "result/track/bw/{sample}_bam2wig.log",
    params:
        wig_pfx = "result/track/bw/{sample}",
        strand = "",
    shell:"""
        set +u; source ~/miniconda3/bin/activate seq; set -u
        echo "input bam file: {input.bam}" > {log}
        bam2wig.py -i {input.bam} -s {input.chr_info} \
            -o {params.wig_pfx} -u >> {log} 2>&1
        set +u; conda deactivate; set -u
        """

rule bw2bg:
    input:
        bw = rules.bam2bw.output.bw,
    output:
        bg = temp("result/track/norm2mean/{sample}_nh.bedGraph"),
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        bigWigToBedGraph {input.bw} {output.bg}
        set +u; conda deactivate; set -u
        """

rule normalization:
    input:
        bg = rules.bw2bg.output.bg,
    output:
        norm = temp("result/track/norm2mean/{sample}_nh_norm.bedGraph"),
    params:
        scale = lambda wildcards: reads_nums["normto"] / reads_nums[wildcards.sample],
    run:
        with open(input.bg, 'r') as fi, open(output.norm, 'w') as fo:
            for line in fi:
                record = line.strip().split('\t')
                normed = float(record[-1]) * params.scale
                record[-1] = format(abs(normed), ".2f")
                print(*record, sep='\t', file=fo)

rule set_color:
    input:
        raw = rules.normalization.output.norm,
    output:
        header = temp("result/track/norm2mean/{sample}_norm.header.txt"),
        gz = "result/track/norm2mean/{sample}_norm.bedGraph.gz",
    params:
        color = lambda wildcards: colors[wildcards.sample],
    run:
        with open(output.header, 'w') as f:
            f.write("track type=bedGraph ")
            f.write('name="%s" '%(wildcards.sample))
            f.write("visibility=full ")
            f.write("color=%s "%(params.color))
            f.write("AutoScale=on ")
            f.write("alwaysZero=on ")
            f.write("maxHeightPixels=50:50:50\n")
        shell("cat {output.header} {input.raw} | gzip -c - > {output.gz}")
