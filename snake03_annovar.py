import os

configfile: "config.yaml"
workdir: ".."

samples = config["SAMPLES"]
snp_vcf = config["REF"]["snp_vcf_dir"]
java_opt = "-Xmx10G"

rule all:
    input:
        expand("result/annovar/{sample}/table_annovar.hg38_multianno.txt", 
            sample=samples),
        expand("result/annovar/{sample}/geneanno.variant_function", 
            sample=samples),
        expand("result/annovar/{sample}/regionanno.hg38_cytoBand", 
            sample=samples),
        expand("result/annovar/{sample}/filter.hg38_gnomad211_exome_filtered", 
            sample=samples),

rule convert2annovar:
    input:
        vcf = "result/snp_indel/GermlineVcfFilter/{sample}/indel_snp.vcf.gz",
        anv = config["SOFT"]["annovar"] + "/convert2annovar.pl",
    output:
        vcf = "result/annovar/{sample}/convert2annovar.{sample}.avinput",
    params:
        pfx = "result/annovar/{sample}/convert2annovar",
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        {input.anv} --format vcf4 --filter pass {input.vcf} \
            --allsample --outfile {params.pfx}
        set +u; conda deactivate; set -u
        """

rule table_annovar:
    input:
        txt = rules.convert2annovar.output,
        db = config["SOFT"]["annovar"] + "/db/hg38",
        anv = config["SOFT"]["annovar"] + "/table_annovar.pl",
    output:
        txt = "result/annovar/{sample}/table_annovar.hg38_multianno.txt",
    params:
        pfx = "result/annovar/{sample}/table_annovar",
        ptc = "refGene,phastConsElements100way,cytoBand,gnomad211_exome,exac03,avsnp147,dbnsfp30a",
        opt = "g,r,r,f,f,f,f"
    threads: 5
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        {input.anv} {input.txt} {input.db} --buildver hg38 \
            --outfile {params.pfx} --remove --thread {threads} \
            --protocol {params.ptc} --operation {params.opt} \
            -nastring . -polish
        set +u; conda deactivate; set -u
        """

rule annotate_variation_geneanno:
    input:
        txt = rules.convert2annovar.output,
        db = rules.table_annovar.input.db,
        anv = config["SOFT"]["annovar"] + "/annotate_variation.pl",
    output:
        txt = "result/annovar/{sample}/geneanno.variant_function",
    params:
        pfx = "result/annovar/{sample}/geneanno",
    log: "result/annovar/{sample}/geneanno.log",
    threads: 28
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        {input.anv} {input.txt} {input.db} --geneanno --buildver hg38 \
            -dbtype refGene --thread {threads} --maxgenethread {threads} \
            --outfile {params.pfx}
        set +u; conda deactivate; set -u
        """

rule annotate_variation_regionanno:
    input:
        txt = rules.convert2annovar.output,
        db = rules.table_annovar.input.db,
        anv = config["SOFT"]["annovar"] + "/annotate_variation.pl",
    output:
        txt = "result/annovar/{sample}/regionanno.hg38_cytoBand",
    params:
        pfx = "result/annovar/{sample}/regionanno",
    log: "result/annovar/{sample}/regionanno.log",
    threads: 28
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        {input.anv} {input.txt} {input.db} --regionanno --buildver hg38 \
            -dbtype cytoBand --thread {threads} --maxgenethread {threads} \
            --outfile {params.pfx}
        set +u; conda deactivate; set -u
        """

rule annotate_variation_filter:
    input:
        txt = rules.convert2annovar.output,
        db = rules.table_annovar.input.db,
        anv = config["SOFT"]["annovar"] + "/annotate_variation.pl",
    output:
        txt = "result/annovar/{sample}/filter.hg38_gnomad211_exome_filtered",
    params:
        pfx = "result/annovar/{sample}/filter",
    log: "result/annovar/{sample}/filter.log",
    threads: 28
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        {input.anv} {input.txt} {input.db} --filter --buildver hg38 \
            -dbtype gnomad211_exome --thread {threads} --maxgenethread {threads} \
            --outfile {params.pfx}
        set +u; conda deactivate; set -u
        """