import os

configfile: "config.yaml"
workdir: ".."

samples = config["SAMPLES"]
snp_vcf = config["REF"]["snp_vcf_dir"]
java_opt = "-Xmx10G"

rule all:
    input:
        expand"result/GATK/{sample}/indel_snp.vcf.gz", sample=samples),

rule call_haplotype:
    input:
        bam = "result/mapping/{sample}/ApplyBqsr.bam",
        fa = config["REF"]["fa"],
        snp_db = snp_vcf + "/dbsnp_146.hg38.vcf.gz",
    output:
        tmp_dir = temp(directory("result/snp_indel/GermlineVcf/{sample}/_tmp")),
        vcf = "result/GATK/{sample}/HaplotypeCaller.vcf.gz",
    params:
        java_opt = java_opt,
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        if [[ -e {output.tmp_dir} ]]; then
            rm -rf {output.tmp_dir}
        fi
        mkdir -p {output.tmp_dir}
        gatk --java-options '{params.java_opt} -Djava.io.tmpdir={output.tmp_dir}' \
            HaplotypeCaller \
            -R {input.fa} -I {input.bam} -D {input.snp_db} -O {output.vcf} \
            -A QualByDepth \
            -A RMSMappingQuality \
            -A MappingQualityRankSumTest \
            -A ReadPosRankSumTest \
            -A FisherStrand \
            -A StrandOddsRatio \
            -A Coverage
        set +u; conda deactivate; set -u
        """

rule variant_recalibrator_snp:
    input:
        vcf = rules.call_haplotype.output.vcf,
        fa = config["REF"]["fa"],
        hapmap = snp_vcf + "/hapmap_3.3.hg38.vcf.gz",
        omni = snp_vcf + "/1000G_omni2.5.hg38.vcf.gz",
        phase1 = snp_vcf + "/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        snp_db = snp_vcf + "/dbsnp_146.hg38.vcf.gz",
    output:
        recal = "result/GATK/{sample}/snp.recal",
        tranches = "result/GATK/{sample}/snp.tranches",
        plotr = "result/GATK/{sample}/snp_plot.R",
    params:
        java_opt = java_opt,
        hapmap = "hapmap,known=false,training=true,truth=true,prior=15.0",
        omni = "omni,known=false,training=false,truth=false,prior=12.0",
        G1000 = "1000G,known=false,training=true,truth=false,prior=10.0",
        snp_db = "dbsnp,known=true,training=false,truth=false,prior=2.0",
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        gatk --java-options '{params.java_opt}' VariantRecalibrator \
            -R {input.fa} -V {input.vcf} -mode SNP  -O {output.recal} \
            --resource:{params.hapmap} {input.hapmap} \
            --resource:{params.omni} {input.omni} \
            --resource:{params.G1000} {input.phase1} \
            --resource:{params.snp_db} {input.snp_db} \
            -an QD -an MQ -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
            --tranches-file {output.tranches} \
            --rscript-file {output.plotr}
        set +u; conda deactivate; set -u
        """

rule apply_vqsr_snp:
    input:
        vcf = rules.call_haplotype.output.vcf,
        fa = config["REF"]["fa"],
        recal = rules.variant_recalibrator_snp.output.recal,
        tranches = rules.variant_recalibrator_snp.output.tranches,
    output:
        vcf = "result/GATK/{sample}/snp.vcf.gz",
    params:
        java_opt = java_opt,
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        gatk --java-options '{params.java_opt}' ApplyVQSR \
            -R {input.fa} -V {input.vcf} -O {output.vcf} \
            --truth-sensitivity-filter-level 99.0 \
            --tranches-file {input.tranches} \
            --recal-file {input.recal} \
            -mode SNP
        set +u; conda deactivate; set -u
        """

rule variant_recalibrator_indel:
    input:
        vcf = rules.call_haplotype.output.vcf,
        fa = config["REF"]["fa"],
        mills = snp_vcf + "/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
    output:
        recal = "result/GATK/{sample}/indel.recal",
        tranches = "result/GATK/{sample}/indel.tranches",
        plotr = "result/GATK/{sample}/indel_plot.R",
    params:
        java_opt = java_opt,
        mills = "mills,known=true,training=true,truth=true,prior=12.0",
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        gatk --java-options '{params.java_opt}' VariantRecalibrator \
            -R {input.fa} -V {input.vcf} -O {output.recal} \
            -resource:{params.mills} {input.mills} \
            -an QD -an MQ -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
            -mode INDEL --max-gaussians 6 \
            --tranches-file {output.tranches} \
            --rscript-file {output.plotr}
        set +u; conda deactivate; set -u
        """
    
rule apply_vqsr_indel:
    input:
        vcf = rules.call_haplotype.output.vcf,
        fa = config["REF"]["fa"],
        recal = rules.variant_recalibrator_indel.output.recal,
        tranches = rules.variant_recalibrator_indel.output.tranches,
    output:
        vcf = "result/GATK/{sample}/indel.vcf.gz",
    params:
        java_opt = java_opt,
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        gatk --java-options '{params.java_opt}' ApplyVQSR \
            -R {input.fa} -V {input.vcf} -O {output.vcf} \
            --truth-sensitivity-filter-level 99.0 \
            --tranches-file {input.tranches} \
            --recal-file {input.recal} \
            -mode INDEL
        set +u; conda deactivate; set -u
        """

rule merge_snp_indel:
    input:
        snp = rules.apply_vqsr_snp.output,
        indel = rules.apply_vqsr_indel.output,
    output:
        merge = "result/GATK/{sample}/indel_snp.vcf.gz",
    shell: """
        set +u; source ~/miniconda3/bin/activate seq; set -u
        gatk MergeVcfs -I {input.snp} -I {input.indel} -O {output.merge}
        set +u; conda deactivate; set -u
        """