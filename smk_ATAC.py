
# Chang Liu

# Step 1 Fastqc

rule Fastqc:
    input:
    output:
    shell:
        fastqc -o -t   # o:output, t:thread, can give multiple files simutaneously and fastqc will test them separetely



rule trim_galore:

rule bowtie2:
    input:
    output:
    shell:
        bowtie2 --very-sensitive -X 2000 -p 8 -x /Sshare/home/liuchang/resource/bowtie2/mm10Index/mm10 \ # bowtie2 要在存放index的位置上多加一个物种名称（前缀！）
            -1 trim_galore/TsinghuaGZJ-CMH-GZJ-K2_FKDL210117656-1a-N705-N504_1_val_1.fq.gz \
            -2 trim_galore/TsinghuaGZJ-CMH-GZJ-K2_FKDL210117656-1a-N705-N504_2_val_2.fq.gz \
                -S atac_alignment_test.sam 

rule samtools:
    input:
    output:
    shell:
        samtools view -b -h -F 1028 -f 3 atac_alignment_test.sam -o atac_filter.bam  # -F 过滤read标准，-f 保留read标准，数字是组合，例如 3 = 1+2， 1028 = 4+1024


samtools flagstat atac_filter.bam > atac_filter.bam.stat \
    |samtools view -h -f 2 -q 30 atac_filter.bam.bam \
    |grep -v chrM \
    |samtools sort -O bam -@ 2 -o - > atac_filter.bam

samtools flagstat atac_filter.bam > atac_filter.bam.stat 
samtools view -h atac_filter.bam  \
    |grep -v chrM \
    |samtools sort -O bam  -@ 2 -o ./ > atac.last.bam \
    |samtools index atac.last.bam
samtools flagstat atac.last.bam > atac.last.stat