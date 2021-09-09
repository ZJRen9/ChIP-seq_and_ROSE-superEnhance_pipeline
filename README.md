ChIP-seq_and_ROSE-superEnhance_pipeline
===

step1:mapping
---
    bowtie2-align-s -p 6 --no-mixed --no-discordant --no-unal -x mm10.genome -1 file_chip_R1.fastq -2 file_chip_R2.fastq -S file_chip.sam
    bowtie2-align-s -p 6 --no-mixed --no-discordant --no-unal -x mm10.genome -1 file_input_R1.fastq -2 file_input_R2.fastq -S file_input.sam

step2:mapping quality cut off and duplicate remove
---
    samtools view -h -b -q 30 file_chip.sam >> file_chip.bam
    samtools view -h -b -q 30 file_input.sam >> file_input.bam

    samtools sort -o file_chip.bam.sorted file_chip.bam
    samtools sort -o file_input.bam.sorted file_input.bam

    samtools index file_chip.bam.sorted
    samtools index file_input.bam.sorted

    samtools rmdup -S file_chip.bam.sorted file_chip.bam.sorted.rmdup
    samtools rmdup -S file_input.bam.sorted file_input.bam.sorted.rmdup

    samtools index file_chip.bam.sorted.rmdup
    samtools index file_input.bam.sorted.rmdup

step3:peak calling
---
    macs2 callpeak -t file_chip.chip.bam.sorted.rmdup -c file_input.bam.sorted.rmdup -B -f BAMPE -g mm -q 0.01 -n file

step4:peak select for single replicate
---
    awk '$7 >9 && $8 >3 {print $1"\t"$2"\t"$3"\t"$10}' file_peaks.xls >> file_peaks.select.bed
    bedtools intersect -v -a file_peakks.select.bed -b mm10_black.bed >> file_peaks.select.clean.bed

step5:RPKM count in bed region and HQ peak select
---
    samtools sort -n -o file_chip.bam.sorted.rmdup.namesort file_chip.bam.sorted.rmdup
    
    bedtools bamtobed -bedpe -i file_chip.bam.sorted.rmdup.namesort >> file_chip.bam.sorted.rmdup.PE.bed

    awk '{print $1"\t"$2"\t"$6"\t"$7}' file_chip.bam.sorted.rmdup.PE.bed >> file_chip.bam.sorted.rmdup.bed

    bedSort file_chip.bam.sorted.rmdup.bed file_chip.bam.sorted.rmdup.bed.sorted

    bedtools coverage -a file_peaks.select.clean.bed -b file_chip.bam.sorted.rmdup.bed.sorted >> file_peaks.select.clean.bed.chip

    python RPKM_count.py file_peaks.select.clean.bed.chip file_chip.tatal_read_counts file_peaks.select.clean.bed >> file_peaks.select.clean.bed.chip.RPKM

    awk '$5 >4 {print $1"\t"$2"\t"$3"\t"$4}' file_peaks.select.clean.bed.chip.RPKM >> file_peaks.select.clean.HQ.bed

step6:Super-Enhance calling
---
    python ROSE_main.py -g MM10 -i file_peak.select.clean.HQ.bed -r file_chip.bam.sorted.rmdup -c file_input.bam.sorted.rmdup -o file_Super -t 2500

step7:Super-Enhance annotation
---
    /usr/bin/python ROSE_geneMapper.py -g MM10 -i file_Super/file_peak.select.clean.HQ_SuperEnhancers.table.txt -o file_peak.select.clean.HQ_SuperEnhancers.annotation
step8:Most suitable genen select
---
    python suitable_annotation_gene_select.py file_peak.select.clean.HQ_SuperEnhancers.annotate gene_exp.txt >> file_peak.select.clean.HQ_SuperEnhancers.annotate.suitable
