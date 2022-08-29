#!/usr/bin/bash
#MISET-seq
#Fuqing Jiang

Main=/MISET-seq
sampleId="E15_5-S1"
Data=$Main/Data
Result=$Main/Result
Scripts=$Main/Scripts

ref='/MISET-seq/reference/refdata-cellranger-arc-mm10-2020-A-2.0.0'

if [ ! -d $Data/$sampleId ];then
  mkdir $Data/$sampleId
fi

#Process RNA Barcode
python $Scripts/MISET_Split_BC_RNA.py -i $Data/$sampleId"-RNA.R2.fastq.gz" -o $Data/$sampleId/$sampleId"_gex_S1_L001_R1_001.fastq"

pigz -p 100 -f $Data/$sampleId/$sampleId"_gex_S1_L001_R1_001.fastq" 

cp $Data/$sampleId"-RNA.R1.fastq.gz" $Data/$sampleId/$sampleId"_gex_S1_L001_R2_001.fastq.gz"

#Process ATAC Barcode
python $Scripts/MISET_Split_BC_ATAC.py -i $Data/$sampleId"-ATAC.R2.fastq.gz" -o1 $Data/$sampleId/$sampleId"_atac_S1_L001_R3_001.fastq" -o2 $Data/$sampleId/$sampleId"_atac_S1_L001_R2_001.fastq"

pigz -p 100 -f $Data/$sampleId/$sampleId"_atac_S1_L001_R3_001.fastq"

pigz -p 100 -f $Data/$sampleId/$sampleId"_atac_S1_L001_R2_001.fastq"

cp $Data/$sampleId"-ATAC.R1.fastq.gz" $Data/$sampleId/$sampleId"_atac_S1_L001_R1_001.fastq.gz"

#Create ${sampleId}_libraries.csv file for cellranger-arc
echo -e "fastqs,sample,library_type \n$Data/$sampleId,${sampleId}_gex,Gene Expression \n$Data/$sampleId,${sampleId}_atac,Chromatin Accessibility" > $Data/${sampleId}_libraries.csv

cd $Result

cellranger-arc count --id=$sampleId \
                   --reference=$ref \
                   --libraries=$Data/${sampleId}_libraries.csv \
                   --localcores=50 \
                   --localmem=64
                   #Default is counting of intronic reads, use --gex-exclude-introns to disable.

#Transcform the bam to bw file
cd $Result/$sampleId/outs/

bamCoverage --bam atac_possorted_bam.bam -o $sampleId"_atac.bw" --normalizeUsing CPM --extendReads 170 --binSize 50

bamCoverage --bam gex_possorted_bam.bam -o $sampleId"_gex.bw" --normalizeUsing CPM --extendReads 170 --binSize 50

#draw plot of TSS
gtf=$ref/genes/genes.gtf
computeMatrix reference-point -S $sampleId"_atac.bw" -R $gtf --skipZeros -o  ${sampleId}_atac_tss.mat.gz -p 30 -a 3000 -b 3000 --missingDataAsZero --referencePoint TSS

plotHeatmap --dpi 300 -m ${sampleId}_atac_tss.mat.gz -out ${sampleId}_atac_tss_heatmap.png --sortUsing sum --startLabel "Peak Start" --endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --missingDataColor "#cc2200" --samplesLabel "${sampleId}"

echo [done] `date +"%Y-%m-%d %H:%M:%S"` 