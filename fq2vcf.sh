#!/bin/bash

path_1=$1

red="\033[31;40m"
green="\033[32;40m"
none="\033[0m"

if [ "$1" == "-h" ]; then

  echo
  echo -e $red"Usage:" $green"fq2vcf.sh" $yellow"path"
  echo
  echo -e $yellow"fq2vcf is a shell script to automatically analyze exome data from trimmed fastq file to vcf file \n followed by annotation by ANNOVAR. Fastq files with fq, fastq, fq.gz, fastq.gz extensions are accepted." 
  exit 0
fi

rm -rf fq2vcf_output > /dev/null
mkdir -p fq2vcf_output
mkdir -p fq2vcf_output/sam
mkdir -p fq2vcf_output/bam
mkdir -p fq2vcf_output/statistics
mkdir -p fq2vcf_output/vcf
mkdir -p fq2vcf_output/result

echo
echo -e $green"############################## Alignment ##############################"$none
echo 

for fq in $(ls $path_1 | grep -E ".fq$|.fastq$|.fq.gz$|.fastq.gz$")
do

		echo 
		echo -e $red"Processing "$fq$none
		echo 
		echo
		name=$( echo $fq | sed 's/.fq$//g;s/.fastq$//g;s/.fq.gz$//g;s/.fastq.gz$//g')
		bwa mem -t 8 ref/hg38.fa $path_1/$fq > fq2vcf_output/sam/$name.sam
done

echo 
echo -e $green"############################## post alignment ##############################"$none
echo 

for sam in $(ls fq2vcf_output/sam/*.sam)
do

	if [ "$sam | awk '{print $5}'" != 0 ]; then
		echo
		echo -e $red"Processing "$sam$none
		echo
		echo
		name=$( echo $sam | sed 's/fq2vcf_output\/sam\///g;s/.sam$//g')
		java -Xmx10g -jar tools/picard.jar SortSam  I=$sam O=fq2vcf_output/bam/$name.bam SO=coordinate
		java -Xmx10g -jar tools/picard.jar MarkDuplicates I=fq2vcf_output/bam/$name.bam O=fq2vcf_output/bam/$name.dedup.bam M=fq2vcf_output/statistics/$name.duplication_metrics.txt
		rm fq2vcf_output/bam/$name.bam
		java -Xmx10g -jar tools/picard.jar AddOrReplaceReadGroups  I=fq2vcf_output/bam/$name.dedup.bam O=fq2vcf_output/bam/$name.RG.bam ID=f1.l1 LB=li1 PL=ILLUMINA PU=fbi SM=$name
		rm fq2vcf_output/bam/$name.dedup.bam
		java -Xmx10g -jar tools/picard.jar BuildBamIndex I=fq2vcf_output/bam/$name.RG.bam O=fq2vcf_output/bam/$name.RG.bam.bai

	fi

done

echo 
echo -e $green"############################## Base Recalibration ##############################"$none
echo 

for bam in $(ls fq2vcf_output/bam/*.bam)
do

	if [ "$bam | awk '{print $5}'" != 0 ]; then
		echo
		echo -e $red"Processing "$bam$none
		echo
		echo
		name=$( echo $bam | sed 's/fq2vcf_output\/bam\///g;s/.RG.bam$//g')
		java -Xmx10g -jar tools/gatk.jar BaseRecalibrator -I $bam -R ref/hg38.fa --known-sites ref/hg38.dbsnp138.vcf -O fq2vcf_output/statistics/$name.recal_data.table
		java -Xmx10g -jar tools/gatk.jar ApplyBQSR -I $bam -R ref/hg38.fa --bqsr-recal-file fq2vcf_output/statistics/$name.recal_data.table -O fq2vcf_output/bam/$name.recalibrated.bam
		rm $bam $bam.bai
	fi
done


echo  
echo -e $green"############################## Variant calling ##############################"$none
echo 

for bam in $(ls fq2vcf_output/bam/*.bam)
do

	if [ "$bam | awk '{print $5}'" != 0 ]; then
		echo
		echo -e $red"Processing "$bam$none
		echo
		echo
		name=$( echo $bam | sed 's/fq2vcf_output\/bam\///g;s/.bam$//g;s/.recalibrated//g')
		java -Xmx10g -jar tools/gatk.jar HaplotypeCaller -I $bam -R ref/hg38.fa -D ref/hg38.dbsnp138.vcf -O fq2vcf_output/vcf/$name.vcf -A Coverage -A QualByDepth -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A ReadPosRankSumTest -A RMSMappingQuality
		
	fi

done
echo
echo -e $green"############################## Variant recalibration ##############################"$none
echo 

for vcf in $(ls fq2vcf_output/vcf/*.vcf)
do

	if [ "$vcf | awk '{print $5}'" != 0 ]; then
		echo
		echo -e $red"Processing "$vcf$none
		echo
		echo
		name=$( echo $vcf | sed 's/fq2vcf_output\/vcf\///g;s/.vcf$//g')
		java -Xmx10g -jar tools/gatk.jar VariantRecalibrator -V $vcf -R ref/hg38.fa --resource:omni,known=false,training=true,truth=false,prior=12.0 ref/1000G_omni2.5.hg38.vcf --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ref/hapmap_3.3.hg38.vcf --resource:dbsnp,known=true,training=true,truth=true,prior=2.0 ref/hg38.dbsnp138.vcf -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -tranche 99.95 -mode BOTH -O fq2vcf_output/statistics/$name.output.recal --tranches-file fq2vcf_output/statistics/$name.output.tranches  
		java -Xmx10g -jar tools/gatk.jar ApplyVQSR -V $vcf -R ref/hg38.fa -O fq2vcf_output/vcf/$name.recalibrated.vcf --truth-sensitivity-filter-level 99.0 --tranches-file fq2vcf_output/statistics/$name.output.tranches --recal-file fq2vcf_output/statistics/$name.output.recal -mode BOTH
		rm $vcf $vcf.idx
	fi

done
echo
echo -e $green"############################## Annotation ##############################"$none
echo 

for vcf in $(ls fq2vcf_output/vcf/*.vcf)
do

	if [ "$vcf | awk '{print $5}'" != 0 ]; then
		echo
		echo -e $red"Processing "$vcf$none
		echo
		echo
		name=$( echo $vcf | sed 's/fq2vcf_output\/vcf\///g;s/.recalibrated.vcf$//g')
		perl tools/annovar/convert2annovar.pl -format vcf4 $vcf > fq2vcf_output/vcf/$name.avinput
		perl tools/annovar/table_annovar.pl fq2vcf_output/vcf/$name.avinput tools/annovar/hg38 -buildver hg38 -out fq2vcf_output/result/$name -protocol refGene,knownGene,refGeneWithVer,ensGene,hrcr1,exac03,gnomad_exome,avsnp147,clinvar_20220320,intervar_20180118,dbnsfp30a -operation gx,g,g,g,f,f,f,f,f,f,f --nastring . --thread 10 --polish --remove --csvout -xref tools/annovar/example/gene_fullxref.txt
	fi

done

