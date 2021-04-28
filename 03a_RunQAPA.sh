#!/bin/bash

## Download and installed Qapa v1.3.0
# see https://github.com/morrislab/qapa

## Download annotation file for v1.3.0
# wget https://github.com/morrislab/qapa/releases/download/v1.3.0/qapa_3utrs.gencode_VM22.mm10.bed.gz

## Download gene annotations:
## Ensembl gene metadata table from Biomart.
# mysql --user=anonymous --host=martdb.ensembl.org --port=5316 -A ensembl_mart_89  -e "select stable_id_1023 as 'Gene stable ID', stable_id_1066 as 'Transcript stable ID',  biotype_1020 as 'Gene type', biotype_1064 as 'Transcript type',  display_label_1074 as 'Gene name' from mmusculus_gene_ensembl__transcript__main"   > ensembl_identifiers.txt
## GENCODE gene prediction annotation table
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select * from wgEncodeGencodeBasicVM9" mm10 > gencode.basic.txt

## Download polyA annotations:
## Standard approach using PolyASite and GENCODE poly(A) track (as described in Ha et al., Genome Biology 2018)
# 1) bedfile from PolyASite database http://polyasite.unibas.ch/ (mm10) 
## filtered as in the paper in clusters.mm10.bed (score in 8th column >= 3, annotated TE or DS in 10th column, see https://polyasite.unibas.ch/atlas)
# awk '$8 >= 3 && $10 ~ /TE|DS/' atlas.clusters.2.0.GRCm38.96.bed > clusters.mm10.bed
# 2) GENCODE Download from UCSC Table Browser and reordered as in the following command:
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A  -e "select chrom, txStart, txEnd, name2, score, strand  from wgEncodeGencodePolyaVM9 where name2 = 'polyA_site'" -N mm10 > gencode.polyA_sites.bed


## QAPA - extract 3′ UTRs from annotation
qapa build --db ensembl_identifiers.txt -g gencode.polyA_sites.bed -p clusters.mm10.bed gencode.basic.txt > output_utrs.bed
qapa fasta -f ~/References/mm10/UCSC/mm10.fa output_utrs.bed output_sequences.fa

## Build Salmon index
salmon index -t output_sequences.fa -i utr_library

## Salmon transcript quantification
path="/users/jvalcarcel/cvivori/RNAseq_MEFs2iPS_shCPSF3_shUL1_T7TIA1/FASTQ_Links/";
for f in $path/*read1.fastq.gz; 
do echo "
  Processing $f.."; 
    		filename="${f##*/}";
    		sample="$(cut -d'_' -f1-4 <<<"$filename")";
 			f1="${sample}_read1.fastq.gz";
 			f2="${sample}_read2.fastq.gz";
 		echo "$f1 	$f2";
 		#printf "$sample\n" >> sample_names.txt;
		qsub -q long-sl7 -V -cwd -N salm -pe ompi 8 -e logfiles/ -l virtual_free=24G -l h_rt=48:00:00 -b y salmon quant -i utr_library -l A -1 $path/$f1 -2 $path/$f2 -p 8 -o quants/${sample};
done

## QAPA - quantify 3′ UTR isoform usage
qapa quant --db ensembl_identifiers.txt Salmon/quants/*/quant.sf > pau_results.txt


