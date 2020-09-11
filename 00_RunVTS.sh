#!/bin/bash
# See https://github.com/vastgroup/vast-tools

## VASTTOOLS ALIGN, mm10
# path="/users/jvalcarcel/cvivori/RNAseq_MEFs2iPS_shCPSF3_shUL1_T7TIA1/FASTQ_Links";
# for f in $path/*_read1.fastq.gz; 
# do echo "
# Processing $f.."; 
#    		filename="${f##*/}";
#    		sample="$(echo "$filename" | sed 's/_read1.fastq.gz//g')"; # remove _read1.fastq.gz
# 			f1="${sample}_read1.fastq.gz";
# 			f2="${sample}_read2.fastq.gz";
# 		qsub_job 48 20 15 -N vt222A vast-tools-2.2.2 align $path/$f1 $path/$f2 -sp Mm2 -c 15 -o vast_out --IR_version 2
# done


## RENAME OUTPUTS TO COMBINE
# rename "_read1." "." vast_out/to_combine/*

## no VASTOOLS MERGE

## VASTTOOLS COMBINE
# qsub_job 12 8 2 -N vtC vast-tools-2.2.2 combine -sp Mm2


## VASTTOOLS COMPARE dPSI0 (!)
# qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day00_Ctrl_1,KD_day00_Ctrl_2 -b KD_day12_shSCR_1,KD_day12_shSCR_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
# qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day00_Ctrl_1,KD_day00_Ctrl_2 -b KD_day12_shCPSF3_n1_1,KD_day12_shCPSF3_n1_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
# qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day00_Ctrl_1,KD_day00_Ctrl_2 -b KD_day12_shCPSF3_n5_1,KD_day12_shCPSF3_n5_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
# qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day00_Ctrl_1,KD_day00_Ctrl_2 -b KD_day12_shUL1_n3_1,KD_day12_shUL1_n3_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
# qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day00_Ctrl_1,KD_day00_Ctrl_2 -b KD_day12_shUL1_n4_1,KD_day12_shUL1_n4_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
# qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day12_shSCR_1,KD_day12_shSCR_2 -b KD_day12_shCPSF3_n1_1,KD_day12_shCPSF3_n1_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
# qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day12_shSCR_1,KD_day12_shSCR_2 -b KD_day12_shCPSF3_n5_1,KD_day12_shCPSF3_n5_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
# qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day12_shSCR_1,KD_day12_shSCR_2 -b KD_day12_shUL1_n3_1,KD_day12_shUL1_n3_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
# qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day12_shSCR_1,KD_day12_shSCR_2 -b KD_day12_shUL1_n4_1,KD_day12_shUL1_n4_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets

# qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a OE_day00_Ctrl_1,OE_day00_Ctrl_2 -b OE_day12_Empty_1,OE_day12_Empty_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
# qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a OE_day00_Ctrl_1,OE_day00_Ctrl_2 -b OE_day12_T7TIA1_1,OE_day12_T7TIA1_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
# qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a OE_day12_Empty_1,OE_day12_Empty_2 -b OE_day12_T7TIA1_1,OE_day12_T7TIA1_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets




## VASTTOOLS COMPARE dPSI10
#qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day00_Ctrl_1,KD_day00_Ctrl_2 -b KD_day12_shSCR_1,KD_day12_shSCR_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
#qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day00_Ctrl_1,KD_day00_Ctrl_2 -b KD_day12_shCPSF3_n1_1,KD_day12_shCPSF3_n1_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
#qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day00_Ctrl_1,KD_day00_Ctrl_2 -b KD_day12_shCPSF3_n5_1,KD_day12_shCPSF3_n5_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
#qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day00_Ctrl_1,KD_day00_Ctrl_2 -b KD_day12_shUL1_n3_1,KD_day12_shUL1_n3_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
#qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day00_Ctrl_1,KD_day00_Ctrl_2 -b KD_day12_shUL1_n4_1,KD_day12_shUL1_n4_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
#qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day12_shSCR_1,KD_day12_shSCR_2 -b KD_day12_shCPSF3_n1_1,KD_day12_shCPSF3_n1_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
#qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day12_shSCR_1,KD_day12_shSCR_2 -b KD_day12_shCPSF3_n5_1,KD_day12_shCPSF3_n5_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
#qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day12_shSCR_1,KD_day12_shSCR_2 -b KD_day12_shUL1_n3_1,KD_day12_shUL1_n3_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
#qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a KD_day12_shSCR_1,KD_day12_shSCR_2 -b KD_day12_shUL1_n4_1,KD_day12_shUL1_n4_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets

#qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a OE_day00_Ctrl_1,OE_day00_Ctrl_2 -b OE_day12_Empty_1,OE_day12_Empty_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
#qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a OE_day00_Ctrl_1,OE_day00_Ctrl_2 -b OE_day12_T7TIA1_1,OE_day12_T7TIA1_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
#qsub_job 6 6 6 -N vtCMP vast-tools-2.2.2 compare vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -a OE_day12_Empty_1,OE_day12_Empty_2 -b OE_day12_T7TIA1_1,OE_day12_T7TIA1_2 --min_dPSI 0 --min_range 0  --GO -sp Mm2 --print_dPSI --print_sets
