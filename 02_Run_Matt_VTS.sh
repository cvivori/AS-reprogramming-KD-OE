#!/bin/bash
#See http://matt.crg.eu/

## GENERATE MATT-COMPATIBLE MATRIX FROM VAST-TOOLS COMBINE OUTPUT (all mapped events)
# matt get_vast ../vast_out/INCLUSION_LEVELS_FULL-Mm218.tab -gtf /no_backup/mirimia/genome_annots/ensembl/Mmu10.gtf > Matt_Sets/tmp_INCLUSION.txt


## CREATE TABLES WITH SELECTED EVENTS (CPSF3-, UL1- and TIA1-dependent or -independent)
# matt get_rows Matt_Sets/tmp_INCLUSION.txt EVENT]Sets_DependentEvents/INCL_CPSF3_dep.txt,EVENT[ | matt add_val - GROUP CPSF3dep > Table_AllEvents.txt
# matt get_rows Matt_Sets/tmp_INCLUSION.txt EVENT]Sets_DependentEvents/INCL_CPSF3_indep.txt,EVENT[ | matt add_val - GROUP CPSF3indep | matt add_rows Table_AllEvents.txt -
# matt get_rows Matt_Sets/tmp_INCLUSION.txt EVENT]Sets_DependentEvents/INCL_UL1_dep.txt,EVENT[ | matt add_val - GROUP UL1dep | matt add_rows Table_AllEvents.txt -
# matt get_rows Matt_Sets/tmp_INCLUSION.txt EVENT]Sets_DependentEvents/INCL_UL1_indep.txt,EVENT[ | matt add_val - GROUP UL1indep | matt add_rows Table_AllEvents.txt -
# matt get_rows Matt_Sets/tmp_INCLUSION.txt EVENT]Sets_DependentEvents/INCL_TIA1_dep.txt,EVENT[ | matt add_val - GROUP TIA1dep | matt add_rows Table_AllEvents.txt -
# matt get_rows Matt_Sets/tmp_INCLUSION.txt EVENT]Sets_DependentEvents/INCL_TIA1_indep.txt,EVENT[ | matt add_val - GROUP TIA1indep | matt add_rows Table_AllEvents.txt -


## CREATE CONTROL SETS FOR CPSF3/UL1 AND TIA1
# cat ../vast_out/DiffAS_dPSI10_files_with_dPSI/AS_NC-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl-* | cut -f 2 | sort -u > ../vast_out/DiffAS_dPSI10_files_with_dPSI/ALLuniqueEVENTS_AS_NC-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl.txt
# cat Sets_DependentEvents/INCL_CPSF3_dep*.txt Sets_DependentEvents/INCL_UL1_dep*.txt | sort -u > Sets_DependentEvents/KD_dep_unique.txt
# matt get_rows Matt_Sets/tmp_INCLUSION.txt EVENT]../vast_out/DiffAS_dPSI10_files_with_dPSI/ALLuniqueEVENTS_AS_NC-Mm218-dPSI10-range5-min_ALT_use25_KD_day00_Ctrl.txt,EVENT[ > Table_allASNC_KD.txt
# matt get_rows Table_allASNC_KD.txt '!EVENT]Sets_DependentEvents/KD_dep_unique.txt,EVENT[' > Table_ASNC_KD.txt

# cat ../vast_out/DiffAS_dPSI10_files_with_dPSI/AS_NC-Mm218-dPSI10-range5-min_ALT_use25_OE_day00_Ctrl-* | cut -f 2 | sort -u > ../vast_out/DiffAS_dPSI10_files_with_dPSI/ALLuniqueEVENTS_AS_NC-Mm218-dPSI10-range5-min_ALT_use25_OE_day00_Ctrl.txt
# cat Sets_DependentEvents/INCL_TIA1_dep*.txt | sort -u > Sets_DependentEvents/OE_dep_unique.txt
# matt get_rows Matt_Sets/tmp_INCLUSION.txt EVENT]../vast_out/DiffAS_dPSI10_files_with_dPSI/ALLuniqueEVENTS_AS_NC-Mm218-dPSI10-range5-min_ALT_use25_OE_day00_Ctrl.txt,EVENT[  > Table_allASNC_OE.txt
# matt get_rows Table_allASNC_OE.txt '!EVENT]Sets_DependentEvents/OE_dep_unique.txt,EVENT[' > Table_ASNC_OE.txt


## SEPARATE ALTERNATIVE EXONS
# matt get_rows Table_AllEvents.txt COMPLEX]S,C1,C2,C3,ANN,MIC[ > Table_AltEx.txt

## ADD CONTROL SETS TO ALTERNATIVE EXONS (similar size to the previously generated sets)
# matt col_uniq Table_AltEx.txt GROUP ## to get group sizes
# matt get_rows Table_allASNC_KD.txt COMPLEX]S,C1,C2,C3,ANN,MIC[ | matt rand_rows - 200 5678 | matt add_val - GROUP ASNC_KD | matt add_rows Table_AltEx.txt -
# matt get_rows Table_allASNC_OE.txt COMPLEX]S,C1,C2,C3,ANN,MIC[ | matt rand_rows - 250 5678 | matt add_val - GROUP ASNC_OE | matt add_rows Table_AltEx.txt -


## RUN MATT cmpr_exons on ALTERNATIVE EXONS
# qsub_job 12 8 8 -N MceC matt cmpr_exons Table_AltEx.txt START END SCAFFOLD STRAND GENEID /no_backup/mirimia/genome_annots/ensembl/Mmu10.gtf /no_backup/mirimia/genome_seqs/Mmu10_gDNA.fasta Mmul 150 GROUP[CPSF3dep,CPSF3indep,UL1dep,UL1indep,ASNC_KD] MATT_Analysis/OUT_cmprexons_CEx_KD
# qsub_job 12 8 8 -N MceC matt cmpr_exons Table_AltEx.txt START END SCAFFOLD STRAND GENEID /no_backup/mirimia/genome_annots/ensembl/Mmu10.gtf /no_backup/mirimia/genome_seqs/Mmu10_gDNA.fasta Mmul 150 GROUP[CPSF3dep,CPSF3indep,ASNC_KD] MATT_Analysis/OUT_cmprexons_CEx_C
# qsub_job 12 8 8 -N MceU matt cmpr_exons Table_AltEx.txt START END SCAFFOLD STRAND GENEID /no_backup/mirimia/genome_annots/ensembl/Mmu10.gtf /no_backup/mirimia/genome_seqs/Mmu10_gDNA.fasta Mmul 150 GROUP[UL1dep,UL1indep,ASNC_KD] MATT_Analysis/OUT_cmprexons_CEx_U
# qsub_job 12 8 8 -N MceT matt cmpr_exons Table_AltEx.txt START END SCAFFOLD STRAND GENEID /no_backup/mirimia/genome_annots/ensembl/Mmu10.gtf /no_backup/mirimia/genome_seqs/Mmu10_gDNA.fasta Mmul 150 GROUP[TIA1dep,TIA1indep,ASNC_OE] MATT_Analysis/OUT_cmprexons_CEx_T

## RUN MATT RNAmaps CISBP
# qsub_job 24 12 1 -N Mrnam_T matt rna_maps_cisbp Table_AltEx.txt UPSTRM_EX_BORDER START END DOSTRM_EX_BORDER SCAFFOLD STRAND GROUP[TIA1dep,ASNC_OE,TIA1indep] 31 50 150 /no_backup/mirimia/genome_seqs/Mmu10_gDNA.fasta cisbprna_regexps -d MATT_Analysis/OUT_rnamaps_cisbp_CEx_T

## RUN MATT RNAmaps CISBP on TIA1 motifs with pvalue"
# qsub_job 24 12 1 -N Mrnam_T matt rna_maps_cisbp Table_AltEx.txt UPSTRM_EX_BORDER START END DOSTRM_EX_BORDER SCAFFOLD STRAND GROUP[TIA1dep,TIA1indep,ASNC_OE] 31 50 150 /no_backup/mirimia/genome_seqs/Mmu10_gDNA.fasta cisbprna_regexps -names TIA1 -d MATT_Analysis/OUT_rnamaps_cisbp_CEx_T_pvalT001_ASNC -p 0.001 10000


