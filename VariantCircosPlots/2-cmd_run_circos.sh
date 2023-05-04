#!/bin/bash
# you need to have circos installed to run this bash file

circos -conf heatmap_circos_variants.conf -debug_group legend | tee out_variants_withMax_no_het.txt

circos -conf heatmap_circos_snp_STRAINS.conf -debug_group legend | tee out_snp_strain_no_het.txt

circos -conf heatmap_circos_indel_STRAINS.conf -debug_group legend | tee out_indel_strain_no_het.txt

circos -conf heatmap_circos_sv_STRAINS.conf -debug_group legend | tee out_sv_strain.txt
