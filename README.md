# Cutaneous surgical wounds have distinct microbiomes from intact skin
This repository contains the code needed to reproduce figures for for Gupta &amp; Poret et al.'s "Cutaneous surgical wounds have distinct microbiomes from intact skin" (currently in review, raw data can be found at: PRJNA809947). A description of each file uploaded will be included below. 

`SG_AJP_wound_microbiome_figure_creation.m` is the MATLAB script used to generate all figures used in this paper. Associated function files needed to run this script are stored in `\matlab_functions`. This script intakes `batch_1_2_QIIME2_output.xlsx`, `batch_3_QIIME2_output.xlsx`, `clinical_metadata.xlsx`, `distance-matrix_final.xlsx`, and `patient_sequencing_metadata.xlsx` and outputs both images and files that will be saved in a subdirectory `/saved_excel_files`. The excel files are saved in the format `{date}_OTUs_combined_batches.xlsx`, `{date}_all_unique_taxons.xlsx`, `{date}_species_counts.xlsx`, `{date}_unique_genus_counts.xlsx`, and contain raw OTU counts, counts combined by taxon (ex. all *Cutibacterium acnes* labeled taxons are combined in this sheet), counts grouped by genus, and counts grouped by species. All these excel sheets contain data post-quality-control filtering, and the {date} refers to the date of file generation. 

Raw data from QIIME2 is outputted in `batch_1_2_QIIME2_output.xlsx` and `batch_3_QIIME2_output.xlsx`. Clinical metadata (i.e. tumor type, age, sex, etc.) is stored in `clinical_metadata.xlsx` and sequencing metadata is stored in `patient_sequencing_metadata.xlsx`. A unifrac distance matrix generated from filtered all OTUs present in `{date}_OTUs_combined_batches` is written to `distance-matrix_final.xlsx`.
