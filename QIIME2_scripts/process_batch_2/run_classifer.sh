cd /scratch/mit_lieberman/projects/alex_corynebacterium/7_30_21_Neel_sample_batch/
source activate /scratch/mit_lieberman/projects/alex_corynebacterium/QIIME2_new_download/qiime2-2020.1

echo "BEGIN CLASSIFICATION FILE"
echo "CREATED CLASSIFER - BEGINNING CLASSIFICATION"

qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads /scratch/mit_lieberman/projects/alex_corynebacterium/7_30_21_Neel_sample_batch/data.demux_se.trim.sq4.dada2.qza \
  --o-classification taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.dada2.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.dada2.qza \
  --o-visualization taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.dada2.qzv

echo "END CLASSIFCATION"
echo "BEGIN EXPORT"

# get data
#qiime tools export \
#  --input-path taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.dada2.bar-plots.qzv \
#  --output-path taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.dada2.bar-plots_export

# In order to see all your raw OTUs, their classification, classification-confidence-level, and sample counts, the below code can be used. 
# feature-table.tsv will hold all your sample counts with QIIME2 generated OTU-labels (this will look like a long string of random numbers and letters)
# The directory taxonomic_OTU_classification will hold metadata.tsv, which will have QIIME2's OTU label, classification, and classification confidence value
# The directory SILVA_OTU_sequences will contain sequences.fasta that holds your OTU fasta sequences and their QIIME2 label
# All these files can be combined manually to make one descriptive CSV. I've found that QIIME2 has outputted all OTUs in the same order every time, which makes
# it easy to copy and paste everything into one big spreadsheet. I would check this, however, by copy and pasting the OTUs from different files
# into a string array in MATLAB and comparing them. 

# Note: QIIME2 includes a whole range of packages for diversity analysis, phylogenetics, etc. There's a lot more you can do with QIIME2 beyond this script! 

qiime tools export \
  --input-path /scratch/mit_lieberman/projects/alex_corynebacterium/7_30_21_Neel_sample_batch/data.demux_se.trim.sq4.dada2.qzv \
  --output-path SILVA_OTU_sequences
	
qiime tools export \
  --input-path taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.dada2.qzv \
  --output-path taxonomic_OTU_classification
	
qiime tools export \
  --input-path /scratch/mit_lieberman/projects/alex_corynebacterium/7_30_21_Neel_sample_batch/dada2-table.qza \
  --output-path dada2_OTU_export
  
biom convert -i dada2_OTU_export/feature-table.biom -o dada2_OTU_export/feature-table.tsv --to-tsv

echo "TAXA ASSIGNMENT DATA EXPORTED!!!"

