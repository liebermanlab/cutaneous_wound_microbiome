cd /scratch/mit_lieberman/projects/alex_corynebacterium/Neel_3_14_21/
source activate /scratch/mit_lieberman/projects/alex_aro_cacnes/SATIVA_QIIME_cutibacterium_training/QIIME_analysis/qiime2-2020.1

echo "BEGIN CLASSIFICATION FILE"

qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads /scratch/mit_lieberman/projects/alex_corynebacterium/Neel_3_14_21/data.demux_se.trim.sq4.dada2.qza \
  --o-classification taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.dada2.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.dada2.qza \
  --o-visualization taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.dada2.qzv

echo "END CLASSIFCATION"
echo "BEGIN EXPORT"

qiime tools export \
  --input-path /scratch/mit_lieberman/projects/alex_corynebacterium/Neel_3_14_21/data.demux_se.trim.sq4.dada2.qzv \
  --output-path SILVA_OTU_sequences
	
qiime tools export \
  --input-path taxonomy.qiimeres-classifier.data.demux_se.trim.sq4.dada2.qzv \
  --output-path taxonomic_OTU_classification
	
qiime tools export \
  --input-path /scratch/mit_lieberman/projects/alex_corynebacterium/Neel_3_14_21/dada2-table.qza \
  --output-path dada2_OTU_export
  
biom convert -i dada2_OTU_export/feature-table.biom -o dada2_OTU_export/feature-table.tsv --to-tsv

echo "TAXA ASSIGNMENT DATA EXPORTED!!!"

