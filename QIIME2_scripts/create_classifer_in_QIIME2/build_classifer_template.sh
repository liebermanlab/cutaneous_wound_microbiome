cd /scratch/mit_lieberman/projects/alex_corynebacterium/Build_classifier/
source activate /scratch/mit_lieberman/projects/alex_aro_cacnes/SATIVA_QIIME_cutibacterium_training/QIIME_analysis/qiime2-2020.1

echo "BEGIN CLASSIFICATION FILE"

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path /scratch/mit_lieberman/projects/alex_corynebacterium/Build_classifier/classifer_refseqs/classifer_input_seqs.fna \
  --output-path ref-seqs.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path /scratch/mit_lieberman/projects/alex_corynebacterium/Build_classifier/classifer_refseqs/classifer_input_seqs.tax \
  --output-path ref-taxonomy.qza
 
 echo "IMPORTED DATA"

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier classifier.qza