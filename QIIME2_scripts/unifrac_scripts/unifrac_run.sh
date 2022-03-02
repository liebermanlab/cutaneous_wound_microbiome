cd  /scratch/mit_lieberman/projects/alex_corynebacterium/Unifrac_final_run
source activate /scratch/mit_lieberman/projects/alex_corynebacterium/QIIME2_new_download/qiime2-2020.1

echo "Import files"

biom convert -i 07-Feb-2022_AVSs_combined_batches.txt -o cleaned_ASV_table_for_unifrac_hdf5.biom --table-type="OTU table" --to-hdf5

qiime tools import \
  --input-path cleaned_ASV_table_for_unifrac_hdf5.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path feature-table.qza


qiime tools import \
  --input-path unifrac_sequences.txt \
  --output-path /scratch/mit_lieberman/projects/alex_corynebacterium/Unifrac_final_run/sequences.qza \
  --type FeatureData[Sequence]


echo "Make tree"


qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences sequences.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza \
  --verbose

echo "Make unifrac table"

qiime diversity-lib weighted-unifrac \
  --i-table feature-table.qza\
  --i-phylogeny rooted-tree.qza  \
  --p-threads 8 \
  --o-distance-matrix distance_matrix.qza \
  --verbose

  
echo "Export unifrac table"

qiime tools export \
  --input-path distance_matrix.qza \
  --output-path table_output/unifrac_exported_3.txt

echo "TAXA ASSIGNMENT DATA EXPORTED!!!"

