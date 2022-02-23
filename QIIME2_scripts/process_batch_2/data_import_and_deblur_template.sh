cd /scratch/mit_lieberman/projects/alex_corynebacterium/7_30_21_Neel_sample_batch/
source activate /scratch/mit_lieberman/projects/alex_corynebacterium/QIIME2_new_download/qiime2-2020.1
#source activate /scratch/mit_lieberman/projects/alex_aro_cacnes/SATIVA_QIIME_cutibacterium_training/QIIME_analysis/qiime2-2020.1

echo "INITIALIZED"
# IF DEMULTIPLEXED
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path /scratch/mit_lieberman/projects/alex_corynebacterium/7_30_21_Neel_sample_batch/manifest_demux_file_loc.tsv \
  --output-path data.demux_se.qza \
  --input-format SingleEndFastqManifestPhred33V2

echo "IMPORTED"


# check what we got
qiime demux summarize \
  --i-data data.demux_se.qza \
  --o-visualization data.demux_se.qzv

echo "CUTADAPT"
# remove primer (use only SE data here, thus need to remove only one)
qiime cutadapt trim-single \
  --p-cores 8 \
  --i-demultiplexed-sequences data.demux_se.qza \
  --p-front AGAGTTTGATCMTGGCTCAG \
  --o-trimmed-sequences data.demux_se.trim.qza

# check what we got
qiime demux summarize \
  --i-data data.demux_se.trim.qza \
  --o-visualization data.demux_se.trim.qzv


##################
###### Quality trim and denoise using dada2
##################

echo "TRIM AND DENOISE"

# remove low seq quality reads
qiime quality-filter q-score \
  --i-demux data.demux_se.trim.qza \
  --p-min-quality 4 \
  --o-filtered-sequences data.demux_se.trim.sq4.qza \
  --o-filter-stats filter_stats_sq4.qza

 # check what we got
qiime demux summarize \
  --i-data data.demux_se.trim.sq4.qza \
  --o-visualization data.demux_se.trim.sq4.qzv

 
qiime dada2 denoise-single \
  --i-demultiplexed-seqs data.demux_se.trim.sq4.qza \
  --p-trunc-len 180 \
  --p-n-threads 8 \
  --o-representative-sequences data.demux_se.trim.sq4.dada2.qza \
  --o-table dada2-table.qza \
  --o-denoising-stats dada2-stats.qza

# visualize results 
qiime metadata tabulate \
  --m-input-file filter_stats_sq4.qza \
  --o-visualization filter_stats_sq4.qzv

#qiime dada2 visualize-stats \
#  --i-dada2-stats dada2-stats.qza \
#  --o-visualization dada2-stats.qzv

# FeatureTable and FeatureData summaries
qiime feature-table summarize \
  --i-table dada2-table.qza \
  --o-visualization dada2-table.qzv \
  --m-sample-metadata-file metadata_sample_desc.tsv 

qiime feature-table tabulate-seqs \
  --i-data data.demux_se.trim.sq4.dada2.qza \
  --o-visualization data.demux_se.trim.sq4.dada2.qzv
 