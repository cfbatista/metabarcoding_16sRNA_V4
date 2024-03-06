
# IMPORTING DATA
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path /mnt/c/Users/camil/Documents/SMA/dada/raw_data \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path demux.qza

--- Foward ou Reverse reads
qiime tools import \
--type 'SampleData[SequencesWithQuality]' \
--input-path /mnt/c/Users/camil/Documents/diabiose_frango_corte/foward_reads/raw_data \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path demux.qza

qiime demux summarize \
--i-data demux.qza \
--o-visualization demux.qzv

#CUTADAPT
qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux.qza \
--p-front-f GTGYCAGCMGCCGCGGTAA \
--p-front-r GGACTACNVGGGTWTCTAAT \
--p-no-indels \
--p-error-rate 0 \
--o-trimmed-sequences demux_trim.qza \
--verbose

qiime demux summarize \
--i-data demux_trim.qza \
--o-visualization demux_trim.qzv

# DEBLUR
qiime quality-filter q-score \
 --i-demux demux.qza \
 --o-filtered-sequences demux-filtered.qza \
 --o-filter-stats demux-filter-stats.qza

qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-filtered.qza \
  --p-trim-length 250 \
  --o-representative-sequences rep-seqs-deblur.qza \
  --o-table table-deblur.qza \
  --p-sample-stats \
  --o-stats deblur-stats.qza


qiime metadata tabulate \
  --m-input-file demux-filter-stats.qza \
  --o-visualization demux-filter-stats.qzv

qiime deblur visualize-stats \
  --i-deblur-stats deblur-stats.qza \
  --o-visualization deblur-stats.qzv


# DADA
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 250 \
  --p-trunc-len-r 250 \
  --p-trunc-q 15 \
  --p-n-threads 2 \
  --o-table table-dada2.qza \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-denoising-stats denoising-stats_dada.qza \
  --verbose

qiime feature-table summarize \
--i-table table-deblur.qza \
--o-visualization table-dada2.qzv \
--m-sample-metadata-file /mnt/c/Users/camil/Documents/SMA/dada/sample-metadata.tsv

# PHYLOGENETIC TREE
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences rep-seqs-deblur.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

# ALPHA RAREFACTION
qiime diversity alpha-rarefaction \
--i-table table-deblur.qza \
--i-phylogeny rooted-tree.qza \
--p-max-depth 684 \
--m-metadata-file sample-metadata.tsv \
--o-visualization alpha-rarefaction.qzv

# TAXONOMY CLASSIFICATION
qiime feature-classifier classify-sklearn \
--i-classifier /mnt/c/Users/camil/Documents/TCC/covid/teste/silva-138-99-nb-classifier.qza \
--i-reads rep-seqs-dada2.qza \
--o-classification taxonomy-silva.qza \
--p-n-jobs 1

qiime metadata tabulate \
--m-input-file taxonomy-silva.qza \
--o-visualization taxonomy-deion.qzv


qiime taxa barplot \
--i-table table-deblur.qza \
--i-taxonomy taxonomy-gg-99.qza \
--m-metadata-file sample-metadata.tsv \
--o-visualization taxa-bar-plots-gg.qzv


qiime feature-table filter-features \
--i-table table-deblur.qza \
--p-min-samples 2 \
--p-min-frequency 10 \
--o-filtered-table filtered_table_ancom.qza

qiime composition add-pseudocount \
--i-table filtered_table_ancom.qza \
--o-composition-table comp-table.qza

qiime composition ancom \
--i-table comp-table.qza \
--m-metadata-file sample-metadata.tsv \
--m-metadata-column Group \
--p-transform-function 'clr' \
--o-visualization ancom-group.qzv

qiime taxa collapse \
--i-table filtered_table_ancom.qza \
--i-taxonomy taxonomy-gg-99.qza \
--p-level 5 \
--o-collapsed-table l5-filtered-table.qza

qiime composition add-pseudocount \
--i-table l5-filtered-table.qza \
--o-composition-table comp-table-l5.qza

qiime composition ancom \
--i-table comp-table-l5.qza \
--m-metadata-file sample-metadata.tsv \
--m-metadata-column Group \
--p-transform-function 'clr' \
--o-visualization l5-ancom-group.qzv


qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted-tree.qza \
--i-table table-deblur.qza \
--p-sampling-depth 1011 \
--m-metadata-file sample-metadata.tsv \
--output-dir core-metrics-results-depth

conda install -c conda-forge deicode

qiime deicode rpca \
--i-table table-deblur.qza \
--p-min-feature-count 10 \
--p-min-sample-count 500 \
--p-max-iterations 30 \
--o-biplot ordination.qza \
--o-distance-matrix distance.qza

qiime emperor biplot \
--i-biplot ordination.qza \
--m-sample-metadata-file sample-metadata.tsv \
--m-feature-metadata-file taxonomy-gg-99.qza \
--o-visualization biplot.qzv \
--p-number-of-features 10

qiime diversity beta-group-significance \
--i-distance-matrix distance.qza \
--m-metadata-file sample-metadata.tsv \
--m-metadata-column Group \
--p-method permanova \
--o-visualization region_aitchison_significance-group.qzv

qiime gemelli phylogenetic-rpca-with-taxonomy \
--i-table table-deblur.qza \
--i-phylogeny rooted-tree.qza \
--m-taxonomy-file taxonomy-gg-99.qza \
--p-min-feature-count 10 \
--p-min-sample-count 500 \
--o-biplot ordination.qza \
--o-distance-matrix distance.qza \
--o-counts-by-node-tree phylo-tree.qza \
--o-counts-by-node phylo-table.qza \
--o-t2t-taxonomy phylo-taxonomy.qza


!qiime empress community-plot \
--i-tree phylo-tree.qza \
--i-feature-table phylo-table.qza \
--i-pcoa ordination.qza \
--m-sample-metadata-file sample-metadata.tsv \
--m-feature-metadata-file phylo-taxonomy.qza \
--p-filter-missing-features \
--p-number-of-features 10 \
--o-visualization phylo-empress.qzv

qiime diversity beta-group-significance \
    --i-distance-matrix distance.qza \
    --m-metadata-file sample-metadata.tsv \
    --m-metadata-column Group \
    --p-method permanova \
    --o-visualization phylo_significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-depth/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results-depth/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-depth/evenness_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results-depth/evenness-group-significance.qzv