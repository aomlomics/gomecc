## This worksheet describes the bioinformatics processing of 16S and 18S amplicon data for the 2021 GOMECC-4 cruise in the Gulf of Mexico

# Activate QIIME 2
conda activate qiime2-2021.2
source tab-qiime

# Prepare taxonomic classifiers for 18S and 16S rRNA sequences
# Import newest PR2 taxonomic files for 18S
qime tools import --type 'FeatureData[Sequence]' --input-path pr2_version_5.0.0_SSU_mothur.fasta --output-path pr2_seqs.qza
qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path pr2_version_5.0.0_SSU_mothur.tax --output-path pr2_tax.qza

# Extract 18S primer to the V9 region
qiime feature-classifier extract-reads --i-sequences pr2_seqs_new.qza --p-f-primer GTACACACCGCCCGTC --p-r-primer TGATCCTTCTGCAGGTTCACCTAC --p-min-length 50 --p-max-length 400 --verbose  --o-reads refseqs_V9.qza 

# Build classifier for 18S 
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads refseqs_V9.qza --i-reference-taxonomy pr2_tax_new.qza --o-classifier classifier_V9.qza 

# Artifact files for 16S SILVA 138 are available online: https://docs.qiime2.org/2023.2/data-resources/#silva-16s-18s-rrna
# Extract 16S primers to the V4-V5 region
qiime feature-classifier extract-reads --i-sequences silva-138-99-seqs.qza --p-f-primer GTGYCAGCMGCCGCGGTAA --p-r-primer CCGYCAATTYMTTTRAGTTT --p-min-length 50 --p-max-length 500 --verbose  --o-reads refseqs_V4-V5.qza

# Build classifier for 16S
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads refseqs_V4-V5.qza --i-reference-taxonomy silva-138-99-tax.qza --o-classifier classifier_V4-V5.qza

# Navigate to folder - 18S first
cd Gomecc-4
cd 18S

# Import paired reads and visualize profiles
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path paired_reads --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end-18S.qza
qiime demux summarize --i-data demux-paired-end-18S.qza --o-visualization demux-paired-end-18S.qzv
qiime tools view demux-paired-end-16S.qzv

# Trim 18S V9 primers with Cutadapt, view trimmed reads, and export for Tourmaline
qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end-18S.qza --p-adapter-f GTACACACCGCCCGTC...GTAGGTGAACCTGCAGAAGGATCA --p-adapter-r TGATCCTTCTGCAGGTTCACCTAC...GACGGGCGGTGTGTAC --verbose --p-error-rate 0.5 --o-trimmed-sequences demux-trimmed-18S.qza
qiime demux summarize --i-data demux-trimmed-18S.qza --o-visualization demux-trimmed-18S.qzv
qiime tools view demux-trimmed-18S.qzv
qiime tools export --input-path demux-trimmed-18S.qza --output-path 18S-trimmed-reads

# Clone new version of Tourmaline and name folder
git clone https://github.com/aomlomics/tourmaline.git
mv tourmaline tourmaline-18S-09222022

# Create manifest files from the FASTQ files - script provided in Tourmaline
scripts/create_manifest_from_fastq_directory.py /Users/sean.r.anderson/Gomecc-4/tourmaline-18S-09222022/00-data/fastq 00-data/manifest_pe.csv 00-data/manifest_se.csv

# Execute Tourmaline - DADA2 paired end (parameters in config.yaml file)
snakemake dada2_pe_denoise

# Check number of reads that made it through each DADA2 step
cd 02-output-dada2-pe-unfiltered/00-table-repseqs/       
qiime metadata tabulate --m-input-file dada2_stats.qza --o-visualization dada2_stats.qzv
qiime tools view dada2_stats.qzv    

# Exectute Tourmaline to infer taxnomy and diversity
# For taxonomy, copy and paste above 18S V9 classifier into the 01-imported folder and rename to "classifier.qza"
snakemake dada2_pe_taxonomy_unfiltered
snakemake dada2_pe_diversity_unfiltered

## Proceed to 16S sequences
## 16S was sequenced on two separate MiSeq runs (plates 1-2 and 3-6)
# Navigate to 16S plates 1-2 folder
cd ..
cd 16S_plates1-2  

# Import paired reads and visualize profiles
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path paired_reads --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end-16S.qza
qiime demux summarize --i-data demux-paired-end-16S.qza --o-visualization demux-paired-end-16S.qzv
qiime tools view demux-paired-end-16S.qzv

# Trim 16S V4-V5 plate 1-2 primers with Cutadapt, view trimmed reads, and export for Tourmaline
qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end-16S.qza --p-front-f GTGYCAGCMGCCGCGGTAA --p-adapter-f AAACTYAAAKRAATTGRCGG --p-front-r CCGYCAATTYMTTTRAGTTT --p-adapter-r TTACCGCGGCKGCTGRCAC --verbose --p-error-rate 0.5  --o-trimmed-sequences demux-trimmed-16S.qza
qiime demux summarize --i-data demux-trimmed-16S.qza --o-visualization demux-trimmed-16S.qzv
qiime tools view demux-trimmed-16S.qzv
qiime tools export --input-path demux-trimmed-16S.qza --output-path 16S-trimmed-reads

# Navigate to Tourmaline folder
cd tourmaline-16S-09262022

# Create manifest files from the FASTQ files - script provided in Tourmaline
scripts/create_manifest_from_fastq_directory.py /Users/seananderson/Gomecc-4/16S_plates1-2/tourmaline-16S-09262022/00-data/fastq 00-data/manifest_pe.csv 00-data/manifest_se.csv

# Execute DADA2 paired end in Tourmaline and view number of reads that made it through
snakemake dada2_pe_denoise
cd 02-output-dada2-pe-unfiltered/00-table-repseqs/
qiime metadata tabulate --m-input-file dada2_stats.qza --o-visualization dada2_stats.qzv
qiime tools view dada2_stats.qzv

# Proceed to 16S plates 3-6
cd 16S_plates3-6

# Import paired reads and visualize profiles
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path paired_reads --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end-16S.qza
qiime demux summarize --i-data demux-paired-end-16S.qza --o-visualization demux-paired-end-16S.qzv
qiime tools view demux-paired-end-16S.qzv

# Trim 16S V4-V5 plate 3-6 primers with Cutadapt, view trimmed reads, and export for Tourmaline
qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end-16S.qza --p-front-f GTGYCAGCMGCCGCGGTAA --p-adapter-f AAACTYAAAKRAATTGRCGG --p-front-r CCGYCAATTYMTTTRAGTTT --p-adapter-r TTACCGCGGCKGCTGRCAC --verbose --p-error-rate 0.5  --o-trimmed-sequences demux-trimmed-16S.qza
qiime demux summarize --i-data demux-trimmed-16S.qza --o-visualization demux-trimmed-16S.qzv
qiime tools view demux-trimmed-16S.qzv
qiime tools export --input-path demux-trimmed-16S.qza --output-path 16S-trimmed-reads

# Navigate to Tourmaline folder
cd tourmaline-16S-09232022

# Create manifest files from the FASTQ files - script provided in Tourmaline
scripts/create_manifest_from_fastq_directory.py /Users/seananderson/Gomecc-4/16S_plates3-6/tourmaline-16S-09232022/00-data/fastq 00-data/manifest_pe.csv 00-data/manifest_se.csv

# Execute DADA2 paired end in Tourmaline and view number of reads that made it through
snakemake dada2_pe_denoise
cd 02-output-dada2-pe-unfiltered/00-table-repseqs/    
qiime metadata tabulate --m-input-file dada2_stats.qza --o-visualization dada2_stats.qzv
qiime tools view dada2_stats.qzv

## FASTQ files from both 16S sequencing runs were processed in DADA2 using the same parameters, and so, we can merge the ASV tables and infer taxonomy on the combined table. The idea is that some ASVs identified in plates 1-2 will also be in 3-6, mitigating biases and inflation of ASVs. 

# Merge 16S ASV and representative sequence tables from each Tourmaline output - rename the files as repseqsV1 and V2 and tableV1 and V2 to represent plates 1-2 and 3-6, respectively. 
qiime feature-table merge-seqs --i-data repseqsv1.qza --i-data repseqsv2.qza --o-merged-data repseqs-16S-merge.qza
qiime feature-table merge --i-tables tablev1.qza --i-tables tablev2.qza --o-merged-table table-16S-merge.qza   

# Infer taxonomy in QIIME 2 with trained 16S V4-V5 classifier from above
qiime feature-classifier classify-sklearn --i-classifier classifier_V4-V5.qza --i-reads repseqs-16S-merge.qza --o-classification taxonomy-16S-merge.qza





















