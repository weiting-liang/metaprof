# metaprof
A quick profiling about metagenomes data


# environment
  - biopython>=1.76
  - bowtie2>=2.3.5.1
  - fastp>=0.20.1
  - metaphlan>=3.0
  - numpy>=1.18.4
  - pandas>=1.0.3
  - pigz>=2.3.4
  - samtools>=1.9
  - seqkit>=0.12.1
  - snakemake>=5.14.0
  - kraken2>=2.1.1
  - bracken>=2.5

# how to run
git clone https://github.com/weiting-liang/metaprof.git

conda env create -n metaprof -f ./rules/env.yaml
conda activate metaprof

# database prepare
mkdir /path/database
cd database && mkdir humanhost && cd humanhost
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
mkdir bowtie2_index && cd bowtie2_index
gzip ../chm13v2.0.fa.gz -d
bowtie2-build chm13v2.0.fa chm13v2

cd /path/database
mkdir metaphlan
metaphlan --install --index mpa_v30_CHOCOPhlAn_201901 --bowtie2db metaphlan_database

cd /path/database
mkdir kraken2
















