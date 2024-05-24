# metaprof v2.2
A quick profiling about metagenomes data
--update metaphlan4

# environment
  - biopython>=1.76
  - bowtie2>=2.3.5.1
  - fastp>=0.20.1
  - metaphlan>=4.1.0
  - numpy>=1.18.4
  - pandas>=1.0.3
  - pigz>=2.3.4
  - samtools>=1.9
  - seqkit>=0.12.1
  - snakemake>=5.14.0
  - kraken2>=2.1.1
  - bracken>=2.5

# install

```
git clone https://github.com/weiting-liang/metaprof.git

cd metaprof
conda env create -n metaprof -f ./rules/env.yaml
conda activate metaprof
```


#database prepare  
human reference
https://github.com/marbl/CHM13
```
mkdir /path/database
cd database && mkdir humanhost && cd humanhost 
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
gzip chm13v2.0.fa.gz -d
mkdir bowtie2_index && cd bowtie2_index
bowtie2-build ../chm13v2.0.fa chm13v2
```

metaphlan3
```
cd /path/database
mkdir metaphlan && cd metaphlan
metaphlan --install --index mpa_vOct22_CHOCOPhlAnSGB_202212 --bowtie2db metaphlan_database
```

kraken2
```
cd /path/database
mkdir kraken2 && cd kraken_pub
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20210517.tar.gz
tar -xzvf k2_standard_20210517.tar.gz -C ./k2_standard_20210517
```

# run
- change the `samples.txt` to adapt to your data: separator should be "\t": id^Ifq1^Ifq2$
- custom the `config.yaml` database's path and parameters

```
#dry run
snakemake --snakefile rules/profile.smk -n
#test
snakemake --snakefile rules/profile.smk --core 16 2> smk.log &
#cluster: custom the cluster.yaml
nohup sh snakemake.sh &
```

# output

2.results/  
  filter_summary.txt  
  metaphlan4.profile.merge.txt 
  bracken.merged.abundance.profile.*.tsv  

#assay:  
1.assay  
  01.trimming/  
  02.rmhost/  
  03.profile/  
  benchmarks/   #check the cpu's time and max_vms to optimize the cluster's parameters  
  cluster_logs/   
  logs/         #find programs' errors  


# references
https://github.com/ohmeta/metapi/blob/master/metapi/Snakefile

https://github.com/liu930724/meta_profile
