##############################################################
### liangweiting@genomics.cn
### 2023/11/28
### metagenomics profile v3.0
### ref: https://github.com/ohmeta/metapi/blob/master/metapi/Snakefile
### ref: https://github.com/liu930724/meta_profile
##############################################################

import os
import sys
import pandas

shell.executable("bash")
configfile: "config.yaml"

def parse_samples(samples_tsv):
    return pandas.read_csv(samples_tsv, sep='\t').set_index("id", drop=False)

def get_sample_id(sample_df, wildcards, col):
    return sample_df.loc[wildcards.sample, [col]].dropna()[0]

_samples = parse_samples("samples.txt") # fixed file name


if config["params"]["metaphlan4"]["do"] and config["params"]["kraken2"]["do"]:
    rule all:
        input:
            expand("{result_dir}/filter_summary.txt", result_dir = config["results"]),
            expand("{result_dir}/metaphlan4.profile.merge.txt", result_dir = config["results"]),
            expand("{result_dir}/bracken.merged.abundance.profile.{level}.tsv", result_dir = config["results"], level=config["params"]["bracken"]["level"])
elif config["params"]["kraken2"]["do"]:
    rule all:
        input:
            expand("{result_dir}/filter_summary.txt", result_dir = config["results"]),
            expand("{result_dir}/bracken.merged.abundance.profile.{level}.tsv", result_dir = config["results"], level=config["params"]["bracken"]["level"])
elif config["params"]["metaphlan4"]["do"]:
    rule all:
        input:
            expand("{result_dir}/filter_summary.txt", result_dir = config["results"]),
            expand("{result_dir}/metaphlan4.profile.merge.txt", result_dir = config["results"])
else:
    rule all:
        input:
            expand("{result_dir}/filter_summary.txt", result_dir = config["results"])




### step1 trimming & remove host reads
### To reduce disk storage usage, merge trimming and remove host together.

rule filter:
    input:
        r1 = lambda wildcards: get_sample_id(_samples, wildcards, "fq1"),
        r2 = lambda wildcards: get_sample_id(_samples, wildcards, "fq2")
    output:
        trim_r1 = temp(os.path.join(config["assay"]["trimming"], "{sample}.trimmed.1.fq.gz")),
        trim_r2 = temp(os.path.join(config["assay"]["trimming"], "{sample}.trimmed.2.fq.gz")),
        html = os.path.join(config["assay"]["trimming"], "{sample}.fastp.html"),
        json = os.path.join(config["assay"]["trimming"], "{sample}.fastp.json"),
        rmhost_r1 = protected(os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.1.fq.gz")),
        rmhost_r2 = protected(os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.2.fq.gz"))
    params:
        min_len = config["params"]["fastp"]["min_len"],
        index = config["params"]["rmhost"]["bowtie2_index"],
        ad_r1 = config["params"]["fastp"]["adapter_r1"],
        ad_r2 = config["params"]["fastp"]["adapter_r2"]
    threads:
        config["params"]["rmhost"]["threads"]
    benchmark:
        os.path.join(config["benchmarks"]["filter"], "{sample}.filter.benchmark.txt")
    log:
        fastp_log = os.path.join(config["logs"]["trimming"], "{sample}.fastp.log"),
        bowtie2_log = os.path.join(config["logs"]["rmhost"], "{sample}.rmhost.log")
    run:
        shell(
        '''
        fastp -i {input.r1} -I {input.r2} -o {output.trim_r1} -O {output.trim_r2} -w {threads} --length_required {params.min_len} --adapter_sequence={params.ad_r1} --adapter_sequence_r2={params.ad_r2} -j {output.json} -h {output.html} 2> {log.fastp_log}
        bowtie2 --end-to-end --very-sensitive -p {threads} -x {params.index} -1 {output.trim_r1} -2 {output.trim_r2} 2> {log.bowtie2_log} | samtools fastq -N -c 5 -f 12 -1 {output.rmhost_r1} -2 {output.rmhost_r2} -
        ''')

### step2 filter_summary
rule seqkit_stat:
    input:
        expand("{rmhost_log_dir}/{{sample}}.rmhost.{reads}.fq.gz", rmhost_log_dir = config["assay"]["rmhost"], reads = ["1","2"])
    output:
        os.path.join(config["logs"]["rmhost"], "{sample}.rmhost.reads.summary")
    shell:
        "seqkit stat {input} > {output}"

rule filter_summary:
    input:
        trim = expand("{trim_res}/{sample}.fastp.json", trim_res = config["assay"]["trimming"], sample = _samples.index),
        rmhost = expand("{rmhost_res}/{sample}.rmhost.reads.summary", rmhost_res = config["logs"]["rmhost"], sample = _samples.index)
    output:
        protected(os.path.join(config["results"], "filter_summary.txt"))
    params:
        trim_summary = temp(os.path.join(config["results"], "trim_summary.txt")),
        rmhost_summary = temp(os.path.join(config["results"], "rmhost_summary.txt"))
    run:
        shell(
        '''
        python rules/filter_summary.py -t {input.trim} > {params.trim_summary}
        python rules/filter_summary.py -r {input.rmhost} > {params.rmhost_summary}
        python rules/merge_summary.py {params.trim_summary} {params.rmhost_summary} {output}
        rm {params.trim_summary} {params.rmhost_summary}
        ''')

### step3 profile

rule metaphlan4:
    input:
        rmhost_r1 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.1.fq.gz"),
        rmhost_r2 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.2.fq.gz")
    output:
        mpa = protected(os.path.join(config["assay"]["profile"]["metaphlan4"], "{sample}.mp4.profile")),
        bw2 = protected(os.path.join(config["assay"]["profile"]["metaphlan4"], "{sample}.mp4.bw2.bz2")),
        sam = protected(os.path.join(config["assay"]["profile"]["metaphlan4"], "{sample}.sam.bz2"))
    threads:
        config["params"]["metaphlan4"]["threads"]
    params:
        bowtie2db = config["params"]["metaphlan4"]["bowtie2db"],
        index = config["params"]["metaphlan4"]["index"],
    log:
        os.path.join(config["logs"]["profile"], "{sample}.metaphlan4.log")
    benchmark:
        os.path.join(config["benchmarks"]["profile"], "{sample}.mp4.benchmark.txt")
    shell:
        '''
        export PATH=/ldfssz1/ST_META/share/User/chenjh356/anaconda3/envs/mpa/bin:$PATH        
        metaphlan {input.rmhost_r1},{input.rmhost_r2} --bowtie2out {output.bw2} --nproc {threads} --input_type fastq -s {output.sam} -t rel_ab_w_read_stats --bowtie2db {params.bowtie2db} --index {params.index} -o {output.mpa} 2> {log}      
        '''


rule merge_profile_mp4:
    input:
        mpa = expand("{profile_dir}/{sample}.mp4.profile", profile_dir = config["assay"]["profile"]["metaphlan4"], sample = _samples.index),
    output:
        mpa_merge = protected(os.path.join(config['results'], "metaphlan4.profile.merge.txt")),
    shell:
        '''
        merge_metaphlan_tables.py {input.mpa} > {output.mpa_merge}
        '''


rule kraken2:
    input:
        rmhost_r1 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.1.fq.gz"),
        rmhost_r2 = os.path.join(config["assay"]["rmhost"], "{sample}.rmhost.2.fq.gz")
    output:
        result = protected(os.path.join(config["assay"]["profile"]["kraken2"], "{sample}.result")),
        report = protected(os.path.join(config["assay"]["profile"]["kraken2"], "{sample}.report.txt"))
    threads:
        config["params"]["kraken2"]["threads"]
    params:
        db = config["params"]["kraken2"]["db"]
    log:
        os.path.join(config["logs"]["profile"], "{sample}.kraken2.log")
    benchmark:
        os.path.join(config["benchmarks"]["profile"], "{sample}.mp3.kraken2.txt")
    shell:
        '''
        kraken2 --use-names --db {params.db} --threads {threads} --report {output.report} --output {output.result} --gzip-compressed {input.rmhost_r1} {input.rmhost_r2} 2> {log}
        '''

rule bracken:
    input:
      report = os.path.join(config["assay"]["profile"]["kraken2"], "{sample}.report.txt")
    output:
      profile = protected(os.path.join(config["assay"]["profile"]["kraken2"],"{sample}.bracken.{level}.profile")),
      report = protected(os.path.join(config["assay"]["profile"]["kraken2"],"{sample}.bracken.{level}.report"))
    params:
      db = config["params"]["kraken2"]["db"],
      reads_len = config["params"]["bracken"]["reads_len"],
      level = "{level}"
    threads:
      config["params"]["kraken2"]["threads"]
    shell:
        '''
        bracken -d {params.db} -i {input.report} -o {output.profile} -w {output.report} -r {params.reads_len} -l {params.level} -t {threads}
        '''

rule merge_profile_bracken:
    input:
       profile = expand(os.path.join(config["assay"]["profile"]["kraken2"], "{sample}.bracken.{{level}}.profile"), sample = _samples.index.unique())
    output:
       merge_profile = protected(os.path.join(config["results"], "bracken.merged.abundance.profile.{level}.tsv"))
    run:
       shell(
         '''
         rules/combine_bracken_outputs.py --files {input.profile} --names %s --output {output.merge_profile}
         ''' % ",".join(_samples.index.unique())
         )
