## For PacBio PAT data
##Dependencies: requirements.txt; install viridis into R

# Expects one or more "PacBio" prefixed folders of runs (cells)
# a barcodes.fa containing the expected barcodes (by default in the scripts subfolder)
# Needs samples.txt; 2 columns tab delimited with index# (barcodes.fa) and sampleNames
# config.yaml should be updated with apprpriate seetings and path to your transcripts fasta for assignment
# Existing perl scripts used: "parse_fastas_for_lengths_with_U.pl", "len_assign_minimap.pl", "minimap2_parse.pl"
# For PAT length discovery (stats), Full assigned table generation (for R), and Gene report file
# By default these scripts expect to work with GENCODE or Ensembl naming conventions for the transcripts fasta

import os
import csv
configfile: "config.yaml"

#Sample ID file
sampleids = {}
with open('samples.txt','r') as f:
    next(f) # skip headings
    reader=csv.reader(f,delimiter='\t')
    for index,id in reader:
        sampleids.update({id : index})

#Identify runs
subreads = []
folders = []
for root, dirs, files in os.walk(".", followlinks=True):
    for file in files:
        if file.endswith("subreads.bam"):
            path = os.path.join(root, file)
            subreads.append(path[2:])
            folders.append(root[2:])

#Show runs
#print(subreads)
#Show samples
#print(sampleids)
#Show folders
#print(folders)

#Simplified variable names from config
transcripts=config["transcripts"]


rule targets:
    input:
        expand('{folder}/ccs.bam', folder=folders),
        expand('{folder}/split.lima.report', folder=folders),
        expand('merged_bams/{sample}.merged.ccs.bam', sample=sampleids.keys()),
        expand('reads/{sample}.fa.gz', sample=sampleids.keys()),
        expand('{sample}/{sample}.report.txt', sample=sampleids.keys()),
        expand('{sample}/{sample}.full.txt.gz', sample=sampleids.keys())

#Make CCSes of each cell
rule ccs:
    output:
        ccs="PacBio{r}/{cell}/ccs.bam"
    threads: 20
    shell:
        "ccs -j {threads} --report-file PacBio{wildcards.r}/{wildcards.cell}/ccs_report.txt PacBio{wildcards.r}/{wildcards.cell}/*.subreads.bam {output.ccs}"

#Split samples by lima
rule lima:
    input:
        ccs="PacBio{r}/{cell}/ccs.bam",
        bcs="scripts/barcodes.master.fasta"
    output:
        report="PacBio{r}/{cell}/split.lima.report"
    threads: 8
    params: detail="scripts/report_detail.R"
    shell:
        "lima --peek-guess --split-bam-named -s -j {threads} {input.ccs} {input.bcs} PacBio{wildcards.r}/{wildcards.cell}/split.bam;"
        "cd PacBio{wildcards.r}/{wildcards.cell}; Rscript --vanilla {params.detail} split.lima.report pdf"

#Samtools merge the barcode split files by whitelist
rule combine:
    input:
        expand('{folder}/split.lima.report', folder=folders)
    output:
        "merged_bams/{id}.merged.ccs.bam"
    threads: 8
    params: value=lambda wcs: sampleids[wcs.id]
    shell:
        "mkdir -p merged_bams;"
        "samtools merge -f -@ {threads} {output} PacBio*/*/split.lbc{params.value}--lbc{params.value}.bam"

#Extract fasta reads
rule extract:
    input:
        "merged_bams/{id}.merged.ccs.bam"
    output:
        fq="reads/{id}.fastq.gz",
        fa="reads/{id}.fa.gz"
    shell:
        "mkdir -p reads;"
        "samtools fastq {input} | gzip - > {output.fq};"
        "samtools fasta {input} | gzip - > {output.fa};"

#Find lengths
rule lengths:
    input:
        fa="reads/{id}.fa.gz",
        pl="scripts/parse_fastas_for_lengths_withU.pl"
    output:
        "{id}/{id}.lengths.txt.gz"
    log: "{id}/len.{id}.log"
    params:
        fadapter=config["forward_adapter"],
        radapter=config["reverse_adapter"]
    shell:
        "mkdir -p {wildcards.id};"
        "perl {input.pl} {input.fa} {params.fadapter} {params.radapter} 2>{log} | gzip - > {output};"

#Minimap2 of sequences
rule minimap2:
    input:
        read="reads/{id}.fa.gz"
    output:
        "{id}/{id}.mm2.bam"
    shell:
        "mkdir -p {wildcards.id};"
        "minimap2 -a {transcripts} {input.read} | samtools sort -O BAM - > {output};"
        "samtools index {output}"

#Extend Length report
rule extend:
    input:
        len="{id}/{id}.lengths.txt.gz",
        minimap="{id}/{id}.mm2.bam",
        pl="scripts/len_assign_minimap.pl"
    output:
        "{id}/{id}.full.txt.gz"
    params: gencode=config["GENCODE_or_ENSEMBL"]
    shell:
        "perl {input.pl} {input.minimap} {input.len} {params.gencode} | gzip - > {output}"

#Build sample gene report file
rule report:
    input:
        blat="{id}/{id}.mm2.bam",
        len="{id}/{id}.lengths.txt.gz",
        pl="scripts/minimap2_parse.pl"
    output:
        "{id}/{id}.report.txt"
    log: "{id}/{id}.report.summary"
    params: gencode=config["GENCODE_or_ENSEMBL"]
    shell:
        "perl {input.pl} {input.blat} {input.len} {params.gencode} 2>{log} >{output}"

