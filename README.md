# Poly-A tail Sequencing (PATseq) Snakemake example workflow for working with PacBio Sequel subread data to analysis ready report files

## Description

A simple Snakemake file using a conda environment built from requirements.txt and helper scripts.  This takes subread data resulting from a Poly-A tail reading experiment from the PacBio Sequel platform, assigns gene identities to the transcripts based on a provided transcriptome fasta, evaluates reads for poly-A tail lengths, then tabulates a final analysis ready table of results.

## Getting Started

### Dependencies

* Anaconda or miniconda for use of conda environments.

### Installing

* Download the repository
```
git clone https://github.com/ibenjame/PATseq
cd PATseq
```
* Install Mamba to your base conda environment if not already present
```
conda install mamba -n base -c conda-forge
```
* Create the conda environment for PATseq analysis using the provided env.yaml file using the mamba solver
```
mamba env create --name PATseq --file env.yaml
```
* Load the environment
```
conda activate PATseq
```

### Executing program

* For each PacBio Run in the analysis containing multiplexed PATseq reads, create a symbolic link of the run folder to a subfolder of the form PacBio{X}.  Where X is a number counting up from 1 to the total number of runs for the samples.  Similar to below (where r54233_20220815_152625 represents your first run folder).
```
ln -s /path/to/raw_data/r54233_20220815_152625 PacBio1
```

* Each folder should be of a form similar to the below, with subfolders within for each sequencing cell (A-H here in this example) in the run and a subreads bam contained within each cell folder.
```
PacBio1
├── 1_A01
├── 2_B01
├── 3_C01
├── 4_D01
├── 5_E01
├── 6_F01
├── 7_G01
└── 8_H01
```

* Edit a samples.txt file to describe the multiplexed samples in the experiment.  This is in the form of header-containing tab-delimited text similar to below.  Columns should be Index (referring to the barcode index used for the sample; found in scripts/barcodes.master.fasta) and ID (a name for the sample and associated filenames; Do not use any special characters or whitespace)
```
Index   ID
89      WT_mRNA_1
90      WT_mRNA_2
91      WT_mRNA_3
92      MUT_mRNA_1
93      MUT_mRNA_2
94      MUT_mRNA_3
```

* Locate or obtain your transcriptome fasta (ideally these should contain UTR sequence and not be just coding).  Edit config.yaml to reflect the path and filename of your transcriptome fasta file.

* Adjust any other config.yaml settings as necessary.  Ideally you should be using a transcriptome from either GENCODE or Ensembl.  Downstream reporting scripts may not parse other gene name inputs properly and may beed to be adjusted to your specific source format.

* Check the job plan with a snakemake dry run
```
snakemake -n
```

* After seeing a job plan without errors, submit the snakemake job.  If on a job scheduler such as slurm on a cluster, use a wrapper for submitting the and distriubting the jobs appropriately.

### Outputs

* Several folders will be created and populated during the process
  * Folder "merged_bams" will contain individual sample CCS bam files combining all flowcells in the experiment into one.  These are named based on the samples.txt file designations provided of the form {name}.merged.ccs.bam
  * Folder "reads" will contain gzipped fastq and fasta format read files for each sample of the form {name}.fastq.gz and {name}.fa.gz
  * Each sample should have an output folder by name containing several files:
    * A sorted and indexed bam file of minimap2 alignments of the form {name}.mm2.bam and {name}.mm2.bam.bai
    * A log file len.{name}.log summarizing overall tail finding stats
    * A report summary containing overall stats related to both tail and gene assignment rates {name}.report.summary
    * A tab-delimited report text file for the sample containing gene-level poly-A tail statistics {name}.report.txt
    * A length assignment zip file containing only tail identification stats per read {name}.lengths.txt.gz
    * A full tab-delimited zipped report file with per-read tail and gene assignment {name}.full.txt.gz

* The sample table {name}.full.txt.gz is a headerless tab-delimited text file with the following fields:
  1. ReadID – string from the sequencing/ccs generation
  2. Length – PAT length determined from the read
  3. Direction – If CCS read is oriented in forward or reverse compliment direction (random based on initial loop priming)
  4. C – Number of non-A C bases found in the poly-A tail
  5. G – Number of non-A G bases found in the poly-A tail
  6. T – Number of non-A T bases found in the poly-A tail
  7. TerminalUs – Number of terminal Uridyls found AFTER poly-A tail (and before adapter)
  8. Gene – Gene symbol designation for the read found by minimap2 alignment
  9. AlignmentScore – Score for the read alignment to the denoted gene
  10. GencodeString – (Unused if GENCODE is not selected in the config; for Gencode transcriptome alignments the full assignment string rather than just symbol)

* Note that non-A bases are defined as occuring within the tail as long as they are bracketed by a minimum of 2 A bases on either side.  We make no assumption that these are true biological in origin vs introduced at some point during reverse transcription and amplification, similar to possible terminal Uridiylation reported stats.  These are reported for convenience and potential QC.  Be cautious when making use of these statistics.

## Authors

James Iben (ibenjame@nih.gov)


## Version History


* 0.1
    * Initial Release

## License

MIT License

Copyright (c) 2022 James R Iben

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Acknowledgments

We thank HPC staff (hpc.nih.gov) for computational resources related to this work and Ryan Dale for expert advice regarding the workflow
