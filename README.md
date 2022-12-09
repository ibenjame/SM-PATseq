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
* Create the conda environment for PATseq analysis using the provided requirements file:
```
conda create --name PATseq --file requirements.txt
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

* Check the job plan with a snakemake dry run
```
snakemake -n
```

* After seeing a job plan without errors, submit the snakemake job.  If on a job scheduler such as slurm on a cluster, use a wrapper for submitting the and distriubting the jobs appropriately.

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
