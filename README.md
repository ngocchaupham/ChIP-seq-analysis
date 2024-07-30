
# ChIP-seq Analysis Pipeline

This project is a comprehensive pipeline for analyzing ChIP-seq sequencing data. The workflow includes downloading sequencing data, performing quality control, trimming, aligning reads to a reference genome, and calling peaks. The workflow is managed using a Makefile and requires certain software tools and modules to be available in your environment.

## Project Structure

- **WORKDIR**: `/shared/projects/2324_tp_m2_bsg/Chau/`  
  The working directory where the analysis will be conducted.

- **INDEX**: `/shared/bank/homo_sapiens/hg38/bowtie2/hg38`  
  The reference genome index used for alignment with Bowtie2.

- **GENOME_SIZE**: `250e6`  
  The estimated size of the reference genome (adapt this value based on the genome you are working with).

## Input Files

- **sample_list.txt**: A file containing sample identifiers and related metadata. Each line should contain two fields separated by a semicolon `;`, representing paired samples (e.g., experimental and control).

## Workflow Rules

### 1. **Download**

- **Output**: `fastq/{smp}.fq`
- **Description**: Downloads and decompresses sequencing data files from Zenodo.
- **Command**:
  ```shell
  wget https://zenodo.org/records/7646128/files/{wildcards.smp}.fq.gz -O fastq/{wildcards.smp}.fq.gz
  gunzip fastq/{wildcards.smp}.fq.gz
  ```

### 2. **FastQC**

- **Output**: `fastqc/{smp}_fastqc.html`
- **Description**: Performs quality control on the raw FASTQ files using FastQC.
- **Command**:
  ```shell
  fastqc -o fastqc {input}
  ```

### 3. **Sickle**

- **Output**: `sickle/{smp}.fq.gz`
- **Description**: Trims and filters the FASTQ files to remove low-quality bases and reads using Sickle.
- **Command**:
  ```shell
  sickle se -f {input} -o {output} -t sanger -g
  ```

### 4. **Bowtie2 Alignment**

- **Output**: `bowtie/{smp}.bam`
- **Description**: Aligns the trimmed reads to the reference genome using Bowtie2 and converts the output to BAM format.
- **Command**:
  ```shell
  bowtie2 -p 1 -U {input} -x {params.index} | samtools view -b | samtools sort > {output}
  ```

### 5. **BAM Indexing**

- **Output**: `bowtie/{smp}.bam.bai`
- **Description**: Indexes the BAM file using Samtools for downstream analysis.
- **Command**:
  ```shell
  samtools index {input}
  ```

### 6. **BAM Coverage**

- **Output**: `bam_coverage/{smp}.bw`
- **Description**: Generates coverage files in BigWig format using DeepTools.
- **Command**:
  ```shell
  bamCoverage -b {input.bam} -bs 25 -o {output} --region chr1
  ```

### 7. **MACS2 Peak Calling**

- **Output**: `macs2/{smp}_peaks.narrowPeak`
- **Description**: Calls peaks on aligned reads using MACS2 for detecting enriched regions.
- **Command**:
  ```shell
  macs2 callpeak -t {input.chip_bam} -c {input.input_bam} -n {wildcards.smp} -f BAM -g {params.genome_size} --outdir macs2
  ```

## Setup Instructions

1. **Modules Required**:
   - FastQC 0.11.9
   - Sickle 1.33
   - Bowtie2 2.5.1
   - Samtools 1.15.1
   - DeepTools 3.5.0
   - MACS2 2.2.7.1

2. **Run the Pipeline**:
   - Ensure all necessary modules are loaded or installed in your environment.
   - Execute the Makefile using the following command:
     ```bash
     make
     ```
   - The pipeline will execute each rule in sequence and generate the corresponding output files in their respective directories.

3. **File Organization**:
   - `fastq/`: Contains downloaded and decompressed FASTQ files.
   - `fastqc/`: Contains FastQC reports.
   - `sickle/`: Contains trimmed FASTQ files.
   - `bowtie/`: Contains aligned BAM files and indices.
   - `bam_coverage/`: Contains BigWig coverage files.
   - `macs2/`: Contains peak-calling results from MACS2.

## Notes

- Update the `GENOME_SIZE` and `INDEX` paths according to the genome you are analyzing.
- Ensure that the `sample_list.txt` file is properly formatted for accurate parsing and sample identification.

This README provides an overview of the workflow and guides you on setting up and executing the analysis pipeline efficiently.
