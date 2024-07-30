WORKDIR = "/shared/projects/2324_tp_m2_bsg/Chau/"
INDEX  = "/shared/bank/homo_sapiens/hg38/bowtie2/hg38"
SAMPLE_DICT = dict()
#Warning adapt depending on genome
GENOME_SIZE = 250e6

for line in open("sample_list.txt"):
	line = line.rstrip("\n")
	if line != "":
		cols = line.split(";")
		SAMPLE_DICT[cols[0]] = cols[1]

SAMPLE = list(SAMPLE_DICT.keys()) + list(SAMPLE_DICT.values())
SAMPLE = list(set(SAMPLE)) # enlever la redondance 

workdir: WORKDIR

rule final:
	input: expand("fastqc/{smp}_fastqc.html", smp=SAMPLE), \
			expand("bam_coverage/{smp}.bw", smp=SAMPLE), \
			expand("macs2/{smp}_peaks.narrowPeak", smp=list(SAMPLE_DICT.keys()))
	params: cpu="1"

rule download: 
	output: "fastq/{smp}.fq"
	params: cpu="1"
	shell: """
	mkdir -p fastq
	wget https://zenodo.org/records/7646128/files/{wildcards.smp}.fq.gz -O fastq/{wildcards.smp}.fq.gz
	gunzip fastq/{wildcards.smp}.fq.gz
	"""

rule fastqc:
	input: "fastq/{smp}.fq"
	output: "fastqc/{smp}_fastqc.html"
	params: cpu="1"
	shell: """
	mkdir -p fastqc
	module load fastqc/0.11.9
	fastqc -o fastqc {input}
	"""

rule sickle:
	input: "fastq/{smp}.fq"
	output: "sickle/{smp}.fq.gz"
	params: cpu="1"
	shell: """
	mkdir -p sickle
	module load sickle-trim/1.33
	sickle se -f {input} -o {output} -t sanger -g
	"""

rule bowtie:
	input: "sickle/{smp}.fq.gz"
	output: "bowtie/{smp}.bam"
	params: index=INDEX, cpu="1"
	shell: """
	mkdir -p bowtie
	module load bowtie2/2.5.1
	module load samtools/1.15.1
	bowtie2 -p 1 -U {input} -x {params.index} | samtools view -b |samtools sort > {output}
	"""

rule bam_index:
	input: "bowtie/{smp}.bam"
	output: "bowtie/{smp}.bam.bai"
	params: cpu="1"
	shell: """
	module load samtools/1.15.1
	samtools index {input}
	"""

rule bam_coverage: 
	input: bam = "bowtie/{smp}.bam",\
			bai = "bowtie/{smp}.bam.bai"
	output: "bam_coverage/{smp}.bw"
	params: cpu="1"
	shell: """
	mkdir -p bam_coverage
	module load deeptools/3.5.0
	bamCoverage -b {input.bam} -bs 25 -o {output} --region chr1
	"""

def get_input_bam(wildcard): 
	return "bowtie/" + SAMPLE_DICT[wildcard.smp] + ".bam"

def get_input_bai(wildcard): 
	return "bowtie/" + SAMPLE_DICT[wildcard.smp] + ".bam.bai"

rule macs2: 
	input: chip_bam = "bowtie/{smp}.bam", \
			chip_bai = "bowtie/{smp}.bam.bai", \
			input_bam = get_input_bam, \
			input_bai = get_input_bai
	output: "macs2/{smp}_peaks.narrowPeak"
	params : genome_size = GENOME_SIZE, cpu="1"
	shell: """
	mkdir -p macs2
	module load macs2/2.2.7.1
	macs2 callpeak -t {input.chip_bam} -c {input.input_bam} -n {wildcards.smp} -f BAM -g {params.genome_size}\
		--outdir macs2
	"""
