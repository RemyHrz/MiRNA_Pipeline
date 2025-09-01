# Manual for Virus infected plant miRNA Analysis

### Requirement:

A linux distribution (only tested on gentoo)

mirdeep-P2 v144

Targetfinder https://github.com/carringtonlab/TargetFinder 

Fastp https://github.com/OpenGene/fastp

Bowtie 1.3.1 https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.1

Python (tested with 3.12.8)
    - pandas

R 4.4.1 
    - Deseq2 1.44.0
    - EnhancedVolcano 1.22.0

softwares must be accessible from path

### Filenames:

Avoid using "." in filenames prefer camelCase if you need separation in filename

#### fastq files

raw fastq file must contain "_x" to indicate replicate number ex: Mock_1, Mock_2 ,...

No other underscore are allowed

The control replicates must be named "Mock_x"

Fastq file can be compressed in gz format.

### references:

user must download reference genomes of the host and viruses, and preferably name the sequences with a "clear name" instead of the official ID 
those genomes are to be put in reference/genomes in fasta format

User must download the targets file (aka 3 untranslated region, introns) of the host coding dna

Do not put any files in bowtie/reference

user must download mature mirna sequence (mature.fa) from mirbase

