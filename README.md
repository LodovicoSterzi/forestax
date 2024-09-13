# Forestax: define bacterial taxonomy using the RecA protein marker

## What does Forestax do?
Forestax is a fast, easy to use tool to assign bacterial taxonomy starting from raw reads, FASTA nucleotide file (such as genome assemblies) or FASTA protein sequences files. Forestax takes the input, extracts RecA protein sequences and assigns the bacterial taxonomy to each sequence using a Random forest algorithm. If raw reads are given, the tool also provides bacterial abundances using the reads coverage. More information about the tool is available at:

  RecA is a reliable marker for bacterial taxonomy, even in the Candidate Phyla Radiation
Lodovico Sterzi, Simona Panelli, Clara Bonaiti, Stella Papaleo, Giorgia Bettoni, Enza D’Auria, Gianvincenzo Zuccotti, Francesco Comandatore
bioRxiv 2024.06.21.600076; doi: https://doi.org/10.1101/2024.06.21.600076

## How to install Forestax
Forestax is a Linux-based software. The tool will be soon available for installation by conda, but until then it can be used as follows.
### Download:
Download the tool from github and unzip:
```
git clone https://github.com/LodovicoSterzi/forestax
unzip forestax.zip
cd Forestax/
```
### Dependencies
The tool performs a pipeline which requires:
- python (>=3.6.8)
- R (>=3.5.0)
- Prodigal (=2.6.1)
- SPAdes (=3.15.5)
- Diamond (=2.1.8)
- Bowtie2 (=2.4.2)
- Samtools (=1.19.1)

And the following R libraries:
- ranger
- data.table
- stringr
- ape
- phylotools
- peakRAM

The list of dependencies might look scary but only **python, R and the R libraries must be installed manually**. The secondary tools can be installed automatically in the bin/ folder running:
```
python3 install_dependencies.py
```
If you wish to install the tools elsewhere, you already have them or if the automatic installation fails, the correct paths must be manually given to the script ***forestax.py*** (line 26-30) as shown here:
```
        --> forestax.py line 26-30
# Paths to tools: CHANGE HERE IF NEEDED
diamond_path = "my/path/to/the/tool/diamond"
prodigal_path = "my/path/to/the/tool/prodigal"
bowtie_path = "my/path/to/the/tool/bowtie2-2.4.2-sra-linux-x86_64")
samtools_path = "my/path/to/the/tool/samtools"
spades_path = "my/path/to/the/tool/"metaspades.py"
```
### Check:
Once the tool has been downloaded and dependencies are ready, the tool can be run as follows:
```
python path/to/forestax/forestax.py -h
```

## How to run Forestax
### Parameters:
```
  -h, --help            show this help message and exit
  -o, --output          Output folder
  -itype, --input_type  Input type (reads, nucleotide or proteins)
  -r1, --reads_forward  If the input type are reads, path of FASTQ forward reads file
  -r2, --reads_reverse  If the input type are reads, path of FASTQ reverse reads file
  -f, --fasta           If the input type are nucleotide or proteins, path of the FASTA file
  -p, --cpu             Number of threads (default: 1)
  -tax, --tax_to_assign Comma-delimited list of taxonomic levels to assign (default: phylum,class,order,family,genus,species)
  -b, --bowtie_mode     Bowtie2 mode used to align reads to RecA bowtie database (default: very-sensitive)
  -k, --keep_track      Keep track of RAM usage/time for each step. Beware: it will slow down a bit the analysis!
  -q, --quantify        Use reads coverage to infer relative abundances of the taxa
  -r, --reads_paths     If the input is a nucleotide FASTA file (eg. a genome assembly) but wish to quantify the bacterial taxa, 
                        the tool must be given a text file with THREE TAB-DELIMITED columns and NO HEADER: 
                        name-of-input-file | path-of-forward-reads-file | path-of-reverse-reads-file
```
### Basic usage
The tool must be given the type of input (reads,nucleotide or protein), the input file/s and the name of the output folder. 

#### From paired FASTQ reads (default)
```
python ./forestax.py -r1 example_inputs/reads/*1.fastq.gz -r2 example_inputs/reads/*2.fastq.gz -o example_reads
```
#### From FASTA nucleotide
```
python ./forestax.py -itype nucleotide -f example_inputs/nucleotide/*.fasta -o example_nucleotide
```
#### From FASTA proteins
```
python ./forestax.py -itype proteins -f example_inputs/proteins/*.faa -o example_proteins
```

### Other usage options
It is also possible to customise some options. For example, the user can ask the tool to quantify the abundance of each taxon
```
python ./forestax.py -r1 example_inputs/reads/*1.fastq.gz r2 example_inputs/reads/*2.fastq.gz -o example_reads -q
```
or customise the taxonomic depth to reach in the assignment (by default it arrives to the species).
```
python ./forestax.py -r1 example_inputs/reads/*1.fastq.gz r2 example_inputs/reads/*2.fastq.gz -o example_reads -tax phylum,class,order,family
```
If the user wants to estimate taxon abundance but would like to start from a FASTA file, it is possible to do so by adding a tab-delimited plain text file 
with the name of the FASTA input and the paths to the corresponding reads files. 
```
python ./forestax.py -itype nucleotide -f example_inputs/nucleotide/*.fasta -o example_nucleotide -r reads_path_file.tab
```
With reads_path_file.tab written like this:
```
example1.fasta  example_inputs/reads/example1_1.fastq.gz  example_inputs/reads/example1_2.fastq.gz
example2.fasta  example_inputs/reads/example2_1.fastq.gz  example_inputs/reads/example2_2.fastq.gz
```
## The output
Forestax will start by creating a folder for the run. In the main folder, there will be: i) an **Output** folder containing the three output files; ii) an **Analysis** folder containing temporary files for each step of the pipeline;  iii) a **log file** containing information about each step of the pipeline.
The output folder contains three files:
- *RecA_Taxonomy_OUTPUT.tab*: this is the main output file and contains the taxonomic assignment and probability at each taxonomic level for each sequence. If (-q) has been flagged, it also gives the reads mean coverage for each sequence. 
- *RecA_Taxonomy_OUTPUT.extended.tab*: this file contains extra information such as the amino acid sequence, the dayhoff-compressed sequence and the DIAMOND hit coverage values. 
- *RecA_sequences.fasta*: this file the RecA sequences. 











