# MAG_refinement
  
MAG sequences derived from short reads consist of a set of fragmented contigs that are cutted in hundreds of position. Gene prediction is severely affected by such fragmented draft genomes. This script aims to improve MAG contiguity by reassembling bin using a pair of properly mapped HiFi long reads and short reads. <br><br>  
Three isolated genome assembler was used to reassembling bins.
1. Unicycler (hybrid assembly)
2. SPAdes (hybrid assembly)
3. Flye (HiFi long reads assembly)
<br><br>
## Requirements  
- SAMTools
- Bowtie2 
- samclip
- minimap2
- Flye v2.9
- SPAdes v3.15
- Unicycler v5  
<br><br>
## Instalation  

```
git clone https://github.com/kazumaxneo/c.git
mamba env create --file MAG_refinement/env.yaml
conda activate binrefinemnet
python Bin_refinement.py -h
```
<br>
  
## Usage  
Bin_refinement.py mustbed runned where the bin fasta file is located. This script recognises the ".fa" extension.
```
python Bin_refinement.py --reads1 short_R1.fastq.gz --reads2 short_R2.fastq.gz --pacbio HiF_reads.fq.gz --hifi_mapped_dir HIFI_saved_dir --sr_dir short_reads_saved_dir --refined_bin_dir refined_bin_dir
```

<br>
## Options

- **`-h, --help`**  
  Show this help message and exit.

- **`--read1 READ1`**  
  Path to the first paired-end read file.

- **`--read2 READ2`**  
  Path to the second paired-end read file.

- **`--pacbio PACBIO`**  
  Path to the PacBio reads file.

- **`--hifi_mapped_dir HIFI_MAPPED_DIR`**  
  Output directory for HiFi mapped data.

- **`--sr_dir SR_DIR`**  
  Output directory for short-read mapped data.

- **`--refined_bin_dir REFINED_BIN_DIR`**  
  Output directory for refined assemblies.  

<br><br>
## How to cite  
These tool should be cited.<br> 
- Bankevich A, Nurk S, Antipov D, Gurevich AA, Dvorkin M, Kulikov AS, Lesin VM, Nikolenko SI, Pham S, Prjibelski AD, Pyshkin AV, Sirotkin AV, Vyahhi N, Tesler G, Alekseyev MA, Pevzner PA. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. J Comput Biol. 2012 May;19(5):455-77. doi: 10.1089/cmb.2012.0021. Epub 2012 Apr 16. PMID: 22506599; PMCID: PMC3342519.  

- Freire B, Ladra S, Parama JR. Memory-Efficient Assembly Using Flye. IEEE/ACM Trans Comput Biol Bioinform. 2022 Nov-Dec;19(6):3564-3577. doi: 10.1109/TCBB.2021.3108843. Epub 2022 Dec 8. PMID: 34469305.  

- Li H. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics. 2018 Sep 15;34(18):3094-3100. doi: 10.1093/bioinformatics/bty191. PMID: 29750242; PMCID: PMC6137996.  

- Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R; 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009 Aug 15;25(16):2078-9. doi: 10.1093/bioinformatics/btp352. Epub 2009 Jun 8. PMID: 19505943; PMCID: PMC2723002.  

- Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. Nat Methods. 2012 Mar 4;9(4):357-9. doi: 10.1038/nmeth.1923. PMID: 22388286; PMCID: PMC3322381.  

- Wick RR, Judd LM, Gorrie CL, Holt KE. Unicycler: Resolving bacterial genome assemblies from short and long sequencing reads. PLoS Comput Biol. 2017 Jun 8;13(6):e1005595. doi: 10.1371/journal.pcbi.1005595. PMID: 28594827; PMCID: PMC5481147.
<br><br>

## Licence ##

GPL v3.
