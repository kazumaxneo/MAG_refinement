# MAG_refinement

Install Anaconda (Mac OS X, Linux).  

```
git clone https://github.com/kazumaxneo/c.git
mamba env create --file MAG_refinement/env.yaml
conda activate binrefinemnet
python Bin_refinement.py -h
usage: Bin_refinement.py [-h] --read1 READ1 --read2 READ2 --pacbio PACBIO [--hifi_mapped_dir HIFI_MAPPED_DIR] [--sr_dir SR_DIR]
                         [--refined_bin_dir REFINED_BIN_DIR]

Run mapping and assembly pipeline.

options:
  -h, --help            show this help message and exit
  --read1 READ1         Path to the first paired-end read file.
  --read2 READ2         Path to the second paired-end read file.
  --pacbio PACBIO       Path to the PacBio reads file.
  --hifi_mapped_dir HIFI_MAPPED_DIR
                        Output directory for HiFi mapped data.
  --sr_dir SR_DIR       Output directory for short-read mapped data.
  --refined_bin_dir REFINED_BIN_DIR
                        Output directory for refined assemblies.
```

## Usage:  
Bin_refinement.py mustbed runned where the bin fasta file is located.  
This script recognises the ".fa" extension.
```
python Bin_refinement.py --reads1 short_R1.fastq.gz --reads2 short_R2.fastq.gz --pacbio HiF_reads.fq.gz --hifi_mapped_dir HIFI_saved_dir --sr_dir short_reads_saved_dir --refined_bin_dir refined_bin_dir
```
    
