#!/usr/bin/env python
import os
import argparse
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Run mapping and assembly pipeline.")
    parser.add_argument("--read1", required=True, help="Path to the first paired-end read file.")
    parser.add_argument("--read2", required=True, help="Path to the second paired-end read file.")
    parser.add_argument("--pacbio", required=True, help="Path to the PacBio reads file.")
    parser.add_argument("--hifi_mapped_dir", default="Hifi_mapped", help="Output directory for HiFi mapped data.")
    parser.add_argument("--sr_dir", default="sr", help="Output directory for short-read mapped data.")
    parser.add_argument("--refined_bin_dir", default="refined_bin", help="Output directory for refined assemblies.")
    args = parser.parse_args()

    # Ensure directories exist
    os.makedirs(args.hifi_mapped_dir, exist_ok=True)
    os.makedirs(args.sr_dir, exist_ok=True)
    os.makedirs(args.refined_bin_dir, exist_ok=True)

    # Find bin*.fa files in the current directory
    bin_files = [f for f in os.listdir(".") if f.startswith("bin") and f.endswith(".fa")]

    for ref in bin_files:
        base_name = os.path.splitext(ref)[0]
        hifi_output = f"{args.hifi_mapped_dir}/{base_name}_map-Hifi_mapped.fastq.gz"
        log_file = f"log_{base_name}.txt"

        # Index the reference
        subprocess.run(["samtools", "faidx", ref])

        # Mapping and filtering
        minimap_cmd = [
            "minimap2", "-ax", "map-hifi", "-t", "20", ref, args.pacbio
        ]
        samtools_view_cmd = [
            "samtools", "view", "-@4", "-F4", "-h", "-", "-q20"
        ]
        samclip_cmd = [
            "samclip", "--max", "100", "--ref", ref, "-"
        ]
        samtools_sort_cmd = ["samtools", "sort", "-@4"]
        samtools_fastq_cmd = ["samtools", "fastq"]
        pigz_cmd = ["pigz", "-p12"]

        with open(hifi_output, "wb") as outfile:
            p1 = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(samtools_view_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)
            p3 = subprocess.Popen(samclip_cmd, stdin=p2.stdout, stdout=subprocess.PIPE)
            p4 = subprocess.Popen(samtools_sort_cmd, stdin=p3.stdout, stdout=subprocess.PIPE)
            p5 = subprocess.Popen(samtools_fastq_cmd, stdin=p4.stdout, stdout=subprocess.PIPE)
            p6 = subprocess.Popen(pigz_cmd, stdin=p5.stdout, stdout=outfile)
            p6.communicate()

        # Bowtie2 index and mapping
        subprocess.run(["bowtie2-build", "--threads", "20", "-f", ref, "bowtie2_index"])
        bowtie2_cmd = [
            "bowtie2", "-x", "bowtie2_index", "-p", "22", "-1", args.read1, "-2", args.read2,
            "--al-conc", f"map_{base_name}_R%.fastq", "--fr", "-I", "0", "-X", "1500", "-S", "/dev/null"
        ]
        subprocess.run(bowtie2_cmd)
        
        # Compress and move the Bowtie2 output
        subprocess.run(["pigz", "-p8", f"map_{base_name}_R1.fastq", f"map_{base_name}_R2.fastq"])
        subprocess.run(["mv", f"map_{base_name}_R1.fastq.gz", f"{args.sr_dir}/"])
        subprocess.run(["mv", f"map_{base_name}_R2.fastq.gz", f"{args.sr_dir}/"])
        subprocess.run(["rm", "-f", "bowtie2_index*"])

        # Run Flye
        flye_cmd = [
            "flye", "--pacbio-hifi", hifi_output, "--out-dir", f"{base_name}_flye",
            "--threads", "10", "--meta"
        ]
        with open(log_file, "a") as log:
            flye_proc = subprocess.Popen(flye_cmd, stdout=log, stderr=subprocess.STDOUT)
            flye_proc.communicate()
        subprocess.run(["cp", f"{base_name}_flye/assembly.fasta", f"{args.refined_bin_dir}/{base_name}_Flye.fa"])

        # Run SPAdes
        spades_cmd = [
            "spades.py", "-1", f"{args.sr_dir}/map_{base_name}_R1.fastq.gz", "-2", f"{args.sr_dir}/map_{base_name}_R2.fastq.gz",
            "--pacbio", hifi_output, "-o", f"{base_name}_spades-hybrid", "-t", "10"
        ]
        with open(log_file, "a") as log:
            spades_proc = subprocess.Popen(spades_cmd, stdout=log, stderr=subprocess.STDOUT)
            spades_proc.communicate()
        subprocess.run(["cp", f"{base_name}_spades-hybrid/scaffolds.fasta", f"{args.refined_bin_dir}/{base_name}_SPAdes-hybrid.fa"])

        # Run Unicycler
        unicycler_cmd = [
            "unicycler", "-1", f"{args.sr_dir}/map_{base_name}_R1.fastq.gz", "-2", f"{args.sr_dir}/map_{base_name}_R2.fastq.gz",
            "-l", hifi_output, "-o", f"{base_name}_unicycler-hybrid", "-t", "8"
        ]
        with open(log_file, "a") as log:
            unicycler_proc = subprocess.Popen(unicycler_cmd, stdout=log, stderr=subprocess.STDOUT)
            unicycler_proc.communicate()
        subprocess.run(["cp", f"{base_name}_unicycler-hybrid/assembly.fasta", f"{args.refined_bin_dir}/{base_name}_Unicycler-hybrid.fa"])

        print(f"Finished processing {base_name}")

if __name__ == "__main__":
    main()
