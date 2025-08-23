#!/usr/bin/env python3
import os
import argparse
import subprocess
import csv

def fasta_length(fasta_file):
    """Calculate total sequence length of a FASTA file"""
    length = 0
    with open(fasta_file) as f:
        seq = []
        for line in f:
            if line.startswith(">"):
                if seq:
                    length += len("".join(seq))
                    seq = []
            else:
                seq.append(line.strip())
        if seq:
            length += len("".join(seq))
    return length

def count_variants(vcf, vartype):
    """Count variants of a specific type (snps or indels) in VCF"""
    p1 = subprocess.Popen(
        ["bcftools", "view", "-H", "-v", vartype, vcf],
        stdout=subprocess.PIPE
    )
    p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    out, _ = p2.communicate()
    return int(out.decode().strip())

def analyze_bin(bin_fasta, r1, r2, outdir, threads=8):
    """Map reads, call variants, and compute SNP/INDEL diversity"""
    bin_name = os.path.basename(bin_fasta).replace(".fa", "").replace(".fasta", "")

    sam = os.path.join(outdir, f"{bin_name}.sam")
    bam = os.path.join(outdir, f"{bin_name}.bam")
    sorted_bam = os.path.join(outdir, f"{bin_name}.sorted.bam")
    vcf = os.path.join(outdir, f"{bin_name}.vcf.gz")

    # minimap2 mapping
    subprocess.run([
        "minimap2", "-ax", "sr", "-t", str(threads),
        "-o", sam, bin_fasta, r1, r2
    ], check=True)

    # samtools sort & index
    subprocess.run(["samtools", "view", "-bS", sam], stdout=open(bam, "wb"), check=True)
    subprocess.run(["samtools", "sort", "-@", str(threads), "-o", sorted_bam, bam], check=True)
    subprocess.run(["samtools", "index", sorted_bam], check=True)

    # bcftools mpileup & call
    mpileup = subprocess.Popen([
        "bcftools", "mpileup", "-Ou", "-f", bin_fasta, sorted_bam
    ], stdout=subprocess.PIPE)
    call = subprocess.Popen([
        "bcftools", "call", "-mv", "-Oz", "-o", vcf
    ], stdin=mpileup.stdout)
    mpileup.stdout.close()
    call.communicate()
    subprocess.run(["bcftools", "index", vcf], check=True)

    # count SNPs and INDELs
    snp_count = count_variants(vcf, "snps")
    indel_count = count_variants(vcf, "indels")

    # genome size
    genome_size = fasta_length(bin_fasta)

    # calculate ratios
    snp_ratio = snp_count / genome_size if genome_size > 0 else 0
    indel_ratio = indel_count / genome_size if genome_size > 0 else 0

    return {
        "bin": bin_name,
        "snp_count": snp_count,
        "indel_count": indel_count,
        "genome_size": genome_size,
        "snp_ratio": snp_ratio,
        "indel_ratio": indel_ratio
    }

def main():
    parser = argparse.ArgumentParser(description="Assess allelic diversity (SNPs & INDELs) per bin.")
    parser.add_argument("--sr_dir", required=True, help="Directory containing short reads (map_bin.*_R1.fastq.gz, map_bin.*_R2.fastq.gz)")
    parser.add_argument("--ref_dir", required=True, help="Directory containing bin fasta files")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads for mapping")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    results = []
    for ref in sorted(os.listdir(args.ref_dir)):
        if ref.endswith(".fa") or ref.endswith(".fasta"):
            bin_name = ref.replace(".fa", "").replace(".fasta", "")
            r1 = os.path.join(args.sr_dir, f"map_{bin_name}_R1.fastq.gz")
            r2 = os.path.join(args.sr_dir, f"map_{bin_name}_R2.fastq.gz")
            if os.path.exists(r1) and os.path.exists(r2):
                res = analyze_bin(os.path.join(args.ref_dir, ref), r1, r2, args.outdir, threads=args.threads)
                results.append(res)

    # write summary CSV
    out_csv = os.path.join(args.outdir, "allelic_diversity_summary.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "bin", "genome_size",
            "snp_count", "snp_ratio",
            "indel_count", "indel_ratio"
        ])
        writer.writeheader()
        writer.writerows(results)

if __name__ == "__main__":
    main()
