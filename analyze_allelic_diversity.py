#!/usr/bin/env python3
import os
import argparse
import subprocess
import glob
import csv

def fasta_length(fasta_file):
    """FASTA の総塩基数を返す"""
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

def analyze_bin(bin_fasta, r1, r2, outdir, threads=4):
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

    # SAM -> sorted BAM
    subprocess.run(["samtools", "view", "-@",
                    str(threads), "-bS", sam, "-o", bam], check=True)
    subprocess.run(["samtools", "sort", "-@",
                    str(threads), "-o", sorted_bam, bam], check=True)
    subprocess.run(["samtools", "index", sorted_bam], check=True)

    # Variant calling
    mpileup = subprocess.Popen(
        ["bcftools", "mpileup", "-Ou", "-f", bin_fasta, sorted_bam],
        stdout=subprocess.PIPE)
    call = subprocess.Popen(
        ["bcftools", "call", "-mv", "-Oz", "-o", vcf],
        stdin=mpileup.stdout)
    mpileup.stdout.close()
    call.communicate()
    subprocess.run(["bcftools", "index", vcf], check=True)

    # count SNPs
    p1 = subprocess.Popen(["bcftools", "view", "-H", "-v", "snps", vcf],
                          stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    out, _ = p2.communicate()
    snp_count = int(out.decode().strip())

    # genome size
    genome_size = fasta_length(bin_fasta)
    ratio = snp_count / genome_size if genome_size > 0 else 0

    return {"bin": bin_name, "snp_count": snp_count,
            "genome_size": genome_size, "snp_ratio": ratio}

def main():
    parser = argparse.ArgumentParser(
        description="Estimate allelic diversity per bin using minimap2 + bcftools")
    parser.add_argument("--sr_dir", required=True, help="Directory with short reads (map_bin.X_R1.fastq.gz / R2)")
    parser.add_argument("--ref_dir", required=True, help="Directory with bin FASTA files")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--threads", type=int, default=4, help="Threads for minimap2/samtools")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    results = []
    for ref in sorted(glob.glob(os.path.join(args.ref_dir, "*.fa")) +
                      glob.glob(os.path.join(args.ref_dir, "*.fasta"))):
        bin_name = os.path.basename(ref).replace(".fa", "").replace(".fasta", "")
        r1 = os.path.join(args.sr_dir, f"map_bin.{bin_name.split('.')[-1]}_R1.fastq.gz")
        r2 = os.path.join(args.sr_dir, f"map_bin.{bin_name.split('.')[-1]}_R2.fastq.gz")

        if not (os.path.exists(r1) and os.path.exists(r2)):
            print(f"[WARNING] Skipping {bin_name}, reads not found")
            continue

        print(f"[INFO] Processing {bin_name} ...")
        res = analyze_bin(ref, r1, r2, args.outdir, threads=args.threads)
        results.append(res)

    # save results
    out_csv = os.path.join(args.outdir, "allelic_diversity_summary.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["bin", "snp_count", "genome_size", "snp_ratio"])
        writer.writeheader()
        writer.writerows(results)

    print(f"[INFO] Done. Results written to {out_csv}")

if __name__ == "__main__":
    main()
