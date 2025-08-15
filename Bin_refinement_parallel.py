#!/usr/bin/env python
import os
import argparse
import subprocess
import sys
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed

def run_cmd(cmd, quiet=False, **kwargs):
    """Run a command with optional quiet mode and handle errors."""
    if quiet:
        kwargs.setdefault("stdout", subprocess.DEVNULL)
        kwargs.setdefault("stderr", subprocess.DEVNULL)
    try:
        subprocess.run(cmd, check=True, **kwargs)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed: {' '.join(cmd)}", file=sys.stderr)
        return False
    return True

def run_cmd_pipe(cmds, output_file=None, quiet=False):
    """Run piped commands with optional quiet mode."""
    processes = []
    prev_stdout = None
    try:
        for i, cmd in enumerate(cmds):
            stdout_target = output_file if (i == len(cmds) - 1 and output_file) else subprocess.PIPE
            stderr_target = subprocess.DEVNULL if quiet else None
            p = subprocess.Popen(cmd, stdin=prev_stdout, stdout=stdout_target, stderr=stderr_target)
            if prev_stdout:
                prev_stdout.close()
            prev_stdout = p.stdout if stdout_target == subprocess.PIPE else None
            processes.append(p)
        for p in processes:
            p.wait()
            if p.returncode != 0:
                raise subprocess.CalledProcessError(p.returncode, p.args)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Pipeline failed at: {e.cmd}", file=sys.stderr)
        return False
    return True

def process_bin(ref, args, bin_threads):
    try:
        base_name = os.path.splitext(ref)[0]
        hifi_tmp = f"{args.hifi_mapped_dir}/{base_name}.tmp.fastq.gz"
        hifi_output = f"{args.hifi_mapped_dir}/{base_name}_map-Hifi_mapped.fastq.gz"
        log_file = f"log_{base_name}.txt"

        # Index reference
        run_cmd(["samtools", "faidx", ref], quiet=True)

        # HiFi mapping
        ok = run_cmd_pipe(
            [
                ["minimap2", "-ax", args.hifi_preset, "-t", str(bin_threads), ref, args.pacbio],
                ["samtools", "view", f"-@{max(1,bin_threads//5)}", "-F4", "-h", "-", f"-q{args.mapq}"],
                ["samclip", "--max", "100", "--ref", ref, "-"],
                ["samtools", "sort", f"-@{max(1,bin_threads//5)}", "-T", f"tmp_sort_{base_name}"],
                ["samtools", "fastq"],
                ["pigz", f"-p{max(2,bin_threads//2)}"],
            ],
            output_file=open(hifi_tmp, "wb"),
            quiet=True
        )
        if not ok or not os.path.exists(hifi_tmp) or os.path.getsize(hifi_tmp) < 1024:
            print(f"[WARN] HiFi mapping failed or too few reads for {base_name}")
            if os.path.exists(hifi_tmp):
                os.remove(hifi_tmp)
            return None
        os.replace(hifi_tmp, hifi_output)

        # Bowtie2 mapping (short reads)
        idx_prefix = f"bowtie2_index_{base_name}"
        run_cmd(["bowtie2-build", "--threads", str(bin_threads), "-f", ref, idx_prefix], quiet=True)
        bowtie2_cmd = [
            "bowtie2", "-x", idx_prefix, "-p", str(bin_threads),
            "-1", args.read1, "-2", args.read2,
            "--al-conc", f"map_{base_name}_R%.fastq", "--fr", "-I", "0", "-X", "1500", "-S", "/dev/null"
        ]
        run_cmd(bowtie2_cmd, quiet=True)
        run_cmd(["pigz", f"-p{max(2,bin_threads//2)}", f"map_{base_name}_R1.fastq", f"map_{base_name}_R2.fastq"], quiet=True)
        run_cmd(["mv", f"map_{base_name}_R1.fastq.gz", f"{args.sr_dir}/"])
        run_cmd(["mv", f"map_{base_name}_R2.fastq.gz", f"{args.sr_dir}/"])
        for f in os.listdir("."):
            if f.startswith(idx_prefix):
                os.remove(f)

        # --- Flye assembly with fallback ---
        flye_outdir = f"{base_name}_flye"
        with open(log_file, "a") as log:
            # Try HiFi mode first
            flye_cmd_hifi = ["flye", "--pacbio-hifi", hifi_output, "--out-dir", flye_outdir,
                             "--threads", str(max(1, bin_threads//2)), "--meta"]
            print(f"[INFO] Running Flye HiFi for {base_name}")
            run_cmd(flye_cmd_hifi, stdout=log, stderr=subprocess.STDOUT)

        if os.path.exists(f"{flye_outdir}/assembly.fasta"):
            run_cmd(["cp", f"{flye_outdir}/assembly.fasta", f"{args.refined_bin_dir}/{base_name}_Flye.fa"])
        else:
            print(f"[WARN] Flye HiFi failed for {base_name}, retrying with CLR mode...")
            shutil.rmtree(flye_outdir, ignore_errors=True)
            with open(log_file, "a") as log:
                flye_cmd_clr = ["flye", "--pacbio-raw", hifi_output, "--out-dir", flye_outdir,
                                "--threads", str(max(1, bin_threads//2)), "--meta"]
                run_cmd(flye_cmd_clr, stdout=log, stderr=subprocess.STDOUT)
            if os.path.exists(f"{flye_outdir}/assembly.fasta"):
                run_cmd(["cp", f"{flye_outdir}/assembly.fasta", f"{args.refined_bin_dir}/{base_name}_Flye.fa"])
            else:
                print(f"[ERROR] Flye failed in both HiFi and CLR modes for {base_name}")

        shutil.rmtree(flye_outdir, ignore_errors=True)

        # SPAdes
        spades_cmd = [
            "spades.py", "-1", f"{args.sr_dir}/map_{base_name}_R1.fastq.gz", "-2", f"{args.sr_dir}/map_{base_name}_R2.fastq.gz",
            "--pacbio", hifi_output, "-o", f"{base_name}_spades-hybrid", "-t", str(max(1, bin_threads//2))
        ]
        with open(log_file, "a") as log:
            run_cmd(spades_cmd, stdout=log, stderr=subprocess.STDOUT)
        if os.path.exists(f"{base_name}_spades-hybrid/scaffolds.fasta"):
            run_cmd(["cp", f"{base_name}_spades-hybrid/scaffolds.fasta", f"{args.refined_bin_dir}/{base_name}_SPAdes-hybrid.fa"])
        shutil.rmtree(f"{base_name}_spades-hybrid", ignore_errors=True)

        # Unicycler
        unicycler_cmd = [
            "unicycler", "-1", f"{args.sr_dir}/map_{base_name}_R1.fastq.gz", "-2", f"{args.sr_dir}/map_{base_name}_R2.fastq.gz",
            "-l", hifi_output, "-o", f"{base_name}_unicycler-hybrid", "-t", str(max(1, bin_threads//2))
        ]
        with open(log_file, "a") as log:
            run_cmd(unicycler_cmd, stdout=log, stderr=subprocess.STDOUT)
        if os.path.exists(f"{base_name}_unicycler-hybrid/assembly.fasta"):
            run_cmd(["cp", f"{base_name}_unicycler-hybrid/assembly.fasta", f"{args.refined_bin_dir}/{base_name}_Unicycler-hybrid.fa"])
        shutil.rmtree(f"{base_name}_unicycler-hybrid", ignore_errors=True)

        print(f"[INFO] Finished processing {base_name}")
        return base_name

    except Exception as e:
        print(f"[ERROR] Unexpected error with {ref}: {e}", file=sys.stderr)
        return None

def main():
    parser = argparse.ArgumentParser(description="Run mapping and assembly pipeline.")
    parser.add_argument("--read1", required=True)
    parser.add_argument("--read2", required=True)
    parser.add_argument("--pacbio", required=True)
    parser.add_argument("--hifi_mapped_dir", default="Hifi_mapped")
    parser.add_argument("--sr_dir", default="sr")
    parser.add_argument("--refined_bin_dir", default="refined_bin")
    parser.add_argument("--mapq", type=int, default=20)
    parser.add_argument("--hifi_preset", default="map-hifi")
    parser.add_argument("--ext", default=".fa")
    parser.add_argument("--jobs", type=int, default=2, help="Number of bins to process in parallel")
    parser.add_argument("--max_total_threads", type=int, default=24, help="Max total threads used")
    args = parser.parse_args()

    # Create dirs
    os.makedirs(args.hifi_mapped_dir, exist_ok=True)
    os.makedirs(args.sr_dir, exist_ok=True)
    os.makedirs(args.refined_bin_dir, exist_ok=True)

    # Find bins
    bin_files = [f for f in os.listdir(".") if f.startswith("bin") and f.endswith(args.ext)]
    if not bin_files:
        print(f"[WARNING] No bin files with extension {args.ext} found.")
        sys.exit(0)

    bin_threads = max(1, args.max_total_threads // args.jobs)
    print(f"[INFO] Processing {len(bin_files)} bins with {args.jobs} jobs, {bin_threads} threads per bin")

    results = []
    with ProcessPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(process_bin, ref, args, bin_threads): ref for ref in bin_files}
        for future in as_completed(futures):
            results.append(future.result())

    # Seqkit stats
    try:
        run_cmd(["seqkit", "stats", "-T"] + bin_files, quiet=True, stdout=open("bin_original_stats.tsv", "w"))
        refined_files = [os.path.join(args.refined_bin_dir, f) for f in os.listdir(args.refined_bin_dir) if f.endswith(".fa") or f.endswith(".fna")]
        if refined_files:
            run_cmd(["seqkit", "stats", "-T"] + refined_files, quiet=True, stdout=open("bin_refined_stats.tsv", "w"))
        print("[INFO] Saved seqkit stats for original and refined bins.")
    except Exception as e:
        print(f"[ERROR] Failed to generate seqkit stats: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()

