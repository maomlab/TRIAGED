#!/usr/bin/env python3
"""
Run `analyize_boltz_predictions.py` for every receptor found under
`input_files/DUDEZ_benchmark/` using the corresponding
`{RECEPTOR}_combined_ids.csv` as `--reference-data`.

This script discovers receptors, prepares an output directory per receptor,
and invokes the analyzer script. Use `--dry-run` to print the commands
without executing them.
"""
from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path
import os


def find_receptor_csvs(dudez_dir: Path):
    return sorted(dudez_dir.glob("*_combined_ids.csv"))


def main():
    repo_root = Path(__file__).resolve().parents[2]
    default_dudez = repo_root / "input_files" / "DUDEZ_benchmark"
    default_outputs = repo_root / "outputs"
    analyzer = repo_root / "modules" / "triage_boltz" / "analyize_boltz_predictions.py"

    parser = argparse.ArgumentParser(description="Batch analyze Boltz predictions with DUDEZ reference labels")
    parser.add_argument("--dudez-dir", type=Path, default=default_dudez, help="Path to input_files/DUDEZ_benchmark directory")
    parser.add_argument("--outputs-dir", type=Path, default=default_outputs, help="Path to outputs directory containing receptor result folders")
    parser.add_argument("--analyzer-script", type=Path, default=analyzer, help="Path to analyize_boltz_predictions.py script")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running them")
    parser.add_argument("--compute-metrics", action="store_true", help="Pass -m to the analyzer to compute metrics")
    parser.add_argument("--bootstrap", action="store_true", help="Pass -b to the analyzer to run bootstrap metrics")
    args = parser.parse_args()

    if not args.dudez_dir.exists():
        print(f"DUDEZ benchmark dir not found: {args.dudez_dir}")
        sys.exit(1)
    if not args.analyzer_script.exists():
        print(f"Analyzer script not found: {args.analyzer_script}")
        sys.exit(1)

    csvs = find_receptor_csvs(args.dudez_dir)
    if not csvs:
        print(f"No '*_combined_ids.csv' files found in {args.dudez_dir}")
        sys.exit(0)

    print(f"Found {len(csvs)} receptors to process.")
    for csv in csvs:
        receptor = csv.name.replace("_combined_ids.csv", "")
        input_dir = args.outputs_dir / receptor
        if not input_dir.exists():
            print(f"[WARN] outputs directory for receptor '{receptor}' not found at {input_dir}; skipping.")
            continue

        out_dir = input_dir / "dudez_analysis"
        out_dir.mkdir(parents=True, exist_ok=True)

        cmd = [sys.executable, str(args.analyzer_script), "-i", str(input_dir), "-o", str(out_dir), "--reference-data", str(csv)]
        if args.compute_metrics:
            cmd.insert(-2, "-m")
        if args.bootstrap:
            # insert before the reference-data and output args
            cmd.insert(-2, "-b")

        print(f"Processing receptor: {receptor}")
        print(" ", " ".join(cmd))
        if args.dry_run:
            continue

        try:
            completed = subprocess.run(cmd, check=False)
            if completed.returncode != 0:
                print(f"[ERROR] Analyzer exited with code {completed.returncode} for receptor {receptor}")
        except KeyboardInterrupt:
            print("Interrupted by user")
            sys.exit(1)
        except Exception as e:
            print(f"[ERROR] Failed to run analyzer for {receptor}: {e}")


if __name__ == "__main__":
    main()
