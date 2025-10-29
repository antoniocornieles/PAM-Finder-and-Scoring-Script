#!/usr/bin/env python3
"""
PAM Finder & Scoring Script
---------------------------
Finds PAM (Protospacer Adjacent Motif) sequences in a DNA string,
optionally fetched from NCBI or read from a local FASTA file.

Now includes JSON and CSV export support.

Usage:
    python3 pam_filter.py --input local.fasta --out results.json
    python3 pam_filter.py --accession NC_000913.3 --out results.csv
"""

import argparse
import csv
import json
from Bio import SeqIO, Entrez
from tqdm import tqdm
import plotly.graph_objs as go
import numpy as np

# --- CONFIG ---
Entrez.email = "your_email@example.com"  # Replace with your email for NCBI access
PAM_PATTERN = "NGG"  # Cas9 PAM sequence pattern


def fetch_sequence_from_ncbi(accession):
    """Fetch a DNA sequence from NCBI by accession number."""
    print(f"Fetching {accession} from NCBI...")
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    print(f"Retrieved: {record.id} ({len(record.seq)} bp)")
    return str(record.seq).upper()


def read_sequence_from_fasta(filepath):
    """Read a local FASTA file and return its sequence as a string."""
    record = SeqIO.read(filepath, "fasta")
    print(f"Loaded: {record.id} ({len(record.seq)} bp)")
    return str(record.seq).upper()


def pam_score(context):
    """
    Simple PAM scoring:
    - +10 if GG matches perfectly
    - +5 if N is a valid nucleotide (A/T/C/G)
    - + bonus if local GC content > 50%
    """
    score = 0
    if context[1:] == "GG":
        score += 10
    if context[0] in "ATCG":
        score += 5
    gc_content = sum(1 for c in context if c in "GC") / len(context)
    if gc_content > 0.5:
        score += 2
    return score


def find_pam_sites(sequence, pam_pattern="NGG"):
    """Find all PAM sites in the given DNA sequence."""
    results = []
    print(f"Scanning sequence for PAM pattern: {pam_pattern}")

    for i in tqdm(range(len(sequence) - len(pam_pattern) + 1)):
        block = sequence[i:i+len(pam_pattern)]
        if block[1:] == "GG":  # match Cas9 PAM
            guide_start = max(0, i - 20)
            guide_seq = sequence[guide_start:i]
            results.append({
                "position": i,
                "context": block,
                "guide": guide_seq,
                "score": pam_score(block)
            })
    return results


def export_results(results, filepath):
    """Save PAM results as JSON or CSV based on file extension."""
    if filepath.endswith(".json"):
        with open(filepath, "w") as f:
            json.dump(results, f, indent=4)
        print(f"✅ Results saved to {filepath}")
    elif filepath.endswith(".csv"):
        with open(filepath, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"✅ Results saved to {filepath}")
    else:
        print("⚠️ Unsupported file format. Use .json or .csv")

def visualize_sequence_3d(sequence, pam_sites):
    # Convert sequence to numeric values for visualization
    base_map = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    seq_numeric = [base_map.get(base, np.nan) for base in sequence]

    x = list(range(len(seq_numeric)))
    y = seq_numeric
    z = [0.1 if i in pam_sites else 0 for i in range(len(sequence))]

    fig = go.Figure(data=[go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers',
        marker=dict(
            size=4,
            color=z,
            colorscale='Viridis',
            opacity=0.8
        )
    )])

    fig.update_layout(
        title='3D CRISPR PAM Site Visualization',
        scene=dict(
            xaxis_title='Position',
            yaxis_title='Base (A=0, T=1, C=2, G=3)',
            zaxis_title='PAM Score / Presence'
        )
    )

    fig.show()



def main():
    parser = argparse.ArgumentParser(description="Find PAM sequences in DNA.")
    parser.add_argument("--input", help="Path to local FASTA file")
    parser.add_argument("--accession", help="NCBI accession ID")
    parser.add_argument("--out", help="Output file (.json or .csv)")
    args = parser.parse_args()

    # Load DNA
    if args.input:
        seq = read_sequence_from_fasta(args.input)
    elif args.accession:
        seq = fetch_sequence_from_ncbi(args.accession)
    else:
        raise ValueError("Provide either --input or --accession")

    # Analyze sequence
    pam_sites = find_pam_sites(seq)
    print(f"\nFound {len(pam_sites)} PAM sites.\n")

    # Print top 10
    for site in pam_sites[:10]:
        print(f"Pos {site['position']:>8} | PAM: {site['context']} | Score: {site['score']} | Guide: {site['guide'][-10:]}")

    visualize_sequence_3d(seq, pam_sites)

    # Export
    if args.out and pam_sites:
        export_results(pam_sites, args.out)


if __name__ == "__main__":
    main()
