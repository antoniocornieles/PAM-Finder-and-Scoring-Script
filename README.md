# 🧬 PAM Finder & Scoring Script

A lightweight CRISPR/Bioinformatics tool for identifying **PAM (Protospacer Adjacent Motif)** sequences in DNA.  
Supports **FASTA file input** or **direct NCBI accession retrieval**, with **JSON/CSV export** and a simple **PAM scoring system** for ranking potential CRISPR/Cas9 target sites.

---

## 📖 Overview

**PAM (Protospacer Adjacent Motif)** sequences are short DNA patterns recognized by the **Cas9** enzyme during CRISPR gene editing.  
Cas9 scans the genome for these motifs (typically `NGG`), binds to them, and makes precise cuts guided by the **guide RNA (gRNA)** sequence.

This tool helps identify and rank those PAM sites in any DNA sequence.

---

## ⚙️ Features

- 🔍 Detects PAM sequences (default pattern: `NGG`)  
- 🧠 Assigns a **score** based on PAM context and GC content  
- 📁 Reads from **local FASTA** files or fetches from **NCBI databases**  
- 💾 Exports results as **JSON** or **CSV**  
- 🧪 Designed for easy integration with future CRISPR design tools

---

## 🚀 Usage

### **1. Local FASTA Input**
```bash
python3 pam_filter.py --input example.fasta --out results.json
