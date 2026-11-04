#!/usr/bin/env python3
"""
metagenome_16S_pipeline.py

A unified pipeline that:
1. Runs barrnap to find rRNA genes in a metagenome.
2. Uses bedtools to extract only 16S sequences.
3. Assigns nearest-species taxonomy using blastn against 16S_ribosomal_RNA.
4. Queries NCBI Entrez to resolve Phylum, Class, Order, and Family.
5. Outputs a TSV and a stacked barplot of the taxonomic composition.
"""

import argparse
import subprocess
import sys
import os
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt

OUTFMT = "6 qseqid sseqid pident length qlen qcovs evalue bitscore sscinames staxids stitle"
COLS = OUTFMT.split()[1:]

def species_from_hit(sciname, stitle):
    """Prefer BLAST's sscinames; fall back to first two words of the title[cite: 3]."""
    if sciname and sciname not in ("N/A", "NA", ""):
        return sciname
    words = stitle.split()
    return " ".join(words[:2]) if len(words) >= 2 else stitle

def rank_label(pident, qcov, a):
    """Assign taxonomic confidence based on Yarza et al. 2014 thresholds[cite: 3]."""
    if qcov < a.min_qcov:
        return "low_coverage"
    if pident >= a.species_id:
        return "species"
    if pident >= a.genus_id:
        return "genus"
    if pident >= a.family_id:
        return "family"
    return "unclassified_below_family"

def fetch_ncbi_lineage(taxids):
    """Fetch Phylum, Class, Order, Family for a list of NCBI TaxIDs via E-utilities API."""
    unique_taxids = list(set([t for t in taxids if t and t != "N/A" and t != "NA"]))
    if not unique_taxids:
        return {}
    
    print(f"Fetching lineages for {len(unique_taxids)} unique TaxIDs from NCBI...")
    taxids_str = ",".join(unique_taxids)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={taxids_str}&retmode=xml"
    
    try:
        req = urllib.request.Request(url)
        response = urllib.request.urlopen(req).read()
        root = ET.fromstring(response)
        
        tax_dict = {}
        for taxon in root.findall('Taxon'):
            tid = taxon.find('TaxId').text
            lineage = {'phylum': 'Unknown', 'class': 'Unknown', 'order': 'Unknown', 'family': 'Unknown'}
            lineage_ex = taxon.find('LineageEx')
            if lineage_ex is not None:
                for t in lineage_ex.findall('Taxon'):
                    rank = t.find('Rank').text
                    name = t.find('ScientificName').text
                    if rank in lineage:
                        lineage[rank] = name
            tax_dict[tid] = lineage
        return tax_dict
    except Exception as e:
        print(f"Warning: Failed to fetch taxonomy from NCBI: {e}", file=sys.stderr)
        return {}

def plot_taxonomy(df, out_prefix):
    """Generate a stacked bar chart for Phylum, Class, Order, and Family."""
    if df.empty:
        return

    ranks = ['phylum', 'class', 'order', 'family']
    fig, axes = plt.subplots(1, 4, figsize=(18, 6))
    fig.suptitle('16S rRNA Taxonomic Composition', fontsize=16)

    for ax, rank in zip(axes, ranks):
        if rank in df.columns:
            counts = df[rank].value_counts().to_frame().T
            counts.plot(kind='bar', stacked=True, ax=ax)
            ax.set_title(rank.capitalize())
            ax.set_xticks([])
            ax.set_xlabel("Assembly Sequences")
            # Move legend outside the plot frame
            ax.legend(title=rank.capitalize(), bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    plot_file = f"{out_prefix}_tax_barplot.pdf"
    plt.savefig(plot_file, bbox_inches='tight')
    print(f"Saved taxonomic barplots to {plot_file}")

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, help="Input metagenome FASTA")
    ap.add_argument("-o", "--output", default="assembly_16S", help="Output prefix")
    ap.add_argument("--db", default="16S_ribosomal_RNA", help="BLAST db name/path")
    ap.add_argument("-t", "--threads", type=int, default=8)
    ap.add_argument("--max-hits", type=int, default=50)
    ap.add_argument("--species-id", type=float, default=98.65)
    ap.add_argument("--genus-id", type=float, default=94.5)
    ap.add_argument("--family-id", type=float, default=86.5)
    ap.add_argument("--min-qcov", type=float, default=80.0, help="min query coverage ")
    ap.add_argument("--tie-window", type=float, default=0.2, help="list species within this  identity of best hit[cite: 3]")
    a = ap.parse_args()

    out_fasta = f"{a.output}.16S.fasta"
    out_tsv = f"{a.output}.taxonomy.tsv"

    # Step 1: Run barrnap
    print("Running Barrnap to identify rRNAs...")
    barrnap_gff = f"{a.output}.barrnap.gff"
    with open(barrnap_gff, "w") as f_out:
        subprocess.run(["barrnap", "--kingdom", "bac", "--threads", str(a.threads), a.input], 
                       stdout=f_out, stderr=subprocess.DEVNULL)
    
    # Step 2: Extract only 16S regions
    print("Extracting 16S sequences...")
    filtered_gff = f"{a.output}.16S.gff"
    with open(barrnap_gff, "r") as f_in, open(filtered_gff, "w") as f_out:
        for line in f_in:
            if line.startswith("#") or "Name=16S_rRNA" in line:
                f_out.write(line)
    
    subprocess.run(["bedtools", "getfasta", "-fi", a.input, "-bed", filtered_gff, "-fo", out_fasta])
    
    if not os.path.exists(out_fasta) or os.path.getsize(out_fasta) == 0:
        sys.exit("No 16S rRNA genes found in the assembly. Exiting.")

    # Step 3: Run BLAST[cite: 3]
    print("Running BLASTn against 16S database...")
    blast_tmp = f"{a.output}.blast.tmp"
    cmd = ["blastn", "-query", out_fasta, "-db", a.db, "-task", "megablast",
           "-max_target_seqs", str(a.max_hits), "-num_threads", str(a.threads),
           "-evalue", "1e-10", "-outfmt", OUTFMT, "-out", blast_tmp]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"blastn failed:\n{r.stderr}")

    hits = defaultdict(list)
    with open(blast_tmp) as tmp:
        for line in tmp:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11: continue
            q = f[0]
            d = dict(zip(COLS, f))
            # Handle N/A floats safely
            d["pident"] = float(d["pident"]) if d["pident"] != "N/A" else 0.0
            d["qcovs"] = float(d["qcovs"]) if d["qcovs"] != "N/A" else 0.0
            d["bitscore"] = float(d["bitscore"]) if d["bitscore"] != "N/A" else 0.0
            d["species"] = species_from_hit(d["sscinames"], d["stitle"])
            hits[q].append(d)

    qids = []
    with open(out_fasta) as fh:
        for l in fh:
            if l.startswith(">"):
                qids.append(l[1:].split()[0])

    # Step 4: Resolve Lineage & Output Results
    print("Parsing results and assigning lineage...")
    results = []
    taxids_to_fetch = [hits[q][0]["staxids"].split(";")[0] for q in qids if q in hits and hits[q]]
    lineage_map = fetch_ncbi_lineage(taxids_to_fetch)

    for q in qids:
        hs = hits.get(q)
        if not hs:
            results.append({"query_id": q, "nearest_species": "NO_HIT", "taxid": "NA", 
                            "pident": "NA", "qcov": "NA", "confidence": "no_hit",
                            "phylum": "Unknown", "class": "Unknown", "order": "Unknown", "family": "Unknown"})
            continue
        
        hs.sort(key=lambda x: (-x["bitscore"], -x["pident"]))
        best = hs[0]
        ties = sorted({h["species"] for h in hs if h["pident"] >= best["pident"] - a.tie_window} - {best["species"]})
        
        tid = best["staxids"].split(";")[0]
        lin = lineage_map.get(tid, {'phylum': 'Unknown', 'class': 'Unknown', 'order': 'Unknown', 'family': 'Unknown'})
        
        results.append({
            "query_id": q,
            "nearest_species": best["species"],
            "taxid": tid,
            "pident": f'{best["pident"]:.2f}',
            "qcov": f'{best["qcovs"]:.0f}',
            "confidence": rank_label(best["pident"], best["qcovs"], a),
            "tied_species": ";".join(ties) if ties else "NA",
            "phylum": lin['phylum'],
            "class": lin['class'],
            "order": lin['order'],
            "family": lin['family']
        })

    df = pd.DataFrame(results)
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"Results written to {out_tsv}")

    # Step 5: Plotting
    plot_taxonomy(df, a.output)

    # Cleanup temporary files
    for temp_file in [barrnap_gff, filtered_gff, blast_tmp]:
        if os.path.exists(temp_file):
            os.remove(temp_file)

if __name__ == "__main__":
    main()
