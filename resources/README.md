# Resources Directory - Local Annotation Files

This directory contains **local annotation files** used for gene signature scoring.

**MODE: OFFLINE ONLY**  
No queries to external databases (TAIR, Araport, InterPro, Pfam, etc.) are performed.

---

## Required Resources

### 1. Arabidopsis Genome Annotation (GFF3)

**File:** `arabidopsis/Araport11_TAIR10.1_symbol_description.gff3`

**Source:**  
Download from TAIR/Araport11:
- https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/

**Contains:**
- Arabidopsis gene IDs (ATxGxxxxx)
- Official symbols (symbol=)
- Full gene names (full_name=)
- Computational descriptions
- Curator summaries (curator_summary=)

**Usage:**  
Provides functional context for Arabidopsis orthologs.  
Keyword matching on curator summaries and gene names.

---

### 2. Arabidopsis PO Annotation

**File:** `arabidopsis/Arabidopsis_annotation.txt`

**Source:**  
Download from TAIR:
- https://www.arabidopsis.org/download_files/Ontologies/Plant_Ontology/

**Contains:**
- Plant Ontology (PO) terms
- Stage/tissue labels
- Senescence-related developmental stages

**Usage:**  
Detects genes expressed in senescence or cell death contexts.

---

### 3. Protein Domain Annotations (MANDATORY)

**Files:**  
- `domains/brachypodium_pfam.tsv`
- `domains/brachypodium_interpro.tsv`

**Source:**  
Generate locally using:
- **Pfam:** InterProScan or hmmscan against Pfam-A.hmm
- **InterPro:** InterProScan

**Expected format (TSV):**
```
gene    pfam_id    pfam_name    ...
BdiBd21-3.1G0001100.1.p    PF00112    Peptidase_C1    ...
```

**Usage:**  
Domain evidence is the PRIMARY mechanistic filter.  
Expected domains (PCD-related):
- PF00656 (Caspase domain)
- PF14743 (Metacaspase)
- PF00112 (Cysteine protease)
- PF00141 (Peroxidase)
- IPR029058 (Caspase-like)
- IPR011600 (Autophagy)

---

## How to Populate This Directory

### Option 1: Download Pre-Curated Files (Recommended)

Contact your bioinformatics core or check lab shared resources.

### Option 2: Generate Locally

**Arabidopsis GFF3:**
```bash
wget https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz
gunzip Araport11_GFF3_genes_transposons.201606.gff.gz
mv Araport11_GFF3_genes_transposons.201606.gff resources/arabidopsis/Araport11_TAIR10.1_symbol_description.gff3
```

**Arabidopsis PO annotations:**
```bash
wget https://www.arabidopsis.org/download_files/Ontologies/Plant_Ontology/po_annotations.txt
mv po_annotations.txt resources/arabidopsis/Arabidopsis_annotation.txt
```

**Brachypodium domain annotations (InterProScan):**
```bash
# Install InterProScan
# Run on Brachypodium proteome
interproscan.sh -i brachypodium_proteins.fasta -f tsv -o resources/domains/brachypodium_interpro.tsv
```

---

## Citation

If using these resources in publications, cite:

- **Araport11:** Cheng CY et al. (2017) Plant Cell 29:9-10
- **TAIR:** Berardini TZ et al. (2015) Nucleic Acids Res 43:D1009-14
- **InterProScan:** Jones P et al. (2014) Bioinformatics 30:1236-1240
- **Pfam:** Mistry J et al. (2021) Nucleic Acids Res 49:D412-D419

---

## Offline Mode Guarantee

**No external database queries are performed.**  
All functional annotations are retrieved from local files only.

This ensures:
- Reproducibility
- Independence from database availability
- Compliance with offline/air-gapped computing requirements
