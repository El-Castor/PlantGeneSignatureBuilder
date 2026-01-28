# Plant Gene Signature Scorer - USER GUIDE

## Multi-Evidence Scoring System (OFFLINE MODE)

**Version:** 3.0 (Scoring Framework)  
**Mode:** OFFLINE - All annotations from local files only

---

## Overview

This system scores plant gene signatures using **convergent evidence** from multiple layers:

1. **GO annotations** (base requirement)
2. **Protein domains** (Pfam/InterPro) - mechanistic evidence
3. **Orthology** to Arabidopsis - evolutionary conservation
4. **Arabidopsis functional annotations** (GFF3) - curated knowledge
5. **PO/expression context** - developmental/senescence stages

**NO GENES ARE FILTERED OUT** - all GO-matched genes are retained and ranked by total score.

---

## Key Principle: Additive Scoring, Not Filtering

| Evidence Layer | Role | Effect |
|---|---|---|
| **GO terms** | Defines candidate universe | Required (defines base list) |
| **Protein domains** | Mechanistic plausibility | Increases score |
| **Orthology** | Evolutionary support | Increases score |
| **AT annotation** | Functional context | Increases score |
| **PO context** | Developmental relevance | Increases score |

**Absence of any single evidence does NOT exclude a gene.**

---

## Outputs

### 1. Base Gene List
**File:** `<prefix>_genes.base.txt`  
All genes matching GO terms (unfiltered, unranked).

### 2. Scored Table
**File:** `<prefix>_genes.scored.tsv`

Columns:
- `brachy_gene_id` - Brachypodium gene ID
- `go_terms` - Matched GO terms
- `go_score` - GO evidence score
- `matched_domains` - Protein domains (Pfam/InterPro)
- `domain_score` - Domain evidence score
- `arabidopsis_orthologs` - Top 3 AT orthologs
- `best_ortholog` - Best AT ortholog
- `best_ortholog_score` - Orthology confidence
- `orthology_class` - one_to_one / one_to_many / none
- `orthology_score` - Orthology evidence score
- `at_symbol` - Arabidopsis symbol
- `at_annotation` - Arabidopsis functional annotation
- `keyword_hits` - Matched keywords
- `keyword_score` - Keyword evidence score
- `po_terms` - Plant Ontology terms
- `po_context_hits` - Senescence/PCD PO hits
- `po_context_score` - PO evidence score
- **`total_score`** - Combined evidence score

### 3. High-Confidence List
**File:** `<prefix>_genes.high_confidence.txt`  
Top genes selected by score distribution (quantile or knee detection).

**File:** `<prefix>_genes.high_confidence.summary.txt`  
Selection method and statistics.

---

## Usage

### Basic Command
```bash
python rank_gene_signatures.py --config config_scoring.yaml
```

### Example Workflow

1. **GO-only scoring** (baseline):
```yaml
output_prefix: "PCD_baseline"
evidence:
  domain_database:
    enabled: false
  orthology_evidence:
    enabled: false
  arabidopsis_context:
    enabled: false
  po_context:
    enabled: false
```

2. **GO + Orthology**:
```yaml
output_prefix: "PCD_with_ortho"
evidence:
  orthology_evidence:
    enabled: true
    orthology_file: "path/to/inParanoid.tsv"
```

3. **Full multi-evidence** (all layers):
```yaml
output_prefix: "PCD_full"
domain_database:
  enabled: true
evidence:
  orthology_evidence:
    enabled: true
  arabidopsis_context:
    enabled: true
  po_context:
    enabled: true
```

---

## Configuration Guide

### Score Weights

All weights are configurable:

```yaml
score_weights:
  go_core_terms:
    GO:0012501: 3    # programmed cell death
    GO:0008219: 3    # cell death
  
  go_ros_terms:
    GO:0006979: 2    # oxidative stress
  
  orthology:
    best_hit_score_multiplier: 3
    one_to_one_bonus: 2
  
  arabidopsis_keywords:
    keyword_hit: 1
    max_keyword_points: 6
  
  po_context:
    senescence_hit: 2
    max_po_points: 4
```

### Domain Database (MANDATORY for mechanistic evidence)

**Fully configurable** - define ALL domain annotation files:

```yaml
domain_database:
  enabled: true
  
  annotation_files:
    - path: "resources/domains/brachypodium_pfam.tsv"
      database: "Pfam"
      gene_id_column: "gene"
      domain_column: "pfam_id"
    
    - path: "resources/domains/brachypodium_interpro.tsv"
      database: "InterPro"
      gene_id_column: "gene"
      domain_column: "interpro_id"
  
  expected_domains:
    Pfam:
      - "PF00656"   # Caspase
      - "PF14743"   # Metacaspase
      - "PF00112"   # Cysteine protease
    InterPro:
      - "IPR029058"  # Caspase-like
  
  scoring:
    match_score: 3
    max_domain_score: 12
```

### Arabidopsis Context (OFFLINE MODE)

```yaml
evidence:
  arabidopsis_context:
    enabled: true
    
    # Local GFF3 file (REQUIRED)
    gff3_file: "resources/arabidopsis/Araport11_TAIR10.1_symbol_description.gff3"
    
    # Keywords to search in annotations
    keywords:
      - "cell death"
      - "PCD"
      - "metacaspase"
      - "senescence"
      - "autophagy"
      - "ROS"
      - "peroxidase"
```

### PO Context (OFFLINE MODE)

```yaml
evidence:
  po_context:
    enabled: true
    po_annotation_file: "resources/arabidopsis/Arabidopsis_annotation.txt"
```

Automatically detects:
- senescence
- sporophyte senescent stage
- cell death
- hypersensitive
- dying, necrosis

### Selection Strategy

**Quantile mode** (recommended):
```yaml
selection:
  mode: "quantile"
  quantile: 0.90    # Top 10%
```

**Knee detection** (automatic):
```yaml
selection:
  mode: "knee"
  min_genes: 150
  max_genes: 350
```

---

## Score Interpretation

### Example Scores

| Gene | GO | Domain | Orthology | Keywords | PO | **Total** | Confidence |
|---|---|---|---|---|---|---|---|
| Gene_A | 3 | 9 | 5 | 3 | 2 | **22** | High |
| Gene_B | 3 | 0 | 5 | 0 | 0 | **8** | Medium |
| Gene_C | 2 | 0 | 0 | 0 | 0 | **2** | Low |

**All three genes are retained**, but Gene_A has the strongest multi-evidence support.

### Typical Score Distributions

| Evidence Layer | Min | Max | Mean |
|---|---|---|---|
| GO score | 2 | 60 | 5-10 |
| Domain score | 0 | 12 | 0-3 |
| Orthology score | 0 | 5 | 0-3 |
| Keyword score | 0 | 6 | 0-2 |
| PO context score | 0 | 4 | 0-1 |
| **Total score** | 2 | 65 | 8-15 |

---

## Offline Mode Guarantee

**NO external database queries are performed.**

All annotations come from:
- Local GO annotation files
- Local domain annotation files (InterProScan/Pfam)
- Local Arabidopsis GFF3 (Araport11)
- Local PO annotation files

This ensures:
✓ Reproducibility  
✓ Independence from database availability  
✓ Compliance with offline/air-gapped requirements  
✓ Version-controlled annotations

---

## Required Resources

See [resources/README.md](resources/README.md) for:
- Arabidopsis GFF3 download instructions
- Arabidopsis PO annotation sources
- Domain annotation generation (InterProScan)
- File format specifications

---

## Comparison: Filtering vs Scoring

### OLD Approach (v2.0 - Binary Filtering)
```
510 genes → [domain filter] → 156 genes → [orthology filter] → 45 genes
```
**Problem:** Loses genes with partial evidence.

### NEW Approach (v3.0 - Additive Scoring)
```
510 genes → ALL retained → scored by evidence → top 64 genes (90th percentile)
```
**Advantage:** Genes with strong GO + orthology but no domain info are still ranked high.

---

## Example Use Cases

### Use Case 1: PCD Gene Discovery
**Goal:** Find plant programmed cell death genes.

**Config:**
- GO terms: PCD core (GO:0012501, GO:0008219)
- Domains: caspase, metacaspase, cysteine protease
- Keywords: "cell death", "PCD", "apoptosis"
- Selection: top 10% by score

### Use Case 2: ROS/Oxidative Stress
**Goal:** Identify ROS metabolism genes.

**Config:**
- GO terms: oxidative stress, ROS metabolic process
- Domains: peroxidase, NADPH oxidase
- Keywords: "ROS", "oxidative", "peroxidase"
- PO context: senescence stages

### Use Case 3: Mechanistic Focus
**Goal:** Prioritize genes with structural evidence.

**Config:**
- Increase domain score weight to 5
- Set keyword weight to 0.5
- Require min domain_score > 3 for high-confidence

---

## Citation

If using this pipeline in publications, cite:

**Arabidopsis annotations:**
- Cheng CY et al. (2017) Araport11: a complete reannotation of the Arabidopsis thaliana reference genome. Plant Cell 29:9-10

**Domain databases:**
- Jones P et al. (2014) InterProScan 5. Bioinformatics 30:1236-1240
- Mistry J et al. (2021) Pfam: The protein families database in 2021. Nucleic Acids Res 49:D412-D419

**Methodology:**
Include in methods:
> "Gene signatures were scored using a multi-evidence framework integrating GO annotations, protein domain evidence (Pfam/InterPro), orthology to Arabidopsis, and curated functional annotations. All evidence was retrieved exclusively from local annotation resources (offline mode)."

---

## Troubleshooting

### No domain scores?
- Check `domain_database.enabled: true`
- Verify annotation file paths exist
- Check column names match config

### No orthology scores?
- Verify orthology file path
- Check regex patterns extract correct IDs
- Run with debug: inspect ortho_data dictionary

### Low keyword scores?
- Check GFF3 file loaded successfully
- Verify keywords are in Arabidopsis annotations
- Try broader keywords

### All genes have same score?
- Enable more evidence layers
- Adjust score weights
- Check GO term weights are varied

---

## Advanced: Custom Scoring Functions

Future extensions could add:
- Expression-based scoring (RNA-seq data)
- Phylogenetic conservation scores
- Network centrality scores
- GWAS association scores

All following the same principle: **additive evidence, no hard filtering**.
