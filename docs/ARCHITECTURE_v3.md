# PlantGeneSignatureBuilder - System Architecture

## Overview

Multi-evidence gene signature scoring system for plant genomics research.  
**Mode: OFFLINE** - All annotations from local files only.

---

## System Components

### 1. Core Scripts

| Script | Purpose | Mode |
|---|---|---|
| `create_list_from_GAF.py` | GO-based gene extraction | Simple filtering |
| `build_pcd_list.py` | v2.0 - Binary filtering pipeline | DEPRECATED |
| **`rank_gene_signatures.py`** | **v3.0 - Multi-evidence scoring** | **CURRENT** |

### 2. Configuration System

All behavior controlled via YAML configs:

| Template | Purpose |
|---|---|
| `config_scoring_template.yaml` | Basic scoring (GO + orthology) |
| `config_scoring_full_template.yaml` | Complete pipeline (all evidence layers) |
| `config_PCD_ROS_scoring.yaml` | Example: PCD/ROS genes |

### 3. Resource Directory

```
resources/
├── arabidopsis/
│   ├── Araport11_TAIR10.1_symbol_description.gff3
│   └── Arabidopsis_annotation.txt (PO terms)
├── domains/
│   ├── brachypodium_pfam.tsv
│   └── brachypodium_interpro.tsv
└── README.md (download instructions)
```

---

## Scoring Framework

### Evidence Layers

```
┌─────────────────────────────────────────────┐
│  GO ANNOTATIONS (required)                  │
│  ↓                                           │
│  Base gene list (510 genes)                 │
└─────────────────────────────────────────────┘
                    ↓
    ┌───────────────────────────────────┐
    │   MULTI-EVIDENCE SCORING          │
    │   (all genes retained)            │
    ├───────────────────────────────────┤
    │ 1. GO score         (2-60 pts)    │
    │ 2. Domain score     (0-12 pts)    │
    │ 3. Orthology score  (0-5 pts)     │
    │ 4. Keyword score    (0-6 pts)     │
    │ 5. PO context score (0-4 pts)     │
    │                                   │
    │ TOTAL SCORE:        2-65 pts      │
    └───────────────────────────────────┘
                    ↓
    ┌───────────────────────────────────┐
    │  SELECTION BY DISTRIBUTION        │
    │  - Quantile (top 10%)             │
    │  - Knee detection (automatic)     │
    └───────────────────────────────────┘
                    ↓
            High-confidence genes (64)
```

### Score Calculation

For each gene:
```python
total_score = (
    go_score +           # GO term weights
    domain_score +       # Domain matches × weight
    orthology_score +    # Best hit × multiplier + 1:1 bonus
    keyword_score +      # Keyword matches (capped)
    po_context_score     # Senescence context (capped)
)
```

---

## Data Flow

### Input Files

1. **GO annotations** (TSV)
   - Columns: gene (protein ID), GO, level
   - Example: `BdiBd21-3.1G0001100.1.p`, `GO:0012501`, `BP`

2. **Orthology file** (TSV - inParanoid format)
   - Columns: OrtoA, OrtoB
   - Example: `BdistachyonBd21-3:BdiBd21-3.1G0001100.1 1.000`, `Athalianacolumbia:AT1G01010.1 1.000`

3. **Domain annotations** (TSV - configurable)
   - User-defined columns
   - Example: gene, pfam_id, interpro_id

4. **Arabidopsis GFF3** (GFF3 format)
   - Columns: seqid, source, type, start, end, score, strand, phase, attributes
   - Required: type="gene", attributes with symbol, curator_summary, etc.

5. **Arabidopsis PO annotations** (TSV)
   - Columns: locus_name, term_name, evidence
   - Example: `AT1G01010`, `sporophyte senescent stage`, `IEP`

### Output Files

1. **`<prefix>_genes.base.txt`**
   - One gene ID per line
   - All GO-matched genes (unfiltered)

2. **`<prefix>_genes.scored.tsv`**
   - Tab-separated table
   - 18 columns with all evidence and scores
   - Sorted by total_score (descending)

3. **`<prefix>_genes.high_confidence.txt`**
   - Top genes by score distribution
   - Selection method in summary file

4. **`<prefix>_genes.high_confidence.summary.txt`**
   - Selection statistics
   - Score distribution info

---

## ID Mapping

### Brachypodium Gene IDs

```
Protein ID: BdiBd21-3.1G0001100.1.p
     ↓ (regex: protein_to_core_regex)
Core ID: BdiBd21-3.1G0001100
     ↓ (add: seurat_suffix)
Seurat ID: BdiBd21-3.1G0001100.v1.2
```

### Arabidopsis Gene IDs

```
Orthology file: Athalianacolumbia:AT1G01010.1 1.000
     ↓ (regex: arabidopsis_id_regex_extract)
Locus: AT1G01010
```

---

## Offline Mode Implementation

### GFF3 Parser
- Reads gene features only (NOT mRNA/CDS)
- Extracts: symbol, curator_summary, full_name, computational_description
- Priority: curator_summary > full_name > comp_desc > Note
- URL-decodes common entities (%2C, %20, etc.)

### PO Parser
- Reads TSV with locus_name and term_name
- Detects senescence keywords:
  - senescence, senescent, cell death, hypersensitive, dying, necrosis
- Stores PO terms and senescence hits per gene

### Domain Parser
- Fully configurable file paths and column names
- Supports multiple databases (Pfam, InterPro, etc.)
- Matches expected domains against annotations
- Scores by number of matches with optional cap

---

## Comparison: v2.0 vs v3.0

| Feature | v2.0 (Binary Filtering) | v3.0 (Scoring) |
|---|---|---|
| **Philosophy** | Pass/fail filters | Additive evidence |
| **Gene loss** | High (filtering cascade) | None (all retained) |
| **Orthology** | Required (hard filter) | Contributes to score |
| **Domain** | Optional hard filter | Contributes to score |
| **Output** | Filtered list | Ranked list + scores |
| **Flexibility** | Low (on/off only) | High (weight tuning) |
| **Transparency** | Opaque (genes discarded) | Clear (all evidence visible) |

### Example

**Input:** 510 genes from GO terms

**v2.0 Output:**
- Base: 510 genes
- After domain filter: 156 genes (354 lost)
- After orthology filter: 45 genes (111 lost)
- **Final: 45 genes** (465 excluded, 91% loss)

**v3.0 Output:**
- Base: 510 genes (all retained)
- Scored: 510 genes (full evidence table)
- High-confidence: 64 genes (top 12.5% by score)
- **All genes accessible in scored table**

---

## Configuration Best Practices

### 1. Start Simple
```yaml
# Baseline: GO only
evidence:
  orthology_evidence:
    enabled: false
```

### 2. Add Orthology
```yaml
# GO + conservation
evidence:
  orthology_evidence:
    enabled: true
```

### 3. Full Pipeline
```yaml
# All evidence layers
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

### 4. Tune Weights

Increase domain importance:
```yaml
domain_database:
  scoring:
    match_score: 5  # default: 3
    max_domain_score: 20  # default: 12
```

Decrease keyword importance:
```yaml
score_weights:
  arabidopsis_keywords:
    keyword_hit: 0.5  # default: 1
```

---

## Validation

### Test Cases

1. **Gene with all evidence**
   - GO: programmed cell death (3 pts)
   - Domain: metacaspase (9 pts)
   - Orthology: 1:1 AT ortholog (5 pts)
   - Keyword: "cell death" hit (1 pt)
   - PO: senescence stage (2 pts)
   - **Total: 20 pts** ✓ High confidence

2. **Gene with GO + orthology only**
   - GO: oxidative stress (2 pts)
   - Orthology: 1:1 AT ortholog (5 pts)
   - **Total: 7 pts** ✓ Medium confidence

3. **Gene with GO only**
   - GO: cell death (3 pts)
   - **Total: 3 pts** ✓ Low confidence

**All three genes retained**, ranked by evidence strength.

---

## Future Extensions

### Planned Features
- Expression data integration (RNA-seq scores)
- Co-expression network scores
- Phylogenetic conservation depth
- GWAS association scores

### Extensibility

Add new evidence layer:
```python
def _score_custom_evidence(self, genes):
    """Score custom evidence layer"""
    for gene in genes:
        custom_score = compute_custom_metric(gene)
        self.gene_scores[gene]['custom_score'] = custom_score
```

Update config:
```yaml
score_weights:
  custom_evidence:
    weight: 2
    max_points: 10
```

---

## Citation and Reproducibility

### Publication Statement

> "Gene signatures were ranked using a multi-evidence scoring framework 
> integrating Gene Ontology annotations, protein domain evidence 
> (Pfam/InterPro), orthology to Arabidopsis thaliana, curated functional 
> annotations (Araport11), and Plant Ontology expression context. All 
> annotations were retrieved exclusively from local resources (offline mode) 
> to ensure reproducibility. Score weights: [specify weights used]."

### Version Control

Tag releases:
```bash
git tag -a v3.0 -m "Multi-evidence scoring system"
git push origin v3.0
```

Archive resources:
```bash
tar -czf resources_v3.0.tar.gz resources/
```

---

## Support

### Documentation
- `USER_GUIDE_v3_SCORING.md` - User guide
- `resources/README.md` - Resource download instructions
- `config_scoring_full_template.yaml` - Full configuration template

### Example Configs
- `config_PCD_ROS_scoring.yaml` - PCD/ROS genes
- `config_scoring_template.yaml` - Basic template

### Contact
Repository: https://github.com/El-Castor/PlantGeneSignatureBuilder
