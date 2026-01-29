# Plant Gene Signature Builder (PGSB) v3.1

Production-ready multi-evidence gene signature scoring pipeline with reproducible run management.

## Overview

PGSB scores genes using multiple independent lines of evidence:
- **GO term annotations** (with scoring caps to prevent dominance)
- **Protein domain matches** (Pfam, InterPro)
- **Orthology evidence** (InParanoid confidence scores)
- **Arabidopsis TAIR annotations** (keyword matching in gene descriptions)
- **PO (Plant Ontology) context** (developmental stage/tissue annotations)
- **Synergy bonuses** (rewards convergent evidence from independent sources)

All evidence layers are scored, capped, and combined. Every run is tracked in a unique directory with complete provenance.

---

## Quick Start

```bash
# Basic run
conda run -n wot_env python rank_gene_signatures.py --config config_PCD_ROS_scoring.yaml

# With QC diagnostic plots
conda run -n wot_env python rank_gene_signatures.py --config config.yaml --qc

# With custom run name
conda run -n wot_env python rank_gene_signatures.py --config config.yaml --run_name experiment1

# Allow overwriting existing run (use with caution)
conda run -n wot_env python rank_gene_signatures.py --config config.yaml --overwrite
```

---

## Run Directory Structure

Every execution creates a unique, timestamped run directory:

```
results/
├── latest -> runs/2026-01-28_234815__PCD_ROS_scored__1e8a360e__test1/  (symlink)
└── runs/
    └── 2026-01-28_234815__PCD_ROS_scored__1e8a360e__test1/
        ├── config_used.yaml           # Config snapshot for reproducibility
        ├── logs/                       # Log files (future use)
        ├── qc/                         # QC plots and diagnostics
        │   ├── *_go_vs_domain.png
        │   ├── *_total_vs_go.png
        │   ├── *_score_distribution.png
        │   └── *_qc_summary.txt
        └── outputs/
            ├── *_genes.base.txt                    # IDs only (all GO-derived genes)
            ├── *_genes.scored.tsv                  # Complete scored table
            ├── *_genes.ALL_overview.tsv            # ⭐ Enriched ALL genes
            ├── *_genes.high_confidence.txt         # IDs only (selected subset)
            ├── *_genes.HIGH_overview.tsv           # ⭐ Enriched high-conf genes
            ├── *_genes.high_confidence.summary.txt # Selection stats
            └── manifest.json                       # File checksums & metadata
```

### Run Directory Naming

Format: `YYYY-MM-DD_HHMMSS__<prefix>__<config_hash>__<optional_run_name>`

- **Timestamp**: Exact execution time
- **Prefix**: From `config['output_prefix']`
- **Config hash**: First 8 chars of SHA1 hash of config content (ensures reproducibility)
- **Run name**: Optional custom suffix (e.g., `--run_name test1`)

### Finding Your Results

```bash
# Latest run (via symlink)
ls results/latest/outputs/

# Specific run
ls results/runs/2026-01-28_234815__PCD_ROS_scored__1e8a360e__test1/outputs/
```

---

## Output Files Explained

### 1. `*_genes.base.txt`
Plain text list of gene IDs (one per line). All genes matching GO terms.

```
BdiBd21-3.1G0398400.v1.2
BdiBd21-3.3G0279600.v1.2
...
```

### 2. `*_genes.scored.tsv`
Complete scored table with all evidence layers (raw and capped scores).

| brachy_gene_id | go_score_raw | go_score_capped | domain_score_raw | ... | total_score |
|---|---|---|---|---|---|
| BdiBd21-3.1G0398400.v1.2 | 24 | 15 | 0 | ... | 20.0 |

### 3. ⭐ `*_genes.ALL_overview.tsv` (NEW!)
**Enriched overview table for ALL genes** with aggregated annotations:

| Column | Description |
|--------|-------------|
| `brachy_gene_id` | Brachypodium gene ID (.v1.2 suffix) |
| `total_score` | Final combined score (capped) |
| `go_score`, `domain_score`, ... | Individual layer scores (capped) |
| `synergy_bonus` | Bonus points from convergent evidence |
| `go_terms` | Comma-separated list of matched GO terms |
| `domain_hits` | Comma-separated domain IDs (e.g., "Pfam:PF00931") |
| `best_arabidopsis_id` | Best ortholog (e.g., AT3G48090) |
| `best_orth_score` | InParanoid confidence score (0-1) |
| `orth_class` | one_to_one, one_to_many, or none |
| `arabidopsis_hits` | All Arabidopsis orthologs (comma-separated) |
| `tair_symbol` | Arabidopsis gene symbol |
| `tair_annotation` | Gene description from TAIR/Araport11 |
| `po_stage_hits` | Plant Ontology developmental stage matches |
| `evidence_summary` | Human-readable summary (e.g., "PCD_GO+Domain; 1:1_AT3G48090") |

**Example row:**
```
brachy_gene_id: BdiBd21-3.4G0416100.v1.2
total_score: 20.0
go_score: 15
domain_score: 0
orth_score: 5.0
go_terms: GO:0008219, GO:0009626, GO:0012501
best_arabidopsis_id: AT3G48090
orth_class: one_to_one
tair_annotation: Encodes a hypersensitive response marker
evidence_summary: PCD+ROS_GO; 1:1_AT3G48090
```

### 4. ⭐ `*_genes.HIGH_overview.tsv` (NEW!)
Same format as ALL_overview.tsv, but only high-confidence genes.

### 5. `*_genes.high_confidence.txt`
Plain text list of high-confidence gene IDs (subset of base list).

### 6. `manifest.json`
Provenance metadata with SHA256 checksums for all files:

```json
{
  "run_id": "2026-01-28_234815__PCD_ROS_scored__1e8a360e__test1",
  "created": "2026-01-28T23:48:18.696116",
  "files": [
    {
      "path": "outputs/PCD_ROS_scored_genes.base.txt",
      "description": "Base gene list (IDs only)",
      "rows": 845,
      "size_bytes": 21125,
      "sha256": "0434253a034dd87b5e89327f88d6909a..."
    },
    ...
  ]
}
```

---

## Configuration

See `config_scoring_caps_template.yaml` for full template with all options.

### Key Sections

#### Scoring Caps (prevents any layer from dominating)
```yaml
scoring:
  caps:
    go_max: 15        # Max points from GO terms
    domain_max: 10
    orthology_max: 6
    tair_max: 20
    po_max: 4
    synergy_max: 6
```

#### Synergy Rules (rewards convergent evidence)
```yaml
scoring:
  synergy_rules:
    - name: "PCD_GO+domain"
      if_all:
        - "has_PCD_GO"
        - "has_expected_domain"
      bonus: 3
```

---

## Testing

```bash
# Run all tests
conda run -n wot_env python -m pytest tests/ -v

# Run unit tests only
conda run -n wot_env python -m pytest tests/test_run_manager.py -v

# Run smoke test with minimal data
conda run -n wot_env python -m pytest tests/test_smoke.py -v
```

---

## Scientific Rationale

### Why Capped Scoring?

**Problem**: GO annotations are abundant. A gene with many GO matches can score 60+ points from GO alone, dominating genes with strong domain or orthology evidence.

**Solution**: Cap each evidence layer (e.g., `go_max=15`). Forces reliance on multiple lines of evidence.

**Result**: Top-scoring genes have:
- High GO relevance (but capped at 15)
- Domain matches (up to 10 points)
- Strong orthology (up to 6 points)
- TAIR keyword support (up to 20 points)
- Synergy bonuses (up to 6 points)

### Why Synergy Bonuses?

Independent evidence sources agreeing is more reliable than any single source. Examples:
- Gene has PCD GO terms **AND** expected protein domain → +3 bonus
- Gene has PCD GO terms **AND** TAIR keywords like "hypersensitive response" → +2 bonus

---

## Offline Mode

All resources are local. No web queries to:
- TAIR/Araport
- UniProt
- Pfam/InterPro APIs

Required local files:
- GO annotations (TSV)
- InParanoid orthology table (TSV)
- Araport11 GFF3 (optional, for TAIR annotations)
- PO annotations (optional, for developmental context)
- Domain annotations (optional, for Pfam/InterPro)

---

## Example Workflow

```bash
# 1. Prepare config for your gene signature
cp config_scoring_caps_template.yaml config_my_signature.yaml
# Edit: set GO terms, expected domains, keywords, etc.

# 2. Run analysis with QC
conda run -n wot_env python rank_gene_signatures.py \
    --config config_my_signature.yaml \
    --qc \
    --run_name initial_run

# 3. Inspect results
cd results/latest/outputs/
head -20 my_signature_genes.HIGH_overview.tsv

# 4. Check QC plots
open ../qc/*_score_distribution.png

# 5. Verify reproducibility
cat ../config_used.yaml
cat manifest.json
```

---

## Migration from v2.0

v2.0 used **binary filtering** (genes were dropped if they lacked domains/orthologs).

v3.1 uses **additive scoring** (all GO-derived genes retained, scored by available evidence).

**Migrating:**
1. Update config to v3.1 format (add `scoring.caps` section)
2. Adjust GO weights (lower values since capping is now enforced)
3. Run and compare HIGH_overview.tsv with old output

v2.0 files archived in `archives/v2.0/`

---

## Citation & License

When using this pipeline, please cite the source repository and include:
- Run ID (from directory name)
- Config hash (ensures exact reproducibility)
- Manifest checksums (file integrity)

---

## Support

For issues, feature requests, or questions:
- GitHub: https://github.com/El-Castor/PlantGeneSignatureBuilder
- Check `USER_GUIDE_v3_SCORING.md` for detailed documentation
- See `ARCHITECTURE_v3.md` for system design
