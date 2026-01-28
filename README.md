# Plant Gene Signature Builder

![Banner](assets/Banniere_PlantGeneSignatureBuilder.png)

<p align="center">
  <img src="assets/Logo_PlantGeneSignatureBuilder_small.png" alt="PGSB Logo" width="200"/>
</p>

**A production-ready Python toolkit for building scientifically defensible gene signatures using multi-evidence scoring**

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Overview

Plant Gene Signature Builder (PGSB) is a sophisticated pipeline for identifying and scoring candidate genes based on multiple independent evidence layers. Unlike simple GO-term extraction, PGSB integrates:

- **Capped multi-evidence scoring** (prevents GO dominance)
- **Synergy detection** (rewards convergent evidence)
- **Arabidopsis orthology** with annotation enrichment
- **Domain-based evidence** (Pfam/InterPro)
- **Run management** for reproducible analyses
- **QC diagnostics** and sensitivity analysis

## Latest Version: v3.1 - Production Release

### Key Features
- ✅ **Scientifically defensible scoring**: Caps per evidence layer prevent any single source from dominating
- ✅ **Run directory management**: Every analysis creates a unique, timestamped run with full traceability
- ✅ **Enriched outputs**: Overview tables include Arabidopsis symbols, annotations, and evidence summaries
- ✅ **Offline-only**: No web queries - fully reproducible with local resources
- ✅ **QC diagnostics**: Automated plots and sensitivity analysis
- ✅ **Production-ready**: Unit tests, manifest generation, SHA256 checksums

### Output Structure

Each run creates a complete directory:
```
results/
  runs/
    2026-01-29_000919__PCD_ROS_scored__5e7138e7__with_symbols/
      config_used.yaml          # Snapshot of config
      outputs/
        *_genes.base.txt        # All candidate IDs
        *_genes.scored.tsv      # Complete scoring table
        *_genes.ALL_overview.tsv    # ⭐ Enriched table (all genes)
        *_genes.HIGH_overview.tsv   # ⭐ Enriched table (high-confidence)
        *_genes.high_confidence.txt # Selected IDs
        *_genes.high_confidence.summary.txt
        manifest.json           # File checksums & metadata
      qc/                       # QC plots (with --qc flag)
      logs/
  latest -> runs/<most_recent>/  # Symlink to latest run
```

### What's in the Overview Tables?

The `ALL_overview.tsv` and `HIGH_overview.tsv` files provide rich annotations:

| Column | Description |
|--------|-------------|
| `brachy_gene_id` | Brachypodium gene ID (.v1.2) |
| `total_score` | Final score (sum of capped layers) |
| `go_score`, `domain_score`, `orth_score`, `tair_score`, `po_score` | Capped scores per layer |
| `synergy_bonus` | Bonus for convergent evidence |
| `go_terms` | Matched GO terms |
| `domain_hits` | Matched protein domains |
| `best_arabidopsis_id` | Top Arabidopsis ortholog (e.g., AT4G38360) |
| `best_orth_score` | InParanoid confidence score |
| `orth_class` | `one_to_one` / `one_to_many` / `none` |
| `arabidopsis_hits` | All Arabidopsis orthologs |
| `tair_symbol` | **Gene symbol** (e.g., LAZ1, EDS1) |
| `tair_annotation` | **Full description** from Araport11 GFF3 |
| `po_stage_hits` | Plant Ontology context matches |
| `evidence_summary` | Human-readable summary (e.g., "PCD_GO; 1:1_AT4G38360") |

**Example**: Gene `BdiBd21-3.1G0398400.v1.2` → Arabidopsis `AT4G38360` → Symbol `LAZ1` (LAZARUS 1)

---

## Versions

### v3.1 (Current) - Production Multi-Evidence Scoring
**Script**: `rank_gene_signatures.py`  
**Features**:
- Capped scoring system (go_max=15, domain_max=10, orthology_max=6, etc.)
- Synergy bonuses when independent evidence agrees
- Run directory management with config snapshots
- Enriched overview tables with Arabidopsis symbols/annotations
- Manifest generation with SHA256 checksums
- QC plots and diagnostics
- Unit tests (pytest)

### v2.0 (Archived) - Binary Filtering
**Script**: `build_pcd_list.py` (in `archives/v2.0/`)  
**Note**: Removed 91% of genes through hard filters. Replaced by scoring approach.

### v1.0 - GO Terms Only
**Script**: `create_list_from_GAF.py`  
**Use case**: Simple GO-based extraction for basic analyses

---

## Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/El-Castor/PlantGeneSignatureBuilder.git
cd PlantGeneSignatureBuilder

# Install dependencies
conda create -n pgsb_env python=3.9 pandas numpy pyyaml matplotlib seaborn pytest -y
conda activate pgsb_env

# Or with pip
pip install pandas numpy pyyaml matplotlib seaborn pytest
```

### v3.1 - Production Multi-Evidence Scoring (Recommended)

```bash
# Basic run
python rank_gene_signatures.py --config config_PCD_ROS_scoring.yaml

# With QC diagnostics
python rank_gene_signatures.py --config config_PCD_ROS_scoring.yaml --qc

# With custom run name
python rank_gene_signatures.py --config config_PCD_ROS_scoring.yaml --run_name my_analysis --qc

# With conda environment
conda run -n pgsb_env python rank_gene_signatures.py --config config_PCD_ROS_scoring.yaml --qc
```

**Outputs**: Check `results/latest/outputs/` for all files, including:
- `*_ALL_overview.tsv` - Enriched table with Arabidopsis symbols/annotations
- `*_HIGH_overview.tsv` - High-confidence genes with full annotations
- `manifest.json` - File checksums and metadata

### v1.0 - Simple GO Extraction

```bash
python create_list_from_GAF.py my_config.yaml
```

**Use case**: Quick GO-based gene lists without scoring

---

## Configuration

### Scoring System (v3.1)

The configuration file controls all aspects of scoring. Key sections:

Or with pip:
```bash
pip install pyyaml pandas
```

## Usage

### Basic Usage

```bash
# Use default config file (config_PCD_stress.yaml)
python create_list_from_GAF.py

# Use custom config file
python create_list_from_GAF.py my_custom_config.yaml
```

### With Conda Environment

```bash
conda run -n wot_env python create_list_from_GAF.py [config_file.yaml]
```

## Configuration File Format

Create a YAML file with the following structure:

```yaml
# Input/Output
input_go_file: "/path/to/GOannotation.tsv"
output_file: "output_genes.txt"

# Gene ID mapping
gene_id_pattern: "^(BdiBd21-3\\.\\dG\\d{7})"
gene_id_suffix: ".v1.2"

# GO filtering
go_level_filter: "BP"  # Options: BP, MF, CC, or all

# Target GO terms
go_terms:
  - go_id: "GO:0012501"
    description: "programmed cell death"
    category: "PCD"
  
  - go_id: "GO:0008219"
    description: "cell death"
    category: "PCD"
  
  # Add more GO terms as needed...
```

## Configuration Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `input_go_file` | Path to GO annotation TSV file | `"/path/to/GOannotation.tsv"` |
| `output_file` | Output filename for gene list | `"my_genes.txt"` |
| `gene_id_pattern` | Regex to extract core gene ID | `"^(BdiBd21-3\\.\\dG\\d{7})"` |
| `gene_id_suffix` | Suffix to append to gene IDs | `".v1.2"` |
| `go_level_filter` | GO level to keep (BP/MF/CC/all) | `"BP"` |
| `go_terms` | List of GO terms to extract | See example config |

## Example Workflows

### 1. PCD/Stress Genes (Default)

```bash
python create_list_from_GAF.py config_PCD_stress.yaml
```

### 2. Create Custom Gene List

Create a new config file:

```yaml
# config_photosynthesis.yaml
input_go_file: "/path/to/GOannotation.tsv"
output_file: "photosynthesis_genes.txt"
gene_id_pattern: "^(BdiBd21-3\\.\\dG\\d{7})"
gene_id_suffix: ".v1.2"
go_level_filter: "BP"

go_terms:
  - go_id: "GO:0015979"
    description: "photosynthesis"
    category: "Photosynthesis"
  
  - go_id: "GO:0009765"
    description: "photosynthesis, light harvesting"
    category: "Photosynthesis"
```

Then run:
```bash
python create_list_from_GAF.py config_photosynthesis.yaml
```

### 3. Defense Response Genes

```yaml
# config_defense.yaml
input_go_file: "/path/to/GOannotation.tsv"
output_file: "defense_genes.txt"
gene_id_pattern: "^(BdiBd21-3\\.\\dG\\d{7})"
gene_id_suffix: ".v1.2"
go_level_filter: "BP"

go_terms:
  - go_id: "GO:0006952"
    description: "defense response"
    category: "Defense"
  
  - go_id: "GO:0009607"
    description: "response to biotic stimulus"
    category: "Defense"
  
  - go_id: "GO:0002376"
    description: "immune system process"
    category: "Defense"
```

## Input File Format

The GO annotation file should be a TSV with these columns:

```
gene	GO	level
BdiBd21-3.3G0362100.1.p	GO:0008219	BP
BdiBd21-3.3G0362100.1.p	GO:0003824	MF
...
```

## Output

The script generates:
1. A gene list file (one gene ID per line)
2. Console output with statistics and gene counts per GO term

Example output:
```
================================================================================
GO Gene List Extraction Tool
================================================================================
Config file: config_PCD_stress.yaml

Configuration:
  Input file: /path/to/GOannotation.tsv
  Output file: PCD_WOT_genes.txt
  Gene ID suffix: .v1.2
  GO level filter: BP
  Number of GO terms: 10

Loading GO annotations...
  Loaded 123,456 annotations
Mapping protein IDs to gene IDs...
  Mapped 25,000 unique genes
Filtering for target GO terms...
  Found 5,234 annotations matching target GO terms
  Filtered to BP only: 3,456 annotations

✓ SUCCESS: 2,617 genes written to PCD_WOT_genes.txt

Gene count by GO term:
  GO:0006950: 1234 genes - response to stress
  GO:0009628:  987 genes - response to abiotic stimulus
  ...
================================================================================
```

## Tips

1. **Test with small config**: Start with 1-2 GO terms to verify the setup
2. **Check GO IDs**: Use [QuickGO](https://www.ebi.ac.uk/QuickGO/) to find GO term IDs
3. **Backup configs**: Keep different config files for different analyses
4. **Version control**: Track config files alongside your analysis scripts

## Troubleshooting

**No genes extracted**:
- Verify GO IDs exist in your annotation file: `grep "GO:0012501" GOannotation.tsv`
- Check the `go_level_filter` matches your data (BP/MF/CC)
- Verify the `gene_id_pattern` matches your gene ID format

**ID mapping issues**:
- Check gene ID format in your annotation file
- Adjust `gene_id_pattern` regex to match your ID format
- Test regex at [regex101.com](https://regex101.com/)

**Import error**:
```bash
pip install pyyaml pandas numpy
```

---

## Testing

Run the test suite:

```bash
# Run all tests
pytest tests/

# Run specific test file
pytest tests/test_run_manager.py

# With verbose output
pytest tests/ -v

# Smoke test with minimal data
pytest tests/test_smoke.py
```

---

## Run Management & Reproducibility

### Understanding Run Directories

Every execution creates a unique run directory with deterministic naming:

```
2026-01-29_000919__PCD_ROS_scored__5e7138e7__with_symbols
│
├── Timestamp: 2026-01-29_000919
├── Prefix: PCD_ROS_scored (from config)
├── Hash: 5e7138e7 (SHA1 of config content)
└── Custom: with_symbols (from --run_name flag)
```

This ensures:
- **No output collisions** between different analyses
- **Config traceability** via content hash
- **Easy identification** via timestamp and custom name
- **Reproducibility** via config snapshot

### Accessing Results

```bash
# Latest run (via symlink)
cat results/latest/outputs/*_HIGH_overview.tsv

# Specific run
cat results/runs/2026-01-29_000919__*/outputs/*_HIGH_overview.tsv

# List all runs
ls -lt results/runs/
```

### Manifest & Checksums

Each run includes `manifest.json` with SHA256 checksums:

```json
{
  "run_id": "2026-01-29_000919__PCD_ROS_scored__5e7138e7__with_symbols",
  "created": "2026-01-29T00:09:19.123456",
  "files": [
    {
      "path": "outputs/PCD_ROS_scored_genes.ALL_overview.tsv",
      "description": "Enriched overview table for all genes",
      "rows": 845,
      "size_bytes": 93702,
      "sha256": "cfce71f294f072ef6b8549c554787f7f6b6ba99b21f66ce20b7edfc4e6ceea9b"
    }
  ]
}
```

---

## Documentation

- **Quick start**: See above
- **v3.1 usage**: This README + `USER_GUIDE_v3_SCORING.md`
- **v2.0 usage (archived)**: See `archives/v2.0/USER_GUIDE_v2.md`
- **Architecture**: See `ARCHITECTURE_v3.md`
- **Publication checklist**: See `PUBLICATION_CHECKLIST.md`
- **Config templates**:
  - v3.1: `config_scoring_caps_template.yaml`, `config_PCD_ROS_scoring.yaml`
  - v1.0: `config_PCD_ROS.yaml`, `config_PCD_stress.yaml`
  - Custom signatures: `config_custom/*.yaml`

---

## Citation

If you use this tool in your research, please cite:

```bibtex
@software{pgsb2026,
  author = {Your Name},
  title = {Plant Gene Signature Builder: Multi-Evidence Scoring for Gene Signature Discovery},
  year = {2026},
  url = {https://github.com/El-Castor/PlantGeneSignatureBuilder},
  version = {3.1.0}
}
```

---

## License

MIT License - see LICENSE file for details

---

## Example Outputs

### v3.1 Overview Tables
```
PCD_ROS_genes.txt              # Gene list with categories
```

### v2.0
```
PCD_ROS_ortho_genes.base.txt                  # After GO extraction
PCD_ROS_ortho_genes.orthology_filtered.txt    # After orthology filter
PCD_ROS_ortho_genes.final.txt                 # Final output
PCD_ROS_ortho_orthology_filter_report.tsv     # Detailed report
PCD_ROS_ortho_filtering_summary.txt           # Summary statistics
```

## Citation

When using in publications, document all parameters used. See `PUBLICATION_CHECKLIST.md` for guidelines.

## Repository

https://github.com/El-Castor/PlantGeneSignatureBuilder

## bash
# Install PyYAML
conda install -n wot_env pyyaml
# or
pip install pyyaml
```

## License

Free to use and modify for research purposes.
