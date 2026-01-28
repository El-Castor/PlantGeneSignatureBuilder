# Plant Gene Signature Builder v2.0 - User Guide

## Overview

Version 2.0 extends the GO term-based gene extraction with optional filtering modules:
- **Base extraction**: GO term-based gene list (v1.0 functionality)
- **Domain filtering**: Filter by Pfam/InterPro domains (optional)
- **Orthology filtering**: Filter by Arabidopsis orthologues (optional)

## Quick Start

### Basic Usage (GO terms only)

```bash
# Use existing v1.0 configs
conda run -n wot_env python create_list_from_GAF.py config_PCD_ROS.yaml
```

### Advanced Usage (with filtering)

```bash
# v2.0 with optional filters
conda run -n wot_env python build_pcd_list.py --config config_extended.yaml
```

## Installation

Requirements already met if you used v1.0:
```bash
conda install -n wot_env pyyaml pandas
```

## Configuration

### Minimal Config (v1.0 compatible)

Use existing configs with `create_list_from_GAF.py`:
- `config_PCD_ROS.yaml`
- `config_PCD_stress.yaml`
- `config_cell_division.yaml`
- `config_photosynthesis.yaml`

### Extended Config (v2.0)

Create a YAML file with these sections:

```yaml
output_prefix: "my_genes"
input_go_file: "/path/to/GOannotation.tsv"

id_mapping:
  protein_to_core_regex: "^(BdiBd21-3\\.\\dG\\d{7})"
  seurat_suffix: ".v1.2"

go_level_filter: "BP"

go_terms:
  - go_id: "GO:0012501"
    description: "programmed cell death"
    category: "PCD"
  # ... more GO terms

filters:
  domain_filter:
    enabled: true/false
    # ... domain settings
  
  orthology_filter:
    enabled: true/false
    # ... orthology settings
```

See `config_extended_template.yaml` for full documentation.

## Pipeline Stages

### Stage 1: Base GO Extraction (Always)

Extracts genes annotated with specified GO terms.

**Output**: `<prefix>_genes.base.txt`

**Example**:
```yaml
go_terms:
  - go_id: "GO:0012501"
    description: "programmed cell death"
    category: "PCD"
```

### Stage 2: Domain Filtering (Optional)

Keeps only genes with specified Pfam/InterPro domains.

**Output**: `<prefix>_genes.domain_filtered.txt`, `<prefix>_domain_filter_report.tsv`

**Config**:
```yaml
filters:
  domain_filter:
    enabled: true
    annotation_file: "/path/to/domains.tsv"
    id_column: "gene"
    domain_columns: ["pfam", "interpro"]
    allowed_domains:
      pfam: ["PF00931", "PF12613"]
      interpro: ["IPR002182", "IPR001611"]
```

**Report columns**:
- `gene_id`: Gene identifier
- `matched_domain_type`: pfam or interpro
- `matched_domain_id`: Domain ID that matched

### Stage 3: Orthology Filtering (Optional)

Filters based on Arabidopsis orthologues.

**Output**: `<prefix>_genes.orthology_filtered.txt`, `<prefix>_orthology_filter_report.tsv`

**Config**:
```yaml
filters:
  orthology_filter:
    enabled: true
    orthology_file: "/path/to/inParanoid_orthology.tsv"
    require_arabidopsis_in_set: false
    arabidopsis_pcd_gene_list: null  # optional
    prefer_one_to_one: true
    one_to_one_rule:
      max_hits: 1
      min_best_score: 0.8
      min_second_best_gap: 0.2
```

**Parameters**:
- `prefer_one_to_one`: If true, keep only 1:1 orthologues
- `require_arabidopsis_in_set`: If true, require Arabidopsis orthologue in reference list
- `arabidopsis_pcd_gene_list`: File with Arabidopsis gene IDs (one per line)
- `one_to_one_rule`:
  - `max_hits`: Max Arabidopsis hits for 1:1 classification
  - `min_best_score`: Minimum confidence score
  - `min_second_best_gap`: Required score gap if multiple hits

**Report columns**:
- `brachy_gene_id`: Brachypodium gene (.v1.2)
- `brachy_protein_id`: Original protein ID
- `arabidopsis_hits`: Comma-separated Arabidopsis IDs
- `best_hit_score`: Highest confidence score
- `n_hits`: Number of Arabidopsis orthologues
- `classification`: one_to_one / one_to_many / low_confidence
- `in_reference_set`: Whether orthologue is in reference set (if provided)

### Stage 4: Final Output (Always)

**Output**: `<prefix>_genes.final.txt`

The final filtered gene list.

### Stage 5: Summary (Always)

**Output**: `<prefix>_filtering_summary.txt`

Contains:
- Gene counts at each stage
- Number removed by each filter
- Classification statistics

## Example Workflows

### Example 1: GO Terms Only (v1.0 behavior)

```bash
conda run -n wot_env python create_list_from_GAF.py config_PCD_ROS.yaml
```

**Output**: `PCD_ROS_genes.txt` (510 genes)

### Example 2: GO + Orthology Filter

```bash
conda run -n wot_env python build_pcd_list.py --config config_PCD_ROS_with_ortho.yaml
```

**Outputs**:
- `PCD_ROS_ortho_genes.base.txt` (510 genes)
- `PCD_ROS_ortho_genes.orthology_filtered.txt` (211 genes, 1:1 orthologues only)
- `PCD_ROS_ortho_genes.final.txt` (211 genes)
- `PCD_ROS_ortho_orthology_filter_report.tsv`
- `PCD_ROS_ortho_filtering_summary.txt`

### Example 3: GO + Domain + Orthology

```yaml
# config_all_filters.yaml
filters:
  domain_filter:
    enabled: true
    allowed_domains:
      pfam: ["PF00931", "PF12613"]  # NB-ARC, LRR
  
  orthology_filter:
    enabled: true
    prefer_one_to_one: true
```

```bash
conda run -n wot_env python build_pcd_list.py --config config_all_filters.yaml
```

**Processing order**: GO → Domain filter → Orthology filter → Final

### Example 4: Using Arabidopsis Reference Set

Create an Arabidopsis PCD gene list:
```bash
# arabidopsis_pcd_genes.txt
AT1G17210
AT2G34690
AT3G25070
...
```

Configure orthology filter:
```yaml
filters:
  orthology_filter:
    enabled: true
    require_arabidopsis_in_set: true
    arabidopsis_pcd_gene_list: "arabidopsis_pcd_genes.txt"
    prefer_one_to_one: false  # Allow 1:many if in reference set
```

This keeps only Brachypodium genes whose Arabidopsis orthologue is in the reference list.

## File Formats

### Input: GO Annotation (TSV)

```
gene    GO          level
BdiBd21-3.3G0362100.1.p   GO:0008219   BP
BdiBd21-3.3G0362100.1.p   GO:0003824   MF
```

### Input: Orthology (inParanoid TSV)

```
OrtoA                                          OrtoB
BdistachyonBd21-3:BdiBd21-3.1G0256200.1 1.000  Athalianacolumbia:AT2G17930.1 1.000 Athalianacolumbia:AT4G36080.1 0.561
BdistachyonBd21-3:BdiBd21-3.4G0117900.2 1.000  Athalianacolumbia:AT3G02260.1 1.000
```

### Input: Domain Annotation (TSV, flexible columns)

```
gene                    pfam      interpro
BdiBd21-3.1G0001000.1.p PF00931   IPR002182
BdiBd21-3.1G0001000.1.p PF12613   IPR001611
```

### Output: Gene Lists (TXT)

```
BdiBd21-3.1G0018700.v1.2
BdiBd21-3.1G0022000.v1.2
BdiBd21-3.1G0025900.v1.2
...
```

One gene per line, Seurat-compatible format.

## Troubleshooting

### No genes after orthology filter

**Problem**: All genes removed by orthology filter

**Solutions**:
1. Check `prefer_one_to_one` setting - try `false` to allow 1:many
2. Lower `min_best_score` threshold (e.g., 0.7 instead of 0.8)
3. Check orthology file path is correct
4. Verify gene ID regex patterns match your data

### Domain filter removes all genes

**Problem**: No genes have specified domains

**Solutions**:
1. Verify domain IDs are correct (check domain annotation file)
2. Check `id_column` and `domain_columns` match your file
3. Try adding more domain IDs to `allowed_domains`

### Gene ID mapping issues

**Problem**: Genes not found in orthology/domain files

**Solutions**:
1. Check `protein_to_core_regex` pattern
2. Verify `seurat_suffix` is correct
3. Use `grep` to check ID format in files

## Version Comparison

### v1.0 (create_list_from_GAF.py)
- GO term-based extraction
- Category annotation (PCD, ROS, etc.)
- Simple, fast

### v2.0 (build_pcd_list.py)
- All v1.0 features
- Optional domain filtering
- Optional orthology filtering
- Multi-stage outputs
- Detailed reports
- More complex configuration

## Citation

If using in publications, document:
1. GO terms used (include GO IDs)
2. Whether filters were applied
3. Filter parameters (domain IDs, orthology rules)
4. Final gene count

See `PUBLICATION_CHECKLIST.md` for complete guidelines.

## Support

For issues:
1. Check config file syntax (YAML format)
2. Verify all file paths exist
3. Check gene ID formats match expected patterns
4. Review filtering summary for diagnostic info

GitHub: https://github.com/El-Castor/PlantGeneSignatureBuilder
