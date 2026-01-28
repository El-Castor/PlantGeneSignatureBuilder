# Custom Gene Signature Configurations

This directory contains example configuration files for different biological gene signatures.

These configs are designed for the **v1.0 simple extraction system** (`create_list_from_GAF.py`).

For the **v3.0 scoring system** (`rank_gene_signatures.py`), use the templates in the main directory.

---

## Available Signatures

### 1. PCD/ROS Genes
**File:** `config_PCD_ROS.yaml`

**GO Terms:**
- Programmed cell death (GO:0012501)
- Cell death (GO:0008219)
- Regulation of programmed cell death (GO:0043067)
- Regulation of cell death (GO:0010941)
- Response to oxidative stress (GO:0006979)
- ROS metabolic process (GO:0072593)
- Cell redox homeostasis (GO:0045454)

**Use:** PCD and oxidative stress response studies

---

### 2. PCD + General Stress
**File:** `config_PCD_stress.yaml`

Includes PCD/ROS terms plus:
- Response to heat (GO:0009408)
- Response to cold (GO:0009409)
- Response to hydrogen peroxide (GO:0042542)

**Use:** Broader stress response analysis

---

### 3. Cell Division
**File:** `config_cell_division.yaml`

**GO Terms:**
- Cell division (GO:0051301)
- Mitotic cell cycle (GO:0000278)
- Regulation of cell cycle (GO:0051726)

**Use:** Cell cycle and proliferation studies

---

### 4. Photosynthesis
**File:** `config_photosynthesis.yaml`

**GO Terms:**
- Photosynthesis (GO:0015979)
- Photosynthesis, light reaction (GO:0019684)
- Photosynthesis, dark reaction (GO:0019685)

**Use:** Photosynthetic pathway analysis

---

## Usage

### With v1.0 Simple Extraction
```bash
python create_list_from_GAF.py --config config_custom/config_PCD_ROS.yaml
```

### Migrating to v3.0 Scoring System

To use with the scoring system, convert format:

**OLD (v1.0):**
```yaml
go_terms:
  - go_id: "GO:0012501"
    description: "programmed cell death"
    category: "PCD"
```

**NEW (v3.0):**
```yaml
score_weights:
  go_core_terms:
    GO:0012501: 3  # programmed cell death
```

See `config_scoring_template.yaml` in main directory for v3.0 format.

---

## Creating Custom Signatures

1. Copy an existing config:
```bash
cp config_custom/config_PCD_ROS.yaml config_custom/config_my_signature.yaml
```

2. Edit GO terms for your pathway of interest

3. Run extraction:
```bash
python create_list_from_GAF.py --config config_custom/config_my_signature.yaml
```

---

## File Format (v1.0)

```yaml
output_prefix: "my_genes"

input_go_file: "/path/to/GOannotation.tsv"

id_mapping:
  protein_to_core_regex: "^(BdiBd21-3\\.\\dG\\d{7})"
  seurat_suffix: ".v1.2"

go_level_filter: "BP"

go_terms:
  - go_id: "GO:xxxxxxx"
    description: "term description"
    category: "category_name"
```

---

## Converting to v3.0

For multi-evidence scoring with domains, orthology, and Arabidopsis context:

```bash
# Use v3.0 templates
cp config_scoring_template.yaml config_my_scoring.yaml

# Add your GO terms to score_weights section
# Enable evidence layers as needed
python rank_gene_signatures.py --config config_my_scoring.yaml
```
