# Custom Gene Signature Configurations

This directory contains custom configuration files for different biological processes.

Each config file is tailored for the **v3.1 scoring system** (`rank_gene_signatures.py`) and includes:
- Category-specific GO terms and weights
- Expected Pfam domains for each category
- Synergy rules for multi-evidence validation
- Domain scanning integration with hmmscan

---

## Available Configurations

### Development
**File:** `config_embryogenesis.yaml`

**Categories:** Embryo_Development, Zygote_Formation, Cell_Fate

**GO Terms:**
- Embryo development (GO:0009790)
- Zygote formation (GO:0000381)
- Cell fate determination (GO:0001709)
- Pattern specification (GO:0007389)
- And 13 more developmental GO terms

**Use:** Embryogenesis and early development studies

---

### Stress Responses
**File:** `config_abiotic_stress.yaml`

**Categories:** Drought, Heat, Cold, Salt, Oxidative

**GO Terms:**
- Response to water deprivation (GO:0009414)
- Response to heat (GO:0009408)
- Response to cold (GO:0009409)
- Response to salt stress (GO:0009651)
- Response to oxidative stress (GO:0006979)
- And more stress-specific terms

**Expected Domains:**
- HSP70, HSP90 (heat stress)
- LEA proteins (cold/drought)
- Ion transporters (salt stress)
- Peroxidases, catalases (oxidative stress)

**Use:** Abiotic stress tolerance studies

---

### Cell Biology
**File:** `config_cell_cycle.yaml`

**Categories:** Mitosis, Meiosis, Cytokinesis, Cell_Cycle_Regulation

**GO Terms:**
- Mitotic cell cycle (GO:0000278)
- Meiotic cell cycle (GO:0051321)
- Cytokinesis (GO:0000910)
- Regulation of cell cycle (GO:0051726)
- G1/S and G2/M transitions

**Expected Domains:**
- Cyclins (PF00149)
- Protein kinases (CDKs)
- Kinesin motors
- F-box proteins (cell cycle regulators)

**Use:** Cell division and proliferation studies

---

### Metabolism
**File:** `config_photosynthesis.yaml`

**Categories:** Light_Reactions, Carbon_Fixation, Chloroplast, Pigment_Biosynthesis

**GO Terms:**
- Photosynthesis light reactions (GO:0019684)
- Carbon fixation (GO:0015977)
- Chloroplast organization (GO:0009658)
- Chlorophyll biosynthesis (GO:0015995)

**Expected Domains:**
- Photosystem I/II proteins
- RuBisCO
- Cytochrome b5
- Mg-chelatase (chlorophyll biosynthesis)

**Use:** Photosynthesis and chloroplast function studies

---

## Usage

```bash
# Embryogenesis
conda run -n wot_env python rank_gene_signatures.py \
  --config config_custom/config_embryogenesis.yaml \
  --run_name embryogenesis

# Abiotic stress
conda run -n wot_env python rank_gene_signatures.py \
  --config config_custom/config_abiotic_stress.yaml \
  --run_name abiotic_stress

# Cell cycle
conda run -n wot_env python rank_gene_signatures.py \
  --config config_custom/config_cell_cycle.yaml \
  --run_name cell_cycle

# Photosynthesis
conda run -n wot_env python rank_gene_signatures.py \
  --config config_custom/config_photosynthesis.yaml \
  --run_name photosynthesis
```

## Output Files

Each run generates:
- `{prefix}_scored_genes.base.txt` - All genes with target GO terms
- `{prefix}_scored_genes.scored.tsv` - Full scoring breakdown
- `{prefix}_scored_genes.ALL_overview.tsv` - Complete overview
- `{prefix}_scored_genes.HIGH_overview.tsv` - High-confidence genes (top 10%)
- `{prefix}_scored_genes.HIGH_category.tsv` - Gene IDs with category labels

The HIGH_category.tsv file shows which category(ies) each gene belongs to based on GO terms.
