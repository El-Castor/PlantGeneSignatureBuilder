# GO Gene List Extraction Tool

A flexible Python script for extracting gene lists from Gene Ontology (GO) annotations based on user-defined GO terms.

## Features

- **Configurable GO terms**: Define target GO terms in a separate YAML config file
- **Flexible ID mapping**: Customizable regex patterns for protein-to-gene ID conversion
- **GO level filtering**: Filter by Biological Process (BP), Molecular Function (MF), Cellular Component (CC), or all
- **Detailed reporting**: Shows gene counts per GO term and extraction statistics
- **Reusable**: Create multiple config files for different gene sets

## Requirements

```bash
conda install -n wot_env pyyaml pandas
```

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
# Install PyYAML
conda install -n wot_env pyyaml
# or
pip install pyyaml
```

## License

Free to use and modify for research purposes.
