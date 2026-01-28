# Arabidopsis Resources

This directory contains Arabidopsis reference files required for orthology-based symbol enrichment.

## Required Files

### GFF3 Annotation File
**File:** `Araport11_GFF3_genes_transposons.20250813.gff`  
**Size:** ~275 MB (too large for GitHub)  
**Source:** [TAIR/Araport11](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FAraport11_genome_release)

**Download:**
```bash
cd resources/arabidopsis/
wget https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz
gunzip Araport11_GFF3_genes_transposons.201606.gff.gz
```

Or use your own version of the Araport11 GFF3 file and update the path in your config file.

### Plant Ontology Annotation File
**File:** `Arabidopsis_annotation_gene.txt`  
**Source:** User-provided or custom annotation mapping

This file maps Arabidopsis gene IDs to Plant Ontology (PO) terms and descriptions.

## Configuration

Update your config YAML file with the correct paths:

```yaml
arabidopsis_context:
  gff3_file: "resources/arabidopsis/Araport11_GFF3_genes_transposons.20250813.gff"
  po_annotation_file: "resources/arabidopsis/Arabidopsis_annotation_gene.txt"
  enabled: false  # Set to true only if using keyword scoring
```

**Note:** The GFF3 file is loaded automatically when orthology scoring is enabled, even if `arabidopsis_context.enabled=false`.
