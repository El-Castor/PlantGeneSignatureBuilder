#!/usr/bin/env python3
"""
Gene List Extraction from GO Annotations
Reads GO terms from a YAML config file and extracts matching genes.

Usage:
    python create_list_from_GAF.py [config_file.yaml]
    
If no config file is specified, uses 'config_PCD_stress.yaml' by default.
"""

import pandas as pd
import re
import sys
from pathlib import Path
try:
    import yaml
except ImportError:
    print("ERROR: PyYAML is required. Install with: pip install pyyaml")
    sys.exit(1)


def load_config(config_path):
    """Load configuration from YAML file."""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def extract_genes(config):
    """Extract genes based on GO terms specified in config."""
    
    # ====== LOAD CONFIG ======
    GO_TSV = Path(config['input_go_file'])
    OUT = Path(config['output_file'])
    VERSION_SUFFIX = config['gene_id_suffix']
    GENE_PATTERN = config['gene_id_pattern']
    GO_LEVEL = config.get('go_level_filter', 'BP')
    
    # Extract GO term IDs from config
    GO_ROOTS = {term['go_id'] for term in config['go_terms']}
    
    print(f"Configuration:")
    print(f"  Input file: {GO_TSV}")
    print(f"  Output file: {OUT}")
    print(f"  Gene ID suffix: {VERSION_SUFFIX}")
    print(f"  GO level filter: {GO_LEVEL}")
    print(f"  Number of GO terms: {len(GO_ROOTS)}")
    print()

    # ====== LOAD GO ANNOTATIONS ======
    print(f"Loading GO annotations from {GO_TSV}...")
    go = pd.read_csv(GO_TSV, sep="\t")
    print(f"  Loaded {len(go):,} annotations")
    
    # ====== MAP protein-id -> gene-id ======
    # Example: BdiBd21-3.3G0362100.1.p -> BdiBd21-3.3G0362100
    print(f"Mapping protein IDs to gene IDs...")
    core_re = re.compile(GENE_PATTERN)
    go["core"] = go["gene"].astype(str).str.extract(core_re, expand=False)
    
    # Drop unmappable lines
    n_before = len(go)
    go = go.dropna(subset=["core"])
    n_dropped = n_before - len(go)
    if n_dropped > 0:
        print(f"  Warning: {n_dropped} annotations could not be mapped")
    
    # Create Seurat-compatible gene IDs
    go["seurat_gene"] = go["core"] + VERSION_SUFFIX
    print(f"  Mapped {go['seurat_gene'].nunique():,} unique genes")
    
    # ====== FILTER GO TERMS ======
    print(f"Filtering for target GO terms...")
    pcd = go[go["GO"].isin(GO_ROOTS)].copy()
    print(f"  Found {len(pcd):,} annotations matching target GO terms")
    
    # Filter by GO level if specified
    if GO_LEVEL and GO_LEVEL.upper() != "ALL":
        n_before = len(pcd)
        pcd = pcd[pcd["level"] == GO_LEVEL.upper()]
        print(f"  Filtered to {GO_LEVEL} only: {len(pcd):,} annotations ({n_before - len(pcd)} removed)")
    
    # ====== EXTRACT UNIQUE GENES ======
    genes = sorted(set(pcd["seurat_gene"].astype(str)))
    
    # ====== CREATE GENE ANNOTATION TABLE ======
    # Build mapping dictionaries from config
    go_term_dict = {term['go_id']: term['description'] for term in config['go_terms']}
    go_category_dict = {term['go_id']: term.get('category', 'Unclassified') for term in config['go_terms']}
    
    # Create a dataframe with gene -> GO term -> annotation mapping
    gene_annotations = []
    for gene in genes:
        # Get all GO terms for this gene
        gene_go_terms = pcd[pcd['seurat_gene'] == gene]['GO'].unique()
        
        # Collect categories and descriptions
        categories = set()
        descriptions = set()
        go_ids = []
        
        for go_id in gene_go_terms:
            go_ids.append(go_id)
            if go_id in go_category_dict:
                categories.add(go_category_dict[go_id])
            if go_id in go_term_dict:
                descriptions.add(go_term_dict[go_id])
        
        # Create category label
        if categories:
            category_label = ", ".join(sorted(categories))
        else:
            category_label = "Unclassified"
        
        # Create description label
        if descriptions:
            description_label = "; ".join(sorted(descriptions))
        else:
            description_label = "Unknown"
        
        gene_annotations.append({
            'gene': gene,
            'category': category_label,
            'description': description_label
        })
    
    # Create DataFrame
    gene_df = pd.DataFrame(gene_annotations)
    
    # ====== WRITE OUTPUT ======
    # Write as TSV with header
    gene_df.to_csv(OUT, sep="\t", index=False)
    print(f"\n✓ SUCCESS: {len(genes)} genes written to {OUT}")
    print(f"  Columns: gene, category, description")
    
    # Print summary by GO term
    if len(genes) > 0:
        print(f"\nGene count by GO term:")
        go_summary = pcd.groupby('GO')['seurat_gene'].nunique().sort_values(ascending=False)
        for go_id, count in go_summary.items():
            desc = go_term_dict.get(go_id, "Unknown")
            print(f"  {go_id}: {count:4d} genes - {desc}")
        
        print(f"\nGene count by category:")
        category_summary = gene_df.explode('category')['category'].value_counts()
        for cat, count in category_summary.items():
            print(f"  {cat}: {count} genes")
    
    return genes


def main():
    """Main entry point."""
    # Get config file from command line or use default
    if len(sys.argv) > 1:
        config_file = sys.argv[1]
    else:
        config_file = "config_PCD_stress.yaml"
    
    config_path = Path(config_file)
    
    if not config_path.exists():
        print(f"ERROR: Config file not found: {config_path}")
        print(f"\nUsage: python {sys.argv[0]} [config_file.yaml]")
        sys.exit(1)
    
    print("=" * 80)
    print("GO Gene List Extraction Tool")
    print("=" * 80)
    print(f"Config file: {config_path}\n")
    
    # Load config and extract genes
    config = load_config(config_path)
    genes = extract_genes(config)
    
    print("=" * 80)


if __name__ == "__main__":
    main()
