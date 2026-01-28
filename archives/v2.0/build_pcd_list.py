#!/usr/bin/env python3
"""
Plant Gene Signature Builder v2.0
Multi-stage gene list construction with optional domain and orthology filtering

This script extends the GO-based gene extraction with:
1. Base GO term extraction (from create_list_from_GAF.py logic)
2. Optional domain-based filtering (Pfam/InterPro)
3. Optional orthology-based filtering (Arabidopsis orthologues)

Usage:
    python build_pcd_list.py --config config_extended.yaml
"""

import pandas as pd
import re
import sys
import argparse
from pathlib import Path
from collections import Counter, defaultdict
try:
    import yaml
except ImportError:
    print("ERROR: PyYAML is required. Install with: pip install pyyaml")
    sys.exit(1)


class GeneListBuilder:
    """Main class for building filtered gene lists"""
    
    def __init__(self, config_path):
        self.config_path = Path(config_path)
        self.config = self._load_config()
        self.output_prefix = self.config.get('output_prefix', 'gene_list')
        self.stats = {
            'base': 0,
            'domain_filtered': 0,
            'orthology_filtered': 0,
            'final': 0
        }
        self.reports = {}
        
    def _load_config(self):
        """Load and validate configuration"""
        with open(self.config_path, 'r') as f:
            config = yaml.safe_load(f)
        return config
    
    def run(self):
        """Execute the full pipeline"""
        print("="*80)
        print("Plant Gene Signature Builder v2.0")
        print("="*80)
        print(f"Config: {self.config_path}\n")
        
        # Step 1: Base GO extraction
        print("STEP 1: Base GO term extraction")
        print("-" * 80)
        base_genes = self._extract_base_genes()
        self.stats['base'] = len(base_genes)
        self._write_output(base_genes, f"{self.output_prefix}_genes.base.txt")
        print(f"✓ Base genes: {len(base_genes)}\n")
        
        current_genes = base_genes.copy()
        
        # Step 2: Domain filtering (optional)
        if self.config.get('filters', {}).get('domain_filter', {}).get('enabled', False):
            print("STEP 2: Domain-based filtering")
            print("-" * 80)
            current_genes = self._apply_domain_filter(current_genes)
            self.stats['domain_filtered'] = len(current_genes)
            self._write_output(current_genes, f"{self.output_prefix}_genes.domain_filtered.txt")
            print(f"✓ After domain filter: {len(current_genes)}\n")
        else:
            print("STEP 2: Domain filtering [DISABLED]\n")
        
        # Step 3: Orthology filtering (optional)
        if self.config.get('filters', {}).get('orthology_filter', {}).get('enabled', False):
            print("STEP 3: Orthology-based filtering")
            print("-" * 80)
            current_genes = self._apply_orthology_filter(current_genes)
            self.stats['orthology_filtered'] = len(current_genes)
            self._write_output(current_genes, f"{self.output_prefix}_genes.orthology_filtered.txt")
            print(f"✓ After orthology filter: {len(current_genes)}\n")
        else:
            print("STEP 3: Orthology filtering [DISABLED]\n")
        
        # Step 4: Final output
        print("STEP 4: Final output")
        print("-" * 80)
        self.stats['final'] = len(current_genes)
        self._write_output(current_genes, f"{self.output_prefix}_genes.final.txt")
        print(f"✓ Final gene list: {len(current_genes)}\n")
        
        # Step 5: Summary
        self._write_summary()
        print("="*80)
        print("Pipeline completed successfully!")
        print("="*80)
    
    def _extract_base_genes(self):
        """Extract genes based on GO terms (Step 1)"""
        # Load GO annotation
        go_file = Path(self.config['input_go_file'])
        print(f"Loading GO annotations: {go_file}")
        go_df = pd.read_csv(go_file, sep="\t")
        print(f"  Loaded {len(go_df):,} annotations")
        
        # Get GO terms from config
        go_terms = {term['go_id'] for term in self.config['go_terms']}
        print(f"  Target GO terms: {len(go_terms)}")
        
        # ID mapping
        id_mapping = self.config.get('id_mapping', {})
        pattern = id_mapping.get('protein_to_core_regex', r"^(BdiBd21-3\.\dG\d{7})")
        suffix = id_mapping.get('seurat_suffix', '.v1.2')
        
        print(f"  Mapping protein IDs to gene IDs...")
        core_re = re.compile(pattern)
        go_df["core"] = go_df["gene"].astype(str).str.extract(core_re, expand=False)
        go_df = go_df.dropna(subset=["core"])
        go_df["seurat_gene"] = go_df["core"] + suffix
        
        # Filter for GO terms
        go_level = self.config.get('go_level_filter', 'BP')
        filtered = go_df[go_df["GO"].isin(go_terms)].copy()
        
        if go_level and go_level.upper() != "ALL":
            filtered = filtered[filtered["level"] == go_level.upper()]
        
        genes = sorted(set(filtered["seurat_gene"].astype(str)))
        
        # Store annotation info for reports
        self.go_annotations = filtered
        
        return genes
    
    def _apply_domain_filter(self, genes):
        """Filter genes based on domain annotations (Step 2)"""
        domain_config = self.config['filters']['domain_filter']
        
        # Load domain annotations
        domain_file = Path(domain_config['annotation_file'])
        print(f"Loading domain annotations: {domain_file}")
        
        try:
            domain_df = pd.read_csv(domain_file, sep="\t")
        except Exception as e:
            print(f"  ERROR: Could not read domain file: {e}")
            return genes
        
        print(f"  Loaded {len(domain_df):,} domain annotations")
        
        # Get column names
        id_col = domain_config.get('id_column', 'gene')
        domain_cols = domain_config.get('domain_columns', ['pfam', 'interpro'])
        
        # Get allowed domains
        allowed_domains = domain_config.get('allowed_domains', {})
        
        # Map gene IDs if needed
        suffix = self.config['id_mapping']['seurat_suffix']
        pattern = self.config['id_mapping']['protein_to_core_regex']
        core_re = re.compile(pattern)
        
        domain_df["core"] = domain_df[id_col].astype(str).str.extract(core_re, expand=False)
        domain_df = domain_df.dropna(subset=["core"])
        domain_df["seurat_gene"] = domain_df["core"] + suffix
        
        # Filter for genes in our list
        domain_df = domain_df[domain_df['seurat_gene'].isin(genes)]
        
        # Find matches
        matched_genes = set()
        match_report = []
        
        for domain_type, domain_ids in allowed_domains.items():
            if domain_type not in domain_cols:
                continue
            
            for _, row in domain_df.iterrows():
                gene = row['seurat_gene']
                domain_val = str(row.get(domain_type, ''))
                
                if pd.isna(domain_val) or domain_val == 'nan':
                    continue
                
                # Check if any allowed domain is present
                for allowed_id in domain_ids:
                    if allowed_id in domain_val:
                        matched_genes.add(gene)
                        match_report.append({
                            'gene_id': gene,
                            'matched_domain_type': domain_type,
                            'matched_domain_id': allowed_id
                        })
        
        # Save report
        if match_report:
            report_df = pd.DataFrame(match_report)
            report_file = f"{self.output_prefix}_domain_filter_report.tsv"
            report_df.to_csv(report_file, sep="\t", index=False)
            print(f"  Domain filter report: {report_file}")
            self.reports['domain'] = report_df
        
        filtered_genes = sorted(matched_genes)
        print(f"  Genes with matching domains: {len(filtered_genes)} / {len(genes)}")
        print(f"  Removed: {len(genes) - len(filtered_genes)}")
        
        return filtered_genes
    
    def _apply_orthology_filter(self, genes):
        """Filter genes based on orthology to Arabidopsis (Step 3)"""
        ortho_config = self.config['filters']['orthology_filter']
        
        # Load orthology file
        ortho_file = Path(ortho_config['orthology_file'])
        print(f"Loading orthology data: {ortho_file}")
        
        try:
            ortho_df = pd.read_csv(ortho_file, sep="\t")
        except Exception as e:
            print(f"  ERROR: Could not read orthology file: {e}")
            return genes
        
        print(f"  Loaded {len(ortho_df):,} orthology pairs")
        
        # Parse orthology data
        ortho_data = self._parse_orthology_file(ortho_df, ortho_config)
        
        # Load Arabidopsis gene list if provided
        arabidopsis_genes = None
        if ortho_config.get('arabidopsis_pcd_gene_list'):
            arab_file = Path(ortho_config['arabidopsis_pcd_gene_list'])
            if arab_file.exists():
                with open(arab_file, 'r') as f:
                    arabidopsis_genes = set(line.strip() for line in f if line.strip())
                print(f"  Loaded {len(arabidopsis_genes)} Arabidopsis reference genes")
        
        # Filter genes
        filtered_genes = []
        orthology_report = []
        
        require_arab_in_set = ortho_config.get('require_arabidopsis_in_set', False)
        prefer_one_to_one = ortho_config.get('prefer_one_to_one', False)
        one_to_one_rule = ortho_config.get('one_to_one_rule', {})
        
        for gene in genes:
            if gene in ortho_data:
                ortho_info = ortho_data[gene]
                arab_hits = ortho_info['arabidopsis_hits']
                classification = ortho_info['classification']
                
                # Apply filters
                keep = True
                
                # Check if Arabidopsis orthologue is in reference set
                if require_arab_in_set and arabidopsis_genes:
                    if not any(hit in arabidopsis_genes for hit in arab_hits):
                        keep = False
                
                # Check one-to-one preference
                if prefer_one_to_one and classification != 'one_to_one':
                    keep = False
                
                if keep:
                    filtered_genes.append(gene)
                    orthology_report.append({
                        'brachy_gene_id': gene,
                        'brachy_protein_id': ortho_info['protein_id'],
                        'arabidopsis_hits': ', '.join(arab_hits),
                        'best_hit_score': ortho_info['best_score'],
                        'n_hits': len(arab_hits),
                        'classification': classification,
                        'in_reference_set': any(hit in arabidopsis_genes for hit in arab_hits) if arabidopsis_genes else 'N/A'
                    })
        
        # Save report
        if orthology_report:
            report_df = pd.DataFrame(orthology_report)
            report_file = f"{self.output_prefix}_orthology_filter_report.tsv"
            report_df.to_csv(report_file, sep="\t", index=False)
            print(f"  Orthology filter report: {report_file}")
            self.reports['orthology'] = report_df
            
            # Print classification summary
            print(f"\n  Orthology classification:")
            for cls, count in report_df['classification'].value_counts().items():
                print(f"    {cls}: {count}")
        
        filtered_genes = sorted(filtered_genes)
        print(f"\n  Genes passing orthology filter: {len(filtered_genes)} / {len(genes)}")
        print(f"  Removed: {len(genes) - len(filtered_genes)}")
        
        return filtered_genes
    
    def _parse_orthology_file(self, ortho_df, ortho_config):
        """Parse inParanoid orthology file format"""
        mapping_rules = ortho_config.get('mapping_rules', {})
        brachy_regex = mapping_rules.get('brachy_id_regex_extract', r"BdiBd21-3\.\dG\d{7}")
        arab_regex = mapping_rules.get('arabidopsis_id_regex_extract', r"AT[1-5MC]G\d{5}")
        one_to_one_rule = ortho_config.get('one_to_one_rule', {})
        
        brachy_pattern = re.compile(brachy_regex)
        arab_pattern = re.compile(arab_regex)
        suffix = self.config['id_mapping']['seurat_suffix']
        
        ortho_data = {}
        
        for _, row in ortho_df.iterrows():
            # Parse OrtoA (Brachypodium)
            orto_a = str(row['OrtoA'])
            brachy_match = brachy_pattern.search(orto_a)
            if not brachy_match:
                continue
            
            brachy_core = brachy_match.group()
            brachy_gene = brachy_core + suffix
            
            # Extract protein ID
            protein_id = orto_a.split()[0].split(':')[-1] if ':' in orto_a else brachy_core
            
            # Parse OrtoB (Arabidopsis)
            orto_b = str(row['OrtoB'])
            arab_hits = []
            arab_scores = []
            
            # Split by whitespace and extract Arabidopsis IDs and scores
            tokens = orto_b.split()
            for i, token in enumerate(tokens):
                arab_match = arab_pattern.search(token)
                if arab_match:
                    arab_id = arab_match.group()
                    # Next token should be score
                    score = float(tokens[i+1]) if i+1 < len(tokens) else 0.0
                    arab_hits.append(arab_id)
                    arab_scores.append(score)
            
            if not arab_hits:
                continue
            
            # Classify relationship
            best_score = max(arab_scores)
            n_hits = len(arab_hits)
            
            classification = self._classify_orthology(
                n_hits, arab_scores, one_to_one_rule
            )
            
            ortho_data[brachy_gene] = {
                'protein_id': protein_id,
                'arabidopsis_hits': arab_hits,
                'scores': arab_scores,
                'best_score': best_score,
                'n_hits': n_hits,
                'classification': classification
            }
        
        return ortho_data
    
    def _classify_orthology(self, n_hits, scores, rule):
        """Classify orthology relationship as one-to-one, one-to-many, or low confidence"""
        if not scores:
            return 'low_confidence'
        
        max_hits = rule.get('max_hits', 1)
        min_best_score = rule.get('min_best_score', 0.8)
        min_second_best_gap = rule.get('min_second_best_gap', 0.2)
        
        best_score = max(scores)
        
        # Low confidence if best score too low
        if best_score < min_best_score:
            return 'low_confidence'
        
        # One-to-one if single hit and good score
        if n_hits == 1:
            return 'one_to_one'
        
        # Check gap to second best
        if n_hits > 1:
            sorted_scores = sorted(scores, reverse=True)
            gap = sorted_scores[0] - sorted_scores[1]
            
            if gap >= min_second_best_gap and n_hits <= max_hits:
                return 'one_to_one'
        
        return 'one_to_many'
    
    def _write_output(self, genes, filename):
        """Write gene list to file"""
        output_file = Path(filename)
        output_file.write_text("\n".join(genes) + "\n")
    
    def _write_summary(self):
        """Write filtering summary"""
        summary_file = f"{self.output_prefix}_filtering_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("Plant Gene Signature Builder v2.0 - Filtering Summary\n")
            f.write("="*80 + "\n\n")
            
            f.write(f"Config file: {self.config_path}\n")
            f.write(f"Output prefix: {self.output_prefix}\n\n")
            
            f.write("GENE COUNTS BY STAGE:\n")
            f.write("-" * 80 + "\n")
            f.write(f"  Base (GO terms):        {self.stats['base']:5d}\n")
            
            if self.config.get('filters', {}).get('domain_filter', {}).get('enabled', False):
                f.write(f"  After domain filter:    {self.stats['domain_filtered']:5d}")
                f.write(f"  (removed: {self.stats['base'] - self.stats['domain_filtered']})\n")
            
            if self.config.get('filters', {}).get('orthology_filter', {}).get('enabled', False):
                prev_count = self.stats.get('domain_filtered', self.stats['base'])
                f.write(f"  After orthology filter: {self.stats['orthology_filtered']:5d}")
                f.write(f"  (removed: {prev_count - self.stats['orthology_filtered']})\n")
            
            f.write(f"  Final:                  {self.stats['final']:5d}\n\n")
            
            # Domain filter details
            if 'domain' in self.reports:
                f.write("DOMAIN FILTER DETAILS:\n")
                f.write("-" * 80 + "\n")
                domain_counts = self.reports['domain']['matched_domain_type'].value_counts()
                for dtype, count in domain_counts.items():
                    f.write(f"  {dtype}: {count} matches\n")
                f.write("\n")
            
            # Orthology filter details
            if 'orthology' in self.reports:
                f.write("ORTHOLOGY FILTER DETAILS:\n")
                f.write("-" * 80 + "\n")
                ortho_counts = self.reports['orthology']['classification'].value_counts()
                for cls, count in ortho_counts.items():
                    f.write(f"  {cls}: {count} genes\n")
                f.write("\n")
        
        print(f"✓ Summary saved: {summary_file}")
        
        # Print to stdout
        with open(summary_file, 'r') as f:
            print("\n" + f.read())


def main():
    parser = argparse.ArgumentParser(
        description='Plant Gene Signature Builder v2.0 - Multi-stage gene list construction'
    )
    parser.add_argument(
        '--config', '-c',
        required=True,
        help='Path to YAML configuration file'
    )
    
    args = parser.parse_args()
    
    # Run pipeline
    builder = GeneListBuilder(args.config)
    builder.run()


if __name__ == "__main__":
    main()
