#!/usr/bin/env python3
"""
Plant Gene Signature Ranker - Multi-Evidence Scoring System v3.1

SCIENTIFICALLY DEFENSIBLE SCORING with caps + synergy bonuses + run management

Features:
- Capped scoring per evidence layer (prevents GO dominance)
- Synergy bonuses when independent evidence agrees
- QC plots and sensitivity analysis
- Run directory management for reproducible analysis
- Enriched overview outputs with aggregated annotations
- Fully OFFLINE (no web queries)

Usage:
    python rank_gene_signatures.py --config config.yaml [--qc] [--run_name NAME] [--overwrite]
"""

import pandas as pd
import numpy as np
import re
import sys
import argparse
import hashlib
import json
from pathlib import Path
from collections import Counter, defaultdict
try:
    import yaml
except ImportError:
    print("ERROR: PyYAML required. Install: pip install pyyaml")
    sys.exit(1)

# Run management
try:
    from pgsb.io.run_manager import (
        create_run_dir, write_config_snapshot, 
        update_latest_symlink, create_manifest
    )
except ImportError:
    print("ERROR: pgsb package not found. Make sure pgsb/ is in PYTHONPATH or current directory.")
    sys.exit(1)

# Optional plotting
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False


class GeneSignatureRanker:
    """Multi-evidence scoring with caps and synergy bonuses - OFFLINE MODE"""
    
    def __init__(self, config_path, run_name=None, overwrite=False):
        self.config_path = Path(config_path)
        self.config = self._load_config()
        self.output_prefix = self.config.get('output_prefix', 'gene_signatures')
        
        # Initialize run directory
        self.run_info = create_run_dir(
            self.config, 
            self.output_prefix,
            run_name=run_name,
            overwrite=overwrite
        )
        
        # Gene data structures
        self.gene_scores_raw = {}
        self.gene_scores_capped = {}
        self.gene_evidence = {}
        self.synergy_bonuses = {}
        
        # Annotation caches
        self.arabidopsis_annotations = {}
        self.arabidopsis_po_context = {}
        
    def _load_config(self):
        with open(self.config_path, 'r') as f:
            return yaml.safe_load(f)
    
    def run(self, enable_qc=False):
        """Execute scoring pipeline"""
        print("="*80)
        print("Plant Gene Signature Ranker v3.1 - Capped Multi-Evidence Scoring")
        print("="*80)
        print(f"Config: {self.config_path}")
        print(f"Run ID: {self.run_info['run_id']}")
        print(f"Output dir: {self.run_info['outputs']}")
        print(f"Mode: OFFLINE\n")
        
        # Save config snapshot
        write_config_snapshot(self.config, self.run_info['run_dir'])
        
        # Load Arabidopsis annotations for enrichment (if orthology or keyword scoring enabled)
        arab_config = self.config.get('evidence', {}).get('arabidopsis_context', {})
        ortho_enabled = self.config.get('evidence', {}).get('orthology_evidence', {}).get('enabled', False)
        
        if ortho_enabled or arab_config.get('enabled', False):
            print("STEP 0: Loading Arabidopsis annotations (OFFLINE)")
            print("-" * 80)
            self._load_arabidopsis_annotations()
            print()
        
        # Extract base genes
        print("STEP 1: Base GO term extraction")
        print("-" * 80)
        base_genes = self._extract_base_genes()
        print(f"✓ Base genes: {len(base_genes)}\n")
        
        # Initialize scoring
        for gene in base_genes:
            self.gene_scores_raw[gene] = {
                'go': 0, 'domain': 0, 'orthology': 0,
                'tair': 0, 'po': 0, 'synergy': 0
            }
            self.gene_scores_capped[gene] = {
                'go': 0, 'domain': 0, 'orthology': 0,
                'tair': 0, 'po': 0, 'synergy': 0, 'total': 0
            }
            self.gene_evidence[gene] = {
                'go_terms': [], 'go_categories': [],
                'matched_domains': [],
                'arabidopsis_orthologs': [], 'best_ortholog': None,
                'best_ortholog_score': 0.0, 'orthology_class': 'none',
                'at_symbol': '', 'at_annotation': '',
                'tair_keyword_hits': [],
                'po_terms': [], 'po_context_hits': []
            }
            self.synergy_bonuses[gene] = []
        
        # Score each evidence layer
        print("STEP 2: GO term scoring (with cap)")
        print("-" * 80)
        self._score_go_evidence(base_genes)
        print()
        
        print("STEP 3: Domain evidence scoring (with cap)")
        print("-" * 80)
        if self.config.get('domain_database', {}).get('enabled', False):
            self._score_domain_evidence(base_genes)
        else:
            print("  Domain database disabled")
        print()
        
        print("STEP 4: Orthology evidence scoring (with cap)")
        print("-" * 80)
        if self.config.get('evidence', {}).get('orthology_evidence', {}).get('enabled', False):
            self._score_orthology_evidence(base_genes)
        else:
            print("  Orthology disabled")
        print()
        
        print("STEP 5: Arabidopsis annotation keyword scoring (with cap)")
        print("-" * 80)
        if self.config.get('evidence', {}).get('arabidopsis_context', {}).get('enabled', False):
            self._score_arabidopsis_keywords(base_genes)
        else:
            print("  Arabidopsis keywords disabled")
        print()
        
        print("STEP 6: PO context scoring (with cap)")
        print("-" * 80)
        if self.config.get('evidence', {}).get('po_context', {}).get('enabled', False):
            self._score_po_context(base_genes)
        else:
            print("  PO context disabled")
        print()
        
        print("STEP 7: Synergy bonus detection")
        print("-" * 80)
        self._compute_synergy_bonuses(base_genes)
        print()
        
        print("STEP 8: Computing total scores")
        print("-" * 80)
        self._compute_total_scores(base_genes)
        print()
        
        print("STEP 9: Writing outputs")
        print("-" * 80)
        self._write_outputs(base_genes)
        print()
        
        if enable_qc and HAS_PLOTTING:
            print("STEP 10: QC diagnostics")
            print("-" * 80)
            self._generate_qc_plots(base_genes)
            print()
        
        print("="*80)
        print("Scoring completed successfully!")
        print("="*80)
    
    def _extract_base_genes(self):
        """Extract base gene list from GO terms"""
        go_file = Path(self.config['input_go_file'])
        print(f"Loading: {go_file.name}")
        go_df = pd.read_csv(go_file, sep="\t")
        print(f"  {len(go_df):,} annotations")
        
        # Collect all GO terms
        go_terms = set()
        scoring = self.config.get('scoring', {})
        go_weights = scoring.get('weights', {}).get('go', {})
        for category, terms in go_weights.items():
            if isinstance(terms, dict):
                go_terms.update(terms.keys())
        
        print(f"  Target GO terms: {len(go_terms)}")
        
        # ID mapping
        pattern = self.config['id_mapping']['protein_to_core_regex']
        suffix = self.config['id_mapping']['seurat_suffix']
        core_re = re.compile(pattern)
        
        go_df["core"] = go_df["gene"].astype(str).str.extract(core_re, expand=False)
        go_df = go_df.dropna(subset=["core"])
        go_df["seurat_gene"] = go_df["core"] + suffix
        
        # Filter
        go_level = self.config.get('go_level_filter', 'BP')
        filtered = go_df[go_df["GO"].isin(go_terms)].copy()
        if go_level and go_level.upper() != "ALL":
            filtered = filtered[filtered["level"] == go_level.upper()]
        
        self.go_annotations = filtered
        genes = sorted(set(filtered["seurat_gene"]))
        return genes
    
    def _score_go_evidence(self, genes):
        """Score GO evidence with cap"""
        scoring = self.config['scoring']
        caps = scoring['caps']
        weights = scoring['weights']['go']
        
        for gene in genes:
            annots = self.go_annotations[self.go_annotations['seurat_gene'] == gene]
            
            raw_score = 0
            go_terms = []
            categories = []
            
            for _, row in annots.iterrows():
                go_id = row['GO']
                # Check each category
                for category, term_weights in weights.items():
                    if isinstance(term_weights, dict) and go_id in term_weights:
                        raw_score += term_weights[go_id]
                        go_terms.append(go_id)
                        if category not in categories:
                            categories.append(category)
            
            capped_score = min(raw_score, caps['go_max'])
            
            self.gene_scores_raw[gene]['go'] = raw_score
            self.gene_scores_capped[gene]['go'] = capped_score
            self.gene_evidence[gene]['go_terms'] = sorted(set(go_terms))
            self.gene_evidence[gene]['go_categories'] = categories
        
        raw_scores = [self.gene_scores_raw[g]['go'] for g in genes]
        capped_scores = [self.gene_scores_capped[g]['go'] for g in genes]
        print(f"  Raw: min={min(raw_scores):.1f}, max={max(raw_scores):.1f}, mean={np.mean(raw_scores):.1f}")
        print(f"  Capped (max={caps['go_max']}): mean={np.mean(capped_scores):.1f}")
        print(f"  {sum(1 for s in raw_scores if s > caps['go_max'])} genes hit cap")
    
    def _load_arabidopsis_annotations(self):
        """Load Arabidopsis GFF3 and PO annotations"""
        arab_config = self.config.get('evidence', {}).get('arabidopsis_context', {})
        
        gff_file = arab_config.get('gff3_file')
        if gff_file:
            gff_path = Path(gff_file)
            if gff_path.exists():
                print(f"Loading GFF3: {gff_path.name}")
                self.arabidopsis_annotations = self._parse_araport_gff3(gff_path)
                print(f"  {len(self.arabidopsis_annotations)} genes")
            else:
                print(f"  WARNING: GFF3 file not found: {gff_file}")
        
        po_file = arab_config.get('po_annotation_file')
        if po_file:
            po_path = Path(po_file)
            if po_path.exists():
                print(f"Loading PO: {po_path.name}")
                self.arabidopsis_po_context = self._parse_arabidopsis_po(po_path)
                print(f"  {len(self.arabidopsis_po_context)} genes with PO context")
            else:
                print(f"  WARNING: PO file not found: {po_file}")
    
    def _parse_araport_gff3(self, gff_file):
        """Parse GFF3 - gene features only"""
        annotations = {}
        try:
            with open(gff_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 9 or parts[2] != 'gene':
                        continue
                    
                    attrs = parts[8]
                    locus_match = re.search(r'(AT[1-5MC]G\d{5})', attrs)
                    if not locus_match:
                        continue
                    locus = locus_match.group(1)
                    
                    # Extract fields
                    symbol = re.search(r'symbol=([^;]+)', attrs)
                    curator = re.search(r'curator_summary=([^;]+)', attrs)
                    full_name = re.search(r'full_name=([^;]+)', attrs)
                    comp_desc = re.search(r'computational_description=([^;]+)', attrs)
                    note = re.search(r'Note=([^;]+)', attrs)
                    
                    # Priority: curator > full_name > comp_desc > note
                    annotation_text = ''
                    if curator:
                        annotation_text = curator.group(1)
                    elif full_name:
                        annotation_text = full_name.group(1)
                    elif comp_desc:
                        annotation_text = comp_desc.group(1)
                    elif note:
                        annotation_text = note.group(1)
                    
                    # URL decode
                    annotation_text = annotation_text.replace('%2C', ',').replace('%3B', ';')
                    annotation_text = annotation_text.replace('%3D', '=').replace('%20', ' ')
                    
                    annotations[locus] = {
                        'symbol': symbol.group(1) if symbol else '',
                        'annotation_text': annotation_text
                    }
        except Exception as e:
            print(f"  WARNING: GFF3 parse error: {e}")
        return annotations
    
    def _parse_arabidopsis_po(self, po_file):
        """Parse PO annotations for senescence/PCD context"""
        po_context = {}
        keywords = ['senescence', 'senescent', 'cell death', 'hypersensitive', 'dying', 'necrosis']
        
        try:
            df = pd.read_csv(po_file, sep='\t')
            for _, row in df.iterrows():
                gene_id = str(row.get('locus_name', row.get('gene_id', '')))
                locus_match = re.search(r'(AT[1-5MC]G\d{5})', gene_id)
                if not locus_match:
                    continue
                locus = locus_match.group(1)
                
                po_term = str(row.get('term_name', row.get('po_term', '')))
                
                if locus not in po_context:
                    po_context[locus] = {'po_terms': [], 'hits': []}
                
                po_context[locus]['po_terms'].append(po_term)
                
                # Check keywords
                po_lower = po_term.lower()
                for kw in keywords:
                    if kw in po_lower:
                        po_context[locus]['hits'].append(po_term)
                        break
        except Exception as e:
            print(f"  WARNING: PO parse error: {e}")
        return po_context
    
    def _score_domain_evidence(self, genes):
        """Score domain evidence using hmmscan with cap"""
        from pgsb.domains import DomainScanner, parse_domtbl
        
        domain_config = self.config['domain_database']
        
        # Get configuration
        pfam_db = domain_config.get('pfam_database')
        protein_fasta = domain_config.get('protein_fasta')
        evalue = domain_config.get('evalue_threshold', 1e-5)
        cpu = domain_config.get('cpu', 4)
        
        if not pfam_db or not protein_fasta:
            print("  ERROR: pfam_database and protein_fasta required in config")
            return
        
        if not Path(pfam_db).exists():
            print(f"  WARNING: Pfam database not found: {pfam_db}")
            return
        
        if not Path(protein_fasta).exists():
            print(f"  WARNING: Protein FASTA not found: {protein_fasta}")
            return
        
        # Get expected domains per category
        expected_domains_config = domain_config.get('expected_domains', {})
        
        # Determine category for each gene from GO terms
        gene_categories = self._categorize_genes(genes)
        
        # Initialize scanner
        try:
            scanner = DomainScanner(pfam_db, protein_fasta, cpu=cpu, evalue=evalue)
        except FileNotFoundError as e:
            print(f"  ERROR: {e}")
            return
        
        # Run hmmscan
        print(f"  Scanning {len(genes)} genes for Pfam domains...")
        output_dir = Path(self.run_info['run_dir']) / 'domain_scan'
        domtbl_file = scanner.scan_genes(genes, output_dir=str(output_dir))
        
        if not domtbl_file or not Path(domtbl_file).exists():
            print("  ERROR: hmmscan failed")
            return
        
        # Parse domain results
        gene_domains = parse_domtbl(domtbl_file, evalue_threshold=evalue)
        
        print(f"  Found domains in {len(gene_domains)}/{len(genes)} genes")
        
        # Score based on category-specific expected domains
        scoring = self.config['scoring']
        caps = scoring['caps']
        match_score = scoring['weights'].get('domain_match', {}).get('hit', 2)
        
        all_matches = defaultdict(list)
        
        for gene in genes:
            if gene not in gene_domains:
                continue
            
            domains = gene_domains[gene]
            category = gene_categories.get(gene, 'unknown')
            
            # Get expected domains for this category
            if category == 'PCD':
                expected = expected_domains_config.get('PCD', [])
            elif category == 'ROS':
                expected = expected_domains_config.get('ROS', [])
            elif category == 'both':
                expected = list(set(
                    expected_domains_config.get('PCD', []) + 
                    expected_domains_config.get('ROS', [])
                ))
            else:
                # If no category, accept all expected domains
                expected = list(set(
                    expected_domains_config.get('PCD', []) + 
                    expected_domains_config.get('ROS', [])
                ))
            
            # Match domains
            for domain in domains:
                pfam_id = domain['pfam_id']
                if pfam_id in expected:
                    pfam_name = domain['pfam_name']
                    evalue = domain['evalue']
                    all_matches[gene].append(f"{pfam_id}({pfam_name},E={evalue:.1e})")
        
        # Score genes
        for gene in genes:
            matched = list(all_matches[gene])
            n_unique_domains = len(set(m.split('(')[0] for m in matched))
            raw_score = n_unique_domains * match_score
            capped_score = min(raw_score, caps['domain_max'])
            
            self.gene_scores_raw[gene]['domain'] = raw_score
            self.gene_scores_capped[gene]['domain'] = capped_score
            self.gene_evidence[gene]['matched_domains'] = sorted(matched)
        
        raw_scores = [self.gene_scores_raw[g]['domain'] for g in genes]
        capped_scores = [self.gene_scores_capped[g]['domain'] for g in genes]
        print(f"  Genes with expected domains: {sum(1 for s in raw_scores if s > 0)}")
        if raw_scores:
            print(f"  Raw: max={max(raw_scores):.1f}, mean={np.mean(raw_scores):.1f}")
            print(f"  Capped (max={caps['domain_max']}): mean={np.mean(capped_scores):.1f}")
    
    def _categorize_genes(self, genes):
        """Categorize genes based on GO term categories from config"""
        categories = {}
        
        # Get all GO categories from config
        go_weights = self.config.get('scoring', {}).get('weights', {}).get('go', {})
        category_terms = {}
        for category_name, terms_dict in go_weights.items():
            category_terms[category_name] = set(terms_dict.keys())
        
        for gene in genes:
            go_terms = set(self.gene_evidence[gene].get('go_terms', []))
            
            # Find which categories this gene belongs to
            gene_categories = []
            for cat_name, cat_terms in category_terms.items():
                if go_terms & cat_terms:
                    gene_categories.append(cat_name)
            
            # Assign category
            if len(gene_categories) == 0:
                categories[gene] = 'unknown'
            elif len(gene_categories) == 1:
                categories[gene] = gene_categories[0]
            else:
                # Multiple categories: join with '-'
                categories[gene] = '-'.join(sorted(gene_categories))
        
        return categories
    
    def _score_orthology_evidence(self, genes):
        """Score orthology with cap"""
        ortho_config = self.config['evidence']['orthology_evidence']
        ortho_file = Path(ortho_config['orthology_file'])
        
        print(f"  Loading: {ortho_file.name}")
        try:
            ortho_df = pd.read_csv(ortho_file, sep="\t")
        except:
            print("  WARNING: Failed to read orthology file")
            return
        
        ortho_data = self._parse_orthology_file(ortho_df, ortho_config)
        
        scoring = self.config['scoring']
        caps = scoring['caps']
        weights = scoring['weights']['orthology']
        multiplier = weights['best_hit_multiplier']
        bonus = weights['one_to_one_bonus']
        min_score = weights.get('min_best_score_for_points', 0.5)
        
        for gene in genes:
            if gene not in ortho_data:
                continue
            
            info = ortho_data[gene]
            best_score = info['best_score']
            
            if best_score < min_score:
                raw_score = 0
            else:
                raw_score = best_score * multiplier
                if info['classification'] == 'one_to_one':
                    raw_score += bonus
            
            capped_score = min(raw_score, caps['orthology_max'])
            
            self.gene_scores_raw[gene]['orthology'] = raw_score
            self.gene_scores_capped[gene]['orthology'] = capped_score
            self.gene_evidence[gene]['arabidopsis_orthologs'] = info['arabidopsis_hits']
            self.gene_evidence[gene]['best_ortholog'] = info['arabidopsis_hits'][0] if info['arabidopsis_hits'] else None
            self.gene_evidence[gene]['best_ortholog_score'] = best_score
            self.gene_evidence[gene]['orthology_class'] = info['classification']
            
            # Populate symbol/annotation from GFF3 if available
            best_orth = info['arabidopsis_hits'][0] if info['arabidopsis_hits'] else None
            if best_orth and best_orth in self.arabidopsis_annotations:
                at_info = self.arabidopsis_annotations[best_orth]
                self.gene_evidence[gene]['at_symbol'] = at_info.get('symbol', '')
                self.gene_evidence[gene]['at_annotation'] = at_info.get('annotation_text', '')[:200]
        
        raw_scores = [self.gene_scores_raw[g]['orthology'] for g in genes]
        capped_scores = [self.gene_scores_capped[gene]['orthology'] for g in genes]
        print(f"  Genes with orthologs: {sum(1 for s in raw_scores if s > 0)}")
        print(f"  Raw: max={max(raw_scores):.1f}, mean={np.mean(raw_scores):.1f}")
        print(f"  Capped (max={caps['orthology_max']}): mean={np.mean(capped_scores):.1f}")
    
    def _parse_orthology_file(self, ortho_df, config):
        """Parse inParanoid orthology file"""
        mapping_rules = config.get('mapping_rules', {})
        brachy_regex = mapping_rules.get('brachy_id_regex_extract', r"BdiBd21-3\.\dG\d{7}")
        arab_regex = mapping_rules.get('arabidopsis_id_regex_extract', r"AT[1-5MC]G\d{5}")
        
        brachy_pattern = re.compile(brachy_regex)
        arab_pattern = re.compile(arab_regex)
        suffix = self.config['id_mapping']['seurat_suffix']
        
        ortho_data = {}
        for _, row in ortho_df.iterrows():
            brachy_match = brachy_pattern.search(str(row['OrtoA']))
            if not brachy_match:
                continue
            brachy_gene = brachy_match.group() + suffix
            
            orto_b = str(row['OrtoB'])
            arab_hits = []
            arab_scores = []
            
            tokens = orto_b.split()
            for i, token in enumerate(tokens):
                arab_match = arab_pattern.search(token)
                if arab_match:
                    arab_id = arab_match.group()
                    score = float(tokens[i+1]) if i+1 < len(tokens) else 0.0
                    arab_hits.append(arab_id)
                    arab_scores.append(score)
            
            if not arab_hits:
                continue
            
            best_score = max(arab_scores)
            classification = self._classify_orthology(len(arab_hits), arab_scores)
            
            ortho_data[brachy_gene] = {
                'arabidopsis_hits': arab_hits,
                'best_score': best_score,
                'classification': classification
            }
        
        return ortho_data
    
    def _classify_orthology(self, n_hits, scores):
        """Classify orthology relationship"""
        if not scores:
            return 'none'
        best_score = max(scores)
        if best_score < 0.8:
            return 'low_confidence'
        if n_hits == 1:
            return 'one_to_one'
        if n_hits > 1:
            sorted_scores = sorted(scores, reverse=True)
            if sorted_scores[0] - sorted_scores[1] >= 0.2:
                return 'one_to_one'
        return 'one_to_many'
    
    def _score_arabidopsis_keywords(self, genes):
        """Score TAIR keywords with cap"""
        if not self.arabidopsis_annotations:
            print("  No annotations loaded")
            return
        
        arab_config = self.config['evidence']['arabidopsis_context']
        keywords = arab_config.get('keywords', [])
        if not keywords:
            print("  No keywords specified")
            return
        
        keyword_patterns = [re.compile(re.escape(kw), re.IGNORECASE) for kw in keywords]
        
        scoring = self.config['scoring']
        caps = scoring['caps']
        hit_score = scoring['weights'].get('tair_keywords', {}).get('hit', 2)
        
        for gene in genes:
            best_orth = self.gene_evidence[gene]['best_ortholog']
            if not best_orth or best_orth not in self.arabidopsis_annotations:
                continue
            
            at_info = self.arabidopsis_annotations[best_orth]
            symbol = at_info.get('symbol', '')
            annotation = at_info.get('annotation_text', '')
            
            text = f"{symbol} {annotation}"
            hits = []
            for kw, pattern in zip(keywords, keyword_patterns):
                if pattern.search(text):
                    hits.append(kw)
            
            raw_score = len(hits) * hit_score
            capped_score = min(raw_score, caps['tair_max'])
            
            self.gene_scores_raw[gene]['tair'] = raw_score
            self.gene_scores_capped[gene]['tair'] = capped_score
            self.gene_evidence[gene]['at_symbol'] = symbol
            self.gene_evidence[gene]['at_annotation'] = annotation[:200]
            self.gene_evidence[gene]['tair_keyword_hits'] = hits
        
        raw_scores = [self.gene_scores_raw[g]['tair'] for g in genes]
        capped_scores = [self.gene_scores_capped[g]['tair'] for g in genes]
        print(f"  Genes with keyword hits: {sum(1 for s in raw_scores if s > 0)}")
        print(f"  Raw: max={max(raw_scores):.1f}, mean={np.mean(raw_scores):.1f}")
        print(f"  Capped (max={caps['tair_max']}): mean={np.mean(capped_scores):.1f}")
    
    def _score_po_context(self, genes):
        """Score PO context with cap"""
        if not self.arabidopsis_po_context:
            print("  No PO context loaded")
            return
        
        scoring = self.config['scoring']
        caps = scoring['caps']
        hit_score = scoring['weights'].get('po_keywords', {}).get('hit', 1)
        
        for gene in genes:
            best_orth = self.gene_evidence[gene]['best_ortholog']
            if not best_orth or best_orth not in self.arabidopsis_po_context:
                continue
            
            po_info = self.arabidopsis_po_context[best_orth]
            po_hits = po_info.get('hits', [])
            
            raw_score = len(po_hits) * hit_score
            capped_score = min(raw_score, caps['po_max'])
            
            self.gene_scores_raw[gene]['po'] = raw_score
            self.gene_scores_capped[gene]['po'] = capped_score
            self.gene_evidence[gene]['po_terms'] = po_info.get('po_terms', [])[:5]
            self.gene_evidence[gene]['po_context_hits'] = po_hits
        
        raw_scores = [self.gene_scores_raw[g]['po'] for g in genes]
        capped_scores = [self.gene_scores_capped[g]['po'] for g in genes]
        print(f"  Genes with PO context: {sum(1 for s in raw_scores if s > 0)}")
        print(f"  Raw: max={max(raw_scores):.1f}, mean={np.mean(raw_scores):.1f}")
        print(f"  Capped (max={caps['po_max']}): mean={np.mean(capped_scores):.1f}")
    
    def _compute_synergy_bonuses(self, genes):
        """Detect synergy when independent evidence agrees"""
        scoring = self.config.get('scoring', {})
        caps = scoring.get('caps', {})
        rules = scoring.get('synergy_rules', [])
        
        if not rules:
            print("  No synergy rules defined")
            return
        
        for gene in genes:
            total_synergy = 0
            
            for rule in rules:
                conditions = rule.get('if_all', [])
                bonus = rule.get('bonus', 0)
                
                all_met = True
                for cond in conditions:
                    if not self._evaluate_condition(gene, cond):
                        all_met = False
                        break
                
                if all_met:
                    total_synergy += bonus
                    self.synergy_bonuses[gene].append(rule.get('name', 'unnamed'))
            
            # Cap synergy
            capped_synergy = min(total_synergy, caps.get('synergy_max', 6))
            self.gene_scores_raw[gene]['synergy'] = total_synergy
            self.gene_scores_capped[gene]['synergy'] = capped_synergy
        
        raw_scores = [self.gene_scores_raw[g]['synergy'] for g in genes]
        capped_scores = [self.gene_scores_capped[g]['synergy'] for g in genes]
        print(f"  Genes with synergy bonuses: {sum(1 for s in raw_scores if s > 0)}")
        print(f"  Raw: max={max(raw_scores):.1f}, mean={np.mean(raw_scores):.1f}")
        print(f"  Capped (max={caps.get('synergy_max', 6)}): mean={np.mean(capped_scores):.1f}")
    
    def _evaluate_condition(self, gene, condition):
        """Evaluate synergy condition"""
        evidence = self.gene_evidence[gene]
        
        if condition == "has_PCD_GO":
            return 'PCD' in evidence['go_categories']
        elif condition == "has_ROS_GO":
            return 'ROS' in evidence['go_categories']
        elif condition == "has_expected_domain":
            return len(evidence['matched_domains']) > 0
        elif condition.startswith("tair_keyword_hits>="):
            threshold = int(condition.split(">=")[1])
            return len(evidence['tair_keyword_hits']) >= threshold
        elif condition == "has_one_to_one_ortholog":
            return evidence['orthology_class'] == 'one_to_one'
        return False
    
    def _compute_total_scores(self, genes):
        """Sum all capped scores"""
        for gene in genes:
            total = sum(self.gene_scores_capped[gene][k] for k in ['go', 'domain', 'orthology', 'tair', 'po', 'synergy'])
            self.gene_scores_capped[gene]['total'] = total
        
        scores = [self.gene_scores_capped[g]['total'] for g in genes]
        print(f"  Total scores: min={min(scores):.1f}, max={max(scores):.1f}")
        print(f"  Mean={np.mean(scores):.1f}, median={np.median(scores):.1f}")
    
    def _write_outputs(self, genes):
        """Write all output files to run directory"""
        outputs_dir = self.run_info['outputs']
        
        # Base list (IDs only)
        base_file = outputs_dir / f"{self.output_prefix}_genes.base.txt"
        base_file.write_text("\n".join(genes) + "\n")
        print(f"✓ {base_file.name} ({len(genes)} genes)")
        
        # Scored table (complete)
        rows = []
        for gene in genes:
            raw = self.gene_scores_raw[gene]
            capped = self.gene_scores_capped[gene]
            evidence = self.gene_evidence[gene]
            
            rows.append({
                'brachy_gene_id': gene,
                'go_terms': ', '.join(evidence['go_terms']),
                'go_score_raw': raw['go'],
                'go_score_capped': capped['go'],
                'domain_hits': ', '.join(evidence['matched_domains']),
                'domain_score_raw': raw['domain'],
                'domain_score_capped': capped['domain'],
                'best_arabidopsis_id': evidence['best_ortholog'] or '',
                'best_orth_score': evidence['best_ortholog_score'],
                'orth_class': evidence['orthology_class'],
                'orth_score_raw': raw['orthology'],
                'orth_score_capped': capped['orthology'],
                'tair_symbol': evidence['at_symbol'],
                'tair_annotation_text': evidence['at_annotation'],
                'tair_keyword_hits': ', '.join(evidence['tair_keyword_hits']),
                'tair_score_raw': raw['tair'],
                'tair_score_capped': capped['tair'],
                'po_hits': ', '.join(evidence['po_context_hits']),
                'po_score_raw': raw['po'],
                'po_score_capped': capped['po'],
                'synergy_bonus': capped['synergy'],
                'synergy_triggers': ', '.join(self.synergy_bonuses[gene]),
                'total_score': capped['total']
            })
        
        df_scored = pd.DataFrame(rows).sort_values('total_score', ascending=False)
        scored_file = outputs_dir / f"{self.output_prefix}_genes.scored.tsv"
        df_scored.to_csv(scored_file, sep="\t", index=False)
        print(f"✓ {scored_file.name}")
        
        # Create enriched ALL overview table
        self._create_overview_table(df_scored, outputs_dir, suffix="ALL")
        
        # High-confidence selection
        high_confidence_genes, high_conf_df = self._select_high_confidence(genes, df_scored, outputs_dir)
        
        # Create enriched HIGH overview table
        if len(high_confidence_genes) > 0:
            self._create_overview_table(high_conf_df, outputs_dir, suffix="HIGH")
            # Create light HIGH category file
            self._create_light_category_file(high_conf_df, outputs_dir)
        
        # Generate manifest
        self._generate_manifest(outputs_dir, len(genes), len(high_confidence_genes))
    
    def _select_high_confidence(self, genes, df, outputs_dir):
        """Select high-confidence genes and return both list and dataframe"""
        selection = self.config.get('selection', {})
        mode = selection.get('mode', 'knee')
        
        scores = df['total_score'].values
        
        if mode == 'knee':
            threshold = self._find_knee_threshold(scores)
            selected_df = df[df['total_score'] >= threshold]
            selected = selected_df['brachy_gene_id'].tolist()
            method_desc = f"knee detection (threshold: {threshold:.1f})"
        else:
            quantile = selection.get('quantile', 0.90)
            threshold = np.quantile(scores, quantile)
            selected_df = df[df['total_score'] >= threshold]
            selected = selected_df['brachy_gene_id'].tolist()
            method_desc = f"{quantile*100:.0f}th percentile (threshold: {threshold:.1f})"
        
        # Write IDs file
        hc_file = outputs_dir / f"{self.output_prefix}_genes.high_confidence.txt"
        hc_file.write_text("\n".join(selected) + "\n")
        print(f"✓ {hc_file.name} ({len(selected)} genes)")
        
        # Write summary
        summary_file = outputs_dir / f"{self.output_prefix}_genes.high_confidence.summary.txt"
        with open(summary_file, 'w') as f:
            f.write("HIGH-CONFIDENCE GENE SELECTION\n")
            f.write("="*80 + "\n\n")
            f.write(f"Method: {method_desc}\n")
            f.write(f"Total genes: {len(genes)}\n")
            f.write(f"High-confidence: {len(selected)} ({len(selected)/len(genes)*100:.1f}%)\n\n")
            f.write(f"Score stats (all genes):\n")
            f.write(f"  Min: {scores.min():.1f}, Max: {scores.max():.1f}\n")
            f.write(f"  Mean: {scores.mean():.1f}, Median: {np.median(scores):.1f}\n")
        
        print(f"✓ {summary_file.name}")
        
        return selected, selected_df
    
    def _find_knee_threshold(self, scores):
        """Find knee point in score distribution"""
        sorted_scores = np.sort(scores)[::-1]
        if len(sorted_scores) < 10:
            return np.median(scores)
        
        # Normalize
        norm = (sorted_scores - sorted_scores.min()) / (sorted_scores.max() - sorted_scores.min() + 1e-10)
        
        # Second derivative
        grad = np.gradient(norm)
        grad2 = np.gradient(grad)
        
        # Find max curvature in top 20-30%
        start = int(len(sorted_scores) * 0.1)
        end = int(len(sorted_scores) * 0.3)
        knee_idx = start + np.argmax(np.abs(grad2[start:end]))
        
        return sorted_scores[knee_idx]
    
    def _generate_qc_plots(self, genes):
        """Generate QC diagnostic plots"""
        if not HAS_PLOTTING:
            print("  matplotlib not available, skipping plots")
            return
        
        qc_dir = self.run_info['qc']
        
        df = pd.DataFrame([
            {
                'gene': g,
                'go_capped': self.gene_scores_capped[g]['go'],
                'domain_capped': self.gene_scores_capped[g]['domain'],
                'total': self.gene_scores_capped[g]['total'],
                'has_domain': len(self.gene_evidence[g]['matched_domains']) > 0,
                'has_tair_kw': len(self.gene_evidence[g]['tair_keyword_hits']) > 0
            }
            for g in genes
        ])
        
        # Plot 1: GO vs domain
        plt.figure(figsize=(8, 6))
        plt.scatter(df['go_capped'], df['domain_capped'], alpha=0.5)
        plt.xlabel('GO Score (capped)')
        plt.ylabel('Domain Score (capped)')
        plt.title('GO vs Domain Evidence')
        plt.savefig(qc_dir / f"{self.output_prefix}_go_vs_domain.png", dpi=150, bbox_inches='tight')
        plt.close()
        
        # Plot 2: Total vs GO
        plt.figure(figsize=(8, 6))
        plt.scatter(df['total'], df['go_capped'], alpha=0.5)
        plt.xlabel('Total Score')
        plt.ylabel('GO Score (capped)')
        plt.title('Total Score vs GO Evidence')
        plt.savefig(qc_dir / f"{self.output_prefix}_total_vs_go.png", dpi=150, bbox_inches='tight')
        plt.close()
        
        # Plot 3: Score distribution
        plt.figure(figsize=(10, 6))
        plt.hist(df['total'], bins=30, edgecolor='black')
        plt.xlabel('Total Score')
        plt.ylabel('Gene Count')
        plt.title('Total Score Distribution')
        plt.savefig(qc_dir / f"{self.output_prefix}_score_distribution.png", dpi=150, bbox_inches='tight')
        plt.close()
        
        # Text report
        top20 = df.nlargest(20, 'total')
        bottom20 = df.nsmallest(20, 'total')
        
        report_file = qc_dir / f"{self.output_prefix}_qc_summary.txt"
        with open(report_file, 'w') as f:
            f.write("QC DIAGNOSTIC REPORT\n")
            f.write("="*80 + "\n\n")
            f.write(f"Top 20 genes:\n")
            f.write(f"  With expected domains: {top20['has_domain'].sum()}/20 ({top20['has_domain'].sum()/20*100:.0f}%)\n")
            f.write(f"  With TAIR keyword hits: {top20['has_tair_kw'].sum()}/20 ({top20['has_tair_kw'].sum()/20*100:.0f}%)\n\n")
            f.write(f"Bottom 20 genes:\n")
            f.write(f"  Mean GO score: {bottom20['go_capped'].mean():.1f}\n")
            f.write(f"  Mean total score: {bottom20['total'].mean():.1f}\n\n")
            f.write(f"Correlation analysis:\n")
            f.write(f"  Total vs GO: r={df['total'].corr(df['go_capped']):.3f}\n")
            f.write(f"  GO vs Domain: r={df['go_capped'].corr(df['domain_capped']):.3f}\n")
        
        print(f"✓ QC plots in {qc_dir}/")
        print(f"✓ {report_file}")
    
    def _create_light_category_file(self, df_high_conf, outputs_dir):
        """
        Create light HIGH confidence file with just gene ID and category
        
        Args:
            df_high_conf: DataFrame with high-confidence genes
            outputs_dir: Output directory path
        """
        light_rows = []
        
        # Determine category from GO terms
        gene_categories = self._categorize_genes(df_high_conf['brachy_gene_id'].tolist())
        
        for _, row in df_high_conf.iterrows():
            gene = row['brachy_gene_id']
            category = gene_categories.get(gene, 'unknown')
            
            light_rows.append({
                'gene_id': gene,
                'category': category
            })
        
        df_light = pd.DataFrame(light_rows)
        light_file = outputs_dir / f"{self.output_prefix}_genes.HIGH_category.tsv"
        df_light.to_csv(light_file, sep="\t", index=False)
        print(f"✓ {light_file.name} ({len(df_light)} genes)")
        
        # Print category counts
        category_counts = df_light['category'].value_counts()
        print(f"  Categories: {', '.join([f'{cat}={count}' for cat, count in category_counts.items()])}")
        
        return light_file
    
    def _create_overview_table(self, df_scored, outputs_dir, suffix="ALL"):
        """
        Create enriched overview table with aggregated annotations
        
        Args:
            df_scored: DataFrame with all scored genes
            outputs_dir: Output directory path
            suffix: "ALL" or "HIGH" for filename
        """
        overview_rows = []
        for _, row in df_scored.iterrows():
            gene = row['brachy_gene_id']
            evidence = self.gene_evidence[gene]
            
            # Create evidence summary
            evidence_parts = []
            if evidence['go_categories']:
                evidence_parts.append("+".join(evidence['go_categories']) + "_GO")
            if evidence['matched_domains']:
                evidence_parts.append("Domain")
            if evidence['tair_keyword_hits']:
                evidence_parts.append("TAIR_keywords")
            if evidence['orthology_class'] == 'one_to_one':
                evidence_parts.append(f"1:1_{evidence['best_ortholog']}")
            elif evidence['orthology_class'] == 'one_to_many':
                evidence_parts.append(f"1:many_ortho")
            
            evidence_summary = "; ".join(evidence_parts) if evidence_parts else "GO_only"
            
            # Arabidopsis orthologs (all hits, comma-separated)
            arabidopsis_hits = ", ".join(evidence['arabidopsis_orthologs'])
            
            overview_rows.append({
                'brachy_gene_id': gene,
                'total_score': row['total_score'],
                'go_score': row['go_score_capped'],
                'domain_score': row['domain_score_capped'],
                'orth_score': row['orth_score_capped'],
                'tair_score': row['tair_score_capped'],
                'po_score': row['po_score_capped'],
                'synergy_bonus': row['synergy_bonus'],
                'go_terms': row['go_terms'],
                'domain_hits': row['domain_hits'],
                'best_arabidopsis_id': row['best_arabidopsis_id'],
                'best_orth_score': row['best_orth_score'],
                'orth_class': row['orth_class'],
                'arabidopsis_hits': arabidopsis_hits,
                'tair_symbol': row['tair_symbol'],
                'tair_annotation': row['tair_annotation_text'],
                'po_stage_hits': row['po_hits'],
                'evidence_summary': evidence_summary
            })
        
        df_overview = pd.DataFrame(overview_rows)
        overview_file = outputs_dir / f"{self.output_prefix}_genes.{suffix}_overview.tsv"
        df_overview.to_csv(overview_file, sep="\t", index=False)
        print(f"✓ {overview_file.name} ({len(df_overview)} genes)")
        
        return overview_file
    
    def _generate_manifest(self, outputs_dir, n_total, n_high_conf):
        """Generate manifest.json with all output files"""
        files_info = [
            {
                'path': f"outputs/{self.output_prefix}_genes.base.txt",
                'description': 'Base gene list (IDs only)',
                'rows': n_total
            },
            {
                'path': f"outputs/{self.output_prefix}_genes.scored.tsv",
                'description': 'Complete scored table with all evidence',
                'rows': n_total
            },
            {
                'path': f"outputs/{self.output_prefix}_genes.ALL_overview.tsv",
                'description': 'Enriched overview table for all genes',
                'rows': n_total
            },
            {
                'path': f"outputs/{self.output_prefix}_genes.high_confidence.txt",
                'description': 'High-confidence genes (IDs only)',
                'rows': n_high_conf
            },
            {
                'path': f"outputs/{self.output_prefix}_genes.HIGH_overview.tsv",
                'description': 'Enriched overview table for high-confidence genes',
                'rows': n_high_conf
            },
            {
                'path': f"outputs/{self.output_prefix}_genes.high_confidence.summary.txt",
                'description': 'Selection method and statistics'
            },
            {
                'path': 'config_used.yaml',
                'description': 'Configuration snapshot for this run'
            }
        ]
        
        manifest_path = create_manifest(files_info, self.run_info['run_dir'])
        print(f"✓ manifest.json")


def main():
    parser = argparse.ArgumentParser(
        description='Plant Gene Signature Ranker v3.1 - Capped Multi-Evidence Scoring + Run Management'
    )
    parser.add_argument('--config', '-c', required=True, help='YAML config file')
    parser.add_argument('--qc', action='store_true', help='Generate QC plots')
    parser.add_argument('--run_name', type=str, help='Optional custom run name suffix')
    parser.add_argument('--overwrite', action='store_true', help='Allow overwriting existing run directory')
    
    args = parser.parse_args()
    
    ranker = GeneSignatureRanker(args.config, run_name=args.run_name, overwrite=args.overwrite)
    ranker.run(enable_qc=args.qc)
    
    # Update latest symlink
    update_latest_symlink(ranker.run_info['run_dir'])
    
    print(f"\n✓ Run completed: {ranker.run_info['run_id']}")
    print(f"  Outputs: {ranker.run_info['outputs']}")


if __name__ == "__main__":
    main()
