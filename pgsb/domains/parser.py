"""
Parse HMMER domtblout format
"""
from typing import List, Dict, Set
from collections import defaultdict


def parse_domtbl(domtbl_file: str, evalue_threshold: float = 1e-5) -> Dict[str, List[Dict]]:
    """
    Parse hmmscan domain table output
    
    Args:
        domtbl_file: Path to domtblout file
        evalue_threshold: E-value threshold for filtering hits
        
    Returns:
        Dictionary mapping gene_id -> list of domain hits
        Each domain hit is a dict with keys: pfam_id, pfam_name, evalue, score
    """
    # Ensure evalue_threshold is float
    evalue_threshold = float(evalue_threshold)
    
    gene_domains = defaultdict(list)
    
    with open(domtbl_file, 'r') as f:
        for line in f:
            # Skip comments
            if line.startswith('#'):
                continue
            
            parts = line.strip().split()
            if len(parts) < 23:
                continue
            
            # Parse domtblout format (space-separated, some fields may have spaces)
            # Format: target_name accession tlen query_name accession qlen full_evalue ...
            target_name = parts[0]      # Pfam domain name (e.g., Peptidase_C14)
            target_accession = parts[1]  # Pfam ID (e.g., PF00656.26)
            query_name = parts[3]        # Protein ID
            full_evalue = float(parts[6])  # Full sequence E-value
            full_score = float(parts[7])   # Full sequence score
            domain_evalue = float(parts[12])  # Best domain E-value
            domain_score = float(parts[13])   # Best domain score
            
            # Use domain E-value for filtering (more stringent)
            if domain_evalue > evalue_threshold:
                continue
            
            # Extract clean Pfam ID (remove version)
            pfam_id = target_accession.split('.')[0]  # PF00656
            
            # Extract gene ID from protein ID
            # BdiBd21-3.2G0277200.1.p -> BdiBd21-3.2G0277200.v1.2
            gene_id = extract_gene_id(query_name)
            
            domain_hit = {
                'pfam_id': pfam_id,
                'pfam_accession': target_accession,
                'pfam_name': target_name,
                'evalue': domain_evalue,
                'score': domain_score,
                'full_evalue': full_evalue,
                'full_score': full_score
            }
            
            gene_domains[gene_id].append(domain_hit)
    
    return dict(gene_domains)


def extract_gene_id(protein_id: str) -> str:
    """
    Convert protein ID to gene ID
    
    Examples:
        BdiBd21-3.2G0277200.1.p -> BdiBd21-3.2G0277200.v1.2
        BdiBd21-3.2G0277200.1 -> BdiBd21-3.2G0277200.v1.2
    """
    # Remove protein suffix
    base = protein_id.replace('.1.p', '').replace('.p', '')
    
    # Add version suffix if not present
    if not base.endswith('.v1.2'):
        base = f"{base}.v1.2"
    
    return base


def filter_domains_by_list(gene_domains: Dict[str, List[Dict]], 
                           expected_domains: List[str]) -> Dict[str, List[Dict]]:
    """
    Filter domain hits to only include expected domains
    
    Args:
        gene_domains: Output from parse_domtbl
        expected_domains: List of Pfam IDs to keep (e.g., ['PF00656', 'PF00931'])
        
    Returns:
        Filtered dictionary with only expected domains
    """
    filtered = {}
    
    for gene_id, domains in gene_domains.items():
        matched = [d for d in domains if d['pfam_id'] in expected_domains]
        if matched:
            filtered[gene_id] = matched
    
    return filtered


def get_domain_summary(gene_domains: Dict[str, List[Dict]]) -> Dict[str, Set[str]]:
    """
    Get summary of unique Pfam IDs per gene
    
    Args:
        gene_domains: Output from parse_domtbl
        
    Returns:
        Dictionary mapping gene_id -> set of Pfam IDs
    """
    summary = {}
    
    for gene_id, domains in gene_domains.items():
        pfam_ids = {d['pfam_id'] for d in domains}
        summary[gene_id] = pfam_ids
    
    return summary
