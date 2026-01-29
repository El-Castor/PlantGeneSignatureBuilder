"""
Protein domain scanner using HMMER/hmmscan
"""
import os
import subprocess
import tempfile
from pathlib import Path
from typing import List, Dict, Set
import re


class DomainScanner:
    """Scan protein sequences for Pfam domains using hmmscan"""
    
    def __init__(self, pfam_db: str, protein_fasta: str, cpu: int = 4, evalue: float = 1e-5):
        """
        Initialize domain scanner
        
        Args:
            pfam_db: Path to Pfam-A.hmm database
            protein_fasta: Path to genome protein FASTA file
            cpu: Number of CPUs for hmmscan
            evalue: E-value threshold for domain hits
        """
        self.pfam_db = pfam_db
        self.protein_fasta = protein_fasta
        self.cpu = cpu
        self.evalue = evalue
        
        # Validate files exist
        if not os.path.exists(pfam_db):
            raise FileNotFoundError(f"Pfam database not found: {pfam_db}")
        if not os.path.exists(protein_fasta):
            raise FileNotFoundError(f"Protein FASTA not found: {protein_fasta}")
        
        # Check for indexed Pfam database
        required_indices = [f"{pfam_db}.h3{ext}" for ext in ['m', 'i', 'f', 'p']]
        missing = [idx for idx in required_indices if not os.path.exists(idx)]
        if missing:
            raise FileNotFoundError(
                f"Pfam database not indexed. Run: hmmpress {pfam_db}\n"
                f"Missing files: {missing}"
            )
    
    def extract_proteins(self, gene_ids: List[str], output_fasta: str) -> int:
        """
        Extract protein sequences for specific genes from genome FASTA
        
        Args:
            gene_ids: List of gene IDs to extract (e.g., BdiBd21-3.2G0277200.v1.2)
            output_fasta: Path to output FASTA file
            
        Returns:
            Number of proteins extracted
        """
        print(f"    Extracting proteins from {os.path.basename(self.protein_fasta)}...")
        
        # Convert gene IDs to protein ID patterns
        protein_patterns = set()
        for gene_id in gene_ids:
            # Try multiple protein ID patterns
            base_id = gene_id.replace('.v1.2', '')
            protein_patterns.add(f"{base_id}.1.p")  # Standard pattern
            protein_patterns.add(f"{base_id}.1")     # Alternative
            protein_patterns.add(gene_id)             # Exact match
        
        extracted = 0
        current_seq = []
        writing = False
        processed_lines = 0
        
        with open(self.protein_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
            for line in infile:
                processed_lines += 1
                if processed_lines % 100000 == 0:
                    print(f"      Processed {processed_lines//1000}K lines, extracted {extracted} proteins...")
                
                if line.startswith('>'):
                    # Save previous sequence if we were writing
                    if writing and current_seq:
                        outfile.write(''.join(current_seq))
                        extracted += 1
                    
                    # Check if this protein matches our gene list
                    header = line.strip()
                    protein_id = header.split()[0][1:]  # Remove '>'
                    
                    # Check if any pattern matches
                    writing = any(pattern in protein_id for pattern in protein_patterns)
                    
                    if writing:
                        outfile.write(line)
                        current_seq = []
                else:
                    if writing:
                        current_seq.append(line)
            
            # Write last sequence
            if writing and current_seq:
                outfile.write(''.join(current_seq))
                extracted += 1
        
        print(f"    ✓ Extracted {extracted}/{len(gene_ids)} proteins")
        return extracted
    
    def run_hmmscan(self, query_fasta: str, output_domtbl: str) -> bool:
        """
        Run hmmscan on protein sequences
        
        Args:
            query_fasta: Input protein FASTA file
            output_domtbl: Output domain table file
            
        Returns:
            True if successful
        """
        cmd = [
            'hmmscan',
            '--cpu', str(self.cpu),
            '--domtblout', output_domtbl,
            '-E', str(self.evalue),
            '--domE', str(self.evalue),
            self.pfam_db,
            query_fasta
        ]
        
        try:
            # Run hmmscan, suppress stdout
            result = subprocess.run(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            return True
        except subprocess.CalledProcessError as e:
            print(f"Error running hmmscan: {e.stderr}")
            return False
    
    def scan_genes(self, gene_ids: List[str], output_dir: str = None) -> str:
        """
        Complete workflow: extract proteins and scan for domains
        
        Args:
            gene_ids: List of gene IDs to scan
            output_dir: Directory for output files (default: temp directory)
            
        Returns:
            Path to domain table output file
        """
        if output_dir is None:
            output_dir = tempfile.mkdtemp(prefix='hmmscan_')
        else:
            os.makedirs(output_dir, exist_ok=True)
        
        # File paths
        protein_fasta = os.path.join(output_dir, 'query_proteins.fa')
        domtbl_file = os.path.join(output_dir, 'domains.domtbl')
        
        # Extract proteins
        n_extracted = self.extract_proteins(gene_ids, protein_fasta)
        
        if n_extracted == 0:
            print(f"    WARNING: No proteins extracted for {len(gene_ids)} genes")
            return None
        
        # Run hmmscan
        print(f"    Running hmmscan (E-value <= {self.evalue}, {self.cpu} CPUs)...")
        success = self.run_hmmscan(protein_fasta, domtbl_file)
        
        if not success:
            return None
        
        return domtbl_file
