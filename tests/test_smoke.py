"""
Smoke test for complete pipeline with minimal data
"""

import pytest
import tempfile
import shutil
from pathlib import Path
import sys

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from rank_gene_signatures import GeneSignatureRanker


def test_smoke_test_minimal():
    """Run complete pipeline with minimal test data"""
    
    # Use test config
    config_path = Path(__file__).parent / "data" / "config_test.yaml"
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create ranker with temporary base_dir
        # Note: We can't easily override base_dir, so we'll use default "results"
        # and clean up after
        
        ranker = GeneSignatureRanker(config_path, run_name="smoke_test")
        ranker.run(enable_qc=False)
        
        # Check outputs exist
        outputs_dir = ranker.run_info['outputs']
        
        assert (outputs_dir / f"{ranker.output_prefix}_genes.base.txt").exists()
        assert (outputs_dir / f"{ranker.output_prefix}_genes.scored.tsv").exists()
        assert (outputs_dir / f"{ranker.output_prefix}_genes.ALL_overview.tsv").exists()
        assert (outputs_dir / f"{ranker.output_prefix}_genes.high_confidence.txt").exists()
        assert (outputs_dir / f"{ranker.output_prefix}_genes.HIGH_overview.tsv").exists()
        assert (outputs_dir / "manifest.json").exists()
        
        # Check config snapshot
        assert (ranker.run_info['run_dir'] / "config_used.yaml").exists()
        
        # Read base file and verify we got some genes
        base_file = outputs_dir / f"{ranker.output_prefix}_genes.base.txt"
        genes = base_file.read_text().strip().split('\n')
        
        assert len(genes) > 0, "Should have at least some genes"
        assert all('.v1.2' in g for g in genes), "All genes should have .v1.2 suffix"
        
        # Clean up test run
        shutil.rmtree(ranker.run_info['run_dir'])


def test_overview_table_columns():
    """Verify overview tables have required columns"""
    
    config_path = Path(__file__).parent / "data" / "config_test.yaml"
    
    ranker = GeneSignatureRanker(config_path, run_name="smoke_columns")
    ranker.run(enable_qc=False)
    
    import pandas as pd
    
    outputs_dir = ranker.run_info['outputs']
    overview_file = outputs_dir / f"{ranker.output_prefix}_genes.ALL_overview.tsv"
    
    df = pd.read_csv(overview_file, sep='\t')
    
    # Required columns
    required_cols = [
        'brachy_gene_id', 'total_score', 'go_score', 'domain_score',
        'orth_score', 'tair_score', 'po_score', 'synergy_bonus',
        'go_terms', 'domain_hits', 'best_arabidopsis_id',
        'orth_class', 'evidence_summary'
    ]
    
    for col in required_cols:
        assert col in df.columns, f"Missing required column: {col}"
    
    # Clean up
    shutil.rmtree(ranker.run_info['run_dir'])


def test_scoring_caps_applied():
    """Verify scoring caps are enforced"""
    
    config_path = Path(__file__).parent / "data" / "config_test.yaml"
    
    ranker = GeneSignatureRanker(config_path, run_name="smoke_caps")
    ranker.run(enable_qc=False)
    
    # Check that no gene exceeds go_max
    go_max = ranker.config['scoring']['caps']['go_max']
    
    for gene, scores in ranker.gene_scores_capped.items():
        assert scores['go'] <= go_max, f"GO score {scores['go']} exceeds cap {go_max}"
        assert scores['domain'] <= ranker.config['scoring']['caps']['domain_max']
        assert scores['orthology'] <= ranker.config['scoring']['caps']['orthology_max']
    
    # Clean up
    shutil.rmtree(ranker.run_info['run_dir'])


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
