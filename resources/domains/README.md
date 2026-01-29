# Protein Domain Resources

This directory contains protein domain databases for functional annotation.

## Pfam Database

**Location:** `pfam/Pfam-A.hmm`  
**Version:** Pfam 36.0  
**Source:** [EBI Pfam FTP](https://ftp.ebi.ac.uk/pub/databases/Pfam/)  
**HMMs:** 20,795 domain models  
**Status:** ✓ Indexed with hmmpress

### Files
- `Pfam-A.hmm` - HMM profiles (881 MB uncompressed)
- `Pfam-A.hmm.h3m` - Binary HMM file
- `Pfam-A.hmm.h3i` - SSI index
- `Pfam-A.hmm.h3f` - MSV filter profiles
- `Pfam-A.hmm.h3p` - Remaining profiles

### Usage

**Scan proteins against Pfam:**
```bash
hmmscan --cpu 4 --domtblout output.domtbl resources/domains/pfam/Pfam-A.hmm proteins.fa
```

**Expected output format (domtblout):**
- Column 1: Target domain (Pfam ID, e.g., PF00931)
- Column 4: Query protein ID
- Column 7: E-value (full sequence)
- Column 13: E-value (best domain)
- And more...

### Relevant Domains for PCD/ROS Genes

**Programmed Cell Death (PCD):**
- PF00656 - ICE-like protease (caspase/metacaspase)
- PF00931 - NB-ARC domain (plant R genes)
- PF01582 - TIR domain (pathogen recognition)
- PF12799, PF13855 - LRR (leucine-rich repeat)

**Reactive Oxygen Species (ROS):**
- PF00141 - Peroxidase
- PF00199 - Catalase
- PF00067 - Cytochrome P450
- PF01764 - Lipase/lipooxygenase
- PF08030 - Ferric reductase NAD binding

## Running Domain Annotation

The pipeline will automatically run `hmmscan` on selected genes when domain scoring is enabled in the configuration.

See `config_PCD_ROS_scoring.yaml` for configuration options.
