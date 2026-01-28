"""
Run directory management for reproducible analysis runs

Each pipeline execution creates a unique run directory with:
- Timestamp (YYYY-MM-DD_HHMMSS)
- Output prefix from config
- Config content hash (first 8 chars of SHA1)
- Optional custom run name

Example: results/runs/2026-01-28_143501__PCD_ROS__abc12345/
"""

import os
import sys
import hashlib
import json
import shutil
from pathlib import Path
from datetime import datetime
import yaml


def compute_config_hash(config_dict):
    """
    Compute stable hash from config content
    
    Args:
        config_dict: Configuration dictionary
    
    Returns:
        str: First 8 characters of SHA1 hash
    """
    # Serialize config deterministically
    config_str = json.dumps(config_dict, sort_keys=True)
    hash_obj = hashlib.sha1(config_str.encode('utf-8'))
    return hash_obj.hexdigest()[:8]


def create_run_dir(config_dict, output_prefix, run_name=None, overwrite=False, base_dir="results"):
    """
    Create a new run directory with deterministic naming
    
    Args:
        config_dict: Full configuration dictionary
        output_prefix: Prefix from config['output_prefix']
        run_name: Optional custom suffix for run name
        overwrite: If True, allow reusing existing run_dir (use with caution)
        base_dir: Base directory for all runs (default: "results")
    
    Returns:
        dict: Paths {
            'run_dir': Path to run directory,
            'outputs': Path to outputs/ subdirectory,
            'qc': Path to qc/ subdirectory,
            'logs': Path to logs/ subdirectory,
            'run_id': Run identifier string
        }
    """
    # Generate run identifier
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    config_hash = compute_config_hash(config_dict)
    
    # Build run_id: timestamp__prefix__hash[__custom]
    run_id = f"{timestamp}__{output_prefix}__{config_hash}"
    if run_name:
        run_id += f"__{run_name}"
    
    # Create directory structure
    base_path = Path(base_dir)
    runs_dir = base_path / "runs"
    runs_dir.mkdir(parents=True, exist_ok=True)
    
    run_dir = runs_dir / run_id
    
    # Check if exists
    if run_dir.exists() and not overwrite:
        raise FileExistsError(
            f"Run directory already exists: {run_dir}\n"
            f"Use --overwrite to reuse, or change config/run_name to create new run."
        )
    
    # Create subdirectories
    run_dir.mkdir(exist_ok=overwrite)
    outputs_dir = run_dir / "outputs"
    qc_dir = run_dir / "qc"
    logs_dir = run_dir / "logs"
    
    outputs_dir.mkdir(exist_ok=True)
    qc_dir.mkdir(exist_ok=True)
    logs_dir.mkdir(exist_ok=True)
    
    return {
        'run_dir': run_dir,
        'outputs': outputs_dir,
        'qc': qc_dir,
        'logs': logs_dir,
        'run_id': run_id
    }


def write_config_snapshot(config_dict, run_dir):
    """
    Save config used for this run as YAML
    
    Args:
        config_dict: Configuration dictionary
        run_dir: Path to run directory
    """
    config_path = Path(run_dir) / "config_used.yaml"
    with open(config_path, 'w') as f:
        yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)
    return config_path


def update_latest_symlink(run_dir, base_dir="results"):
    """
    Update 'latest' symlink to point to most recent run
    
    On Windows: creates latest_path.txt instead of symlink
    
    Args:
        run_dir: Path to the new run directory
        base_dir: Base directory containing runs/
    """
    base_path = Path(base_dir)
    latest_link = base_path / "latest"
    run_path = Path(run_dir)
    
    # Remove existing symlink/file
    if latest_link.exists() or latest_link.is_symlink():
        latest_link.unlink()
    
    # Try to create symlink (Unix/Mac)
    try:
        # Use relative path for portability
        relative_path = os.path.relpath(run_path, base_path)
        latest_link.symlink_to(relative_path)
    except (OSError, NotImplementedError):
        # Windows fallback: write text file
        latest_txt = base_path / "latest_path.txt"
        with open(latest_txt, 'w') as f:
            f.write(str(run_path.absolute()) + "\n")


def create_manifest(files_info, run_dir):
    """
    Create manifest.json documenting all outputs
    
    Args:
        files_info: List of dicts with keys:
            - path: relative path from run_dir
            - description: file description
            - rows: number of rows (if tabular)
        run_dir: Path to run directory
    
    Returns:
        Path to manifest.json
    """
    run_path = Path(run_dir)
    manifest = {
        'run_id': run_path.name,
        'created': datetime.now().isoformat(),
        'files': []
    }
    
    for file_info in files_info:
        file_path = run_path / file_info['path']
        
        if not file_path.exists():
            continue
        
        # Compute checksum
        sha256_hash = hashlib.sha256()
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b""):
                sha256_hash.update(chunk)
        
        manifest['files'].append({
            'path': file_info['path'],
            'description': file_info.get('description', ''),
            'rows': file_info.get('rows'),
            'size_bytes': file_path.stat().st_size,
            'sha256': sha256_hash.hexdigest()
        })
    
    manifest_path = run_path / "outputs" / "manifest.json"
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    
    return manifest_path


def get_latest_run_dir(base_dir="results"):
    """
    Get path to latest run directory
    
    Args:
        base_dir: Base directory containing runs/
    
    Returns:
        Path to latest run directory, or None if no runs exist
    """
    base_path = Path(base_dir)
    latest_link = base_path / "latest"
    latest_txt = base_path / "latest_path.txt"
    
    # Try symlink first
    if latest_link.is_symlink() or (latest_link.exists() and latest_link.is_dir()):
        return latest_link.resolve()
    
    # Try text file (Windows)
    if latest_txt.exists():
        with open(latest_txt, 'r') as f:
            path_str = f.read().strip()
            return Path(path_str)
    
    # Fallback: find most recent by name
    runs_dir = base_path / "runs"
    if runs_dir.exists():
        run_dirs = sorted([d for d in runs_dir.iterdir() if d.is_dir()], reverse=True)
        if run_dirs:
            return run_dirs[0]
    
    return None
