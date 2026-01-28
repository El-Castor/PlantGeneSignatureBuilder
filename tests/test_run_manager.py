"""
Unit tests for PGSB run manager
"""

import pytest
import tempfile
import shutil
from pathlib import Path
import yaml

from pgsb.io.run_manager import (
    compute_config_hash,
    create_run_dir,
    write_config_snapshot,
    update_latest_symlink,
    create_manifest
)


def test_config_hash_stable():
    """Config hash should be deterministic"""
    config1 = {'a': 1, 'b': 2, 'c': 3}
    config2 = {'c': 3, 'a': 1, 'b': 2}  # Different order
    
    hash1 = compute_config_hash(config1)
    hash2 = compute_config_hash(config2)
    
    assert hash1 == hash2, "Hash should be order-independent"
    assert len(hash1) == 8, "Hash should be 8 characters"


def test_config_hash_changes():
    """Config hash should change when content changes"""
    config1 = {'a': 1, 'b': 2}
    config2 = {'a': 1, 'b': 3}  # Different value
    
    hash1 = compute_config_hash(config1)
    hash2 = compute_config_hash(config2)
    
    assert hash1 != hash2, "Hash should differ for different configs"


def test_create_run_dir():
    """Test run directory creation"""
    with tempfile.TemporaryDirectory() as tmpdir:
        config = {'output_prefix': 'test', 'data': 'value'}
        
        run_info = create_run_dir(
            config, 
            'test_prefix',
            run_name='unittest',
            base_dir=tmpdir
        )
        
        # Check structure
        assert run_info['run_dir'].exists()
        assert run_info['outputs'].exists()
        assert run_info['qc'].exists()
        assert run_info['logs'].exists()
        
        # Check run_id format
        assert 'test_prefix' in run_info['run_id']
        assert 'unittest' in run_info['run_id']
        assert '__' in run_info['run_id']


def test_run_dir_collision():
    """Test that duplicate run_dir raises error without overwrite"""
    with tempfile.TemporaryDirectory() as tmpdir:
        config = {'output_prefix': 'test'}
        
        # Create first run
        run_info1 = create_run_dir(config, 'test', base_dir=tmpdir)
        
        # Try to create same run without overwrite (impossible due to timestamp)
        # But if we manually create the dir...
        (Path(tmpdir) / "runs" / "2099-01-01_000000__test__abcd1234").mkdir(parents=True)
        
        # This should work (different timestamp)
        run_info2 = create_run_dir(config, 'test', base_dir=tmpdir)
        
        # But trying to recreate exact same dir should fail
        with pytest.raises(FileExistsError):
            create_run_dir(config, 'test', base_dir=tmpdir, overwrite=False)


def test_write_config_snapshot():
    """Test config snapshot writing"""
    with tempfile.TemporaryDirectory() as tmpdir:
        config = {'key': 'value', 'number': 42}
        
        config_path = write_config_snapshot(config, tmpdir)
        
        assert config_path.exists()
        
        # Read back and verify
        with open(config_path, 'r') as f:
            loaded = yaml.safe_load(f)
        
        assert loaded == config


def test_create_manifest():
    """Test manifest creation with checksums"""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create some test files
        test_file1 = Path(tmpdir) / "test1.txt"
        test_file1.write_text("content1")
        
        test_file2 = Path(tmpdir) / "test2.txt"
        test_file2.write_text("content2")
        
        files_info = [
            {'path': 'test1.txt', 'description': 'File 1', 'rows': 10},
            {'path': 'test2.txt', 'description': 'File 2'}
        ]
        
        # Ensure outputs dir exists
        outputs_dir = Path(tmpdir) / "outputs"
        outputs_dir.mkdir()
        
        manifest_path = create_manifest(files_info, tmpdir)
        
        assert manifest_path.exists()
        
        # Verify manifest content
        import json
        with open(manifest_path, 'r') as f:
            manifest = json.load(f)
        
        assert 'run_id' in manifest
        assert 'created' in manifest
        assert len(manifest['files']) == 2
        assert manifest['files'][0]['sha256']  # Has checksum
