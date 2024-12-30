"""
Utility functions for MHC-peptide prediction pipeline.

This module provides utilities for:
- File handling and locking
- Resource optimization
- Batch processing
- Progress tracking
- Error handling
"""

import os
import logging
import tempfile
import shutil
import math
from typing import List, Dict, Optional, Any, Tuple, Generator
from pathlib import Path
import pandas as pd
from filelock import FileLock
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import psutil
from tqdm.auto import tqdm

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class ResourceManager:
    """Manages system resources for optimal performance."""
    
    @staticmethod
    def get_system_resources() -> Dict[str, int]:
        """Get available system resources."""
        return {
            'cpu_count': os.cpu_count() or 1,
            'memory_gb': psutil.virtual_memory().total // (1024**3),
            'memory_available_gb': psutil.virtual_memory().available // (1024**3)
        }

    @staticmethod
    def optimize_batch_parameters(
        total_items: int,
        memory_per_item_mb: float = 10,
        min_batch_size: int = 1000,
        max_batch_size: int = 10000
    ) -> Tuple[int, int]:
        """
        Calculate optimal batch size and number of processes.
        
        Args:
            total_items: Total number of items to process
            memory_per_item_mb: Estimated memory per item in MB
            min_batch_size: Minimum batch size
            max_batch_size: Maximum batch size
            
        Returns:
            Tuple of (batch_size, num_processes)
        """
        resources = ResourceManager.get_system_resources()
        
        # Reserve 20% of memory and at least 2 cores for system
        available_memory_mb = (resources['memory_available_gb'] * 1024) * 0.8
        available_cores = max(1, resources['cpu_count'] - 2)
        
        # Calculate maximum items that can fit in memory
        max_items_in_memory = int(available_memory_mb / memory_per_item_mb)
        
        # Calculate optimal batch size
        optimal_batches = available_cores * 4  # Aim for 4 batches per core
        batch_size = max(
            min_batch_size,
            min(
                max_batch_size,
                math.ceil(total_items / optimal_batches)
            )
        )
        
        # Calculate number of processes based on memory constraints
        max_processes_memory = max_items_in_memory // batch_size
        num_processes = min(available_cores, max_processes_memory)
        
        return batch_size, num_processes

class BatchProcessor:
    """Handles batch processing with progress tracking."""
    
    @staticmethod
    def create_batches(
        items: List[Any],
        batch_size: int
    ) -> Generator[List[Any], None, None]:
        """
        Split items into batches.
        
        Args:
            items: List of items to batch
            batch_size: Size of each batch
            
        Yields:
            Batches of items
        """
        for i in range(0, len(items), batch_size):
            yield items[i:i + batch_size]

    @staticmethod
    def process_batches(
        items: List[Any],
        process_func: callable,
        batch_size: int,
        n_jobs: int,
        desc: str = "Processing",
        use_processes: bool = True
    ) -> List[Any]:
        """
        Process items in batches with parallel execution.
        
        Args:
            items: Items to process
            process_func: Function to apply to each batch
            batch_size: Size of each batch
            n_jobs: Number of parallel jobs
            desc: Description for progress bar
            use_processes: Whether to use processes instead of threads
            
        Returns:
            List of processed results
        """
        batches = list(BatchProcessor.create_batches(items, batch_size))
        results = []
        
        executor_class = ProcessPoolExecutor if use_processes else ThreadPoolExecutor
        
        with executor_class(max_workers=n_jobs) as executor:
            futures = [
                executor.submit(process_func, batch)
                for batch in batches
            ]
            
            for future in tqdm(
                futures,
                total=len(batches),
                desc=desc,
                unit='batch'
            ):
                try:
                    result = future.result()
                    if result is not None:
                        results.extend(result)
                except Exception as e:
                    logger.error(f"Error processing batch: {str(e)}")
        
        return results

class FileManager:
    """Handles file operations with safety measures."""
    
    def __init__(self, base_path: str):
        """
        Initialize FileManager.
        
        Args:
            base_path: Base directory for file operations
        """
        self.base_path = Path(base_path)
        self.lock_dir = self.base_path / "locks"
        self.tmp_dir = self.base_path / "tmp"
        
        # Create directories
        self.lock_dir.mkdir(parents=True, exist_ok=True)
        self.tmp_dir.mkdir(parents=True, exist_ok=True)

    def safe_save(
        self,
        data: pd.DataFrame,
        filename: str,
        mode: str = 'csv'
    ) -> bool:
        """
        Safely save data to file with locking.
        
        Args:
            data: DataFrame to save
            filename: Target filename
            mode: Save mode ('csv' or 'parquet')
            
        Returns:
            bool indicating success
        """
        file_path = self.base_path / filename
        lock_path = self.lock_dir / f"{filename}.lock"
        
        with FileLock(str(lock_path)):
            try:
                # Create temporary file
                with tempfile.NamedTemporaryFile(
                    dir=str(self.tmp_dir),
                    delete=False
                ) as tmp_file:
                    tmp_path = Path(tmp_file.name)
                
                # Save to temporary file
                if mode == 'csv':
                    data.to_csv(tmp_path, index=False)
                elif mode == 'parquet':
                    data.to_parquet(tmp_path, index=False)
                else:
                    raise ValueError(f"Unsupported save mode: {mode}")
                
                # Atomic move to final location
                shutil.move(str(tmp_path), str(file_path))
                return True
                
            except Exception as e:
                logger.error(f"Error saving file {filename}: {str(e)}")
                return False
            
            finally:
                # Cleanup temporary file if it exists
                if 'tmp_path' in locals():
                    try:
                        tmp_path.unlink(missing_ok=True)
                    except:
                        pass

    def safe_load(
        self,
        filename: str,
        mode: str = 'csv'
    ) -> Optional[pd.DataFrame]:
        """
        Safely load data from file with locking.
        
        Args:
            filename: File to load
            mode: Load mode ('csv' or 'parquet')
            
        Returns:
            Loaded DataFrame or None if loading fails
        """
        file_path = self.base_path / filename
        lock_path = self.lock_dir / f"{filename}.lock"
        
        if not file_path.exists():
            return None
        
        with FileLock(str(lock_path)):
            try:
                if mode == 'csv':
                    return pd.read_csv(file_path)
                elif mode == 'parquet':
                    return pd.read_parquet(file_path)
                else:
                    raise ValueError(f"Unsupported load mode: {mode}")
            except Exception as e:
                logger.error(f"Error loading file {filename}: {str(e)}")
                return None

class ProgressTracker:
    """Tracks and reports progress of long-running operations."""
    
    def __init__(self, total: int, desc: str = "Progress"):
        """
        Initialize progress tracker.
        
        Args:
            total: Total number of items
            desc: Description for progress bar
        """
        self.total = total
        self.desc = desc
        self.pbar = tqdm(total=total, desc=desc)
        self.completed = 0
        self.failed = 0
        self.start_time = pd.Timestamp.now()

    def update(self, success: bool = True):
        """Update progress."""
        self.completed += 1
        if not success:
            self.failed += 1
        self.pbar.update(1)

    def get_stats(self) -> Dict[str, Any]:
        """Get progress statistics."""
        elapsed = pd.Timestamp.now() - self.start_time
        return {
            'total': self.total,
            'completed': self.completed,
            'failed': self.failed,
            'success_rate': (self.completed - self.failed) / self.total * 100,
            'elapsed': str(elapsed),
            'remaining': str(elapsed / self.completed * (self.total - self.completed))
            if self.completed > 0 else 'unknown'
        }

    def close(self):
        """Close progress bar and print final statistics."""
        self.pbar.close()
        stats = self.get_stats()
        logger.info(
            f"Completed {stats['completed']}/{stats['total']} items "
            f"({stats['success_rate']:.1f}% success) in {stats['elapsed']}"
        )
        if stats['failed'] > 0:
            logger.warning(f"Failed items: {stats['failed']}")