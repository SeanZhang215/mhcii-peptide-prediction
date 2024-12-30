"""
NetMHCII Peptide Binding Predictor

This module provides a Python interface for running netMHCIIpan predictions.
It handles peptide processing, batch execution, and results management.
"""

import os
import pandas as pd
import subprocess
import tempfile
from typing import List, Dict, Optional, Tuple
import logging
from pathlib import Path

from .utils import (
    ResourceManager,
    BatchProcessor,
    FileManager,
    ProgressTracker
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class NetMHCIIPredictor:
    """
    Manages NetMHCIIpan predictions for peptide-MHC Class II binding.
    
    Attributes:
        base_path: Base directory for storing results
        model_path: Path to NetMHCIIpan installation
        file_manager: Handles file operations safely
    """
    
    def __init__(
        self, 
        base_path: str, 
        model_path: str
    ):
        """Initialize the predictor."""
        self.base_path = Path(base_path)
        self.model_path = Path(model_path)
        self.file_manager = FileManager(base_path)
        
        self._setup_environment()

    def _setup_environment(self):
        """Configure environment variables for NetMHCIIpan."""
        os.environ.update({
            'NMHOME': str(self.model_path),
            'TMPDIR': str(self.base_path / 'tmp'),
            'UNIX': os.uname().sysname,
            'AR': os.uname().machine,
        })
        
        platform = f"{os.environ['UNIX']}_{os.environ['AR']}"
        os.environ.update({
            'NETMHCIIpan': os.environ['NMHOME'],
            'NetMHCIIpanPLAT': f"{os.environ['NMHOME']}/{platform}"
        })
        
        os.makedirs(os.environ['TMPDIR'], exist_ok=True)

    def run_prediction(
        self, 
        peptides: List[str], 
        alleles: List[str], 
        context: bool = False
    ) -> Optional[pd.DataFrame]:
        """
        Run NetMHCIIpan prediction for peptides and alleles.
        
        Args:
            peptides: List of peptide sequences
            alleles: List of MHC alleles
            context: Whether to include context information
            
        Returns:
            DataFrame with prediction results or None if prediction fails
        """
        try:
            with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_input:
                for peptide in peptides:
                    temp_input.write(f"{peptide}\n")
                input_file = temp_input.name

            cmd = [
                f"{os.environ['NetMHCIIpanPLAT']}/bin/NetMHCIIpan-4.3",
                '-BA', '1',
                '-f', input_file,
                '-inptype', '1',
                '-inv_all', '1',
                '-a', ','.join(alleles)
            ]
            
            if context:
                cmd.append('-context')

            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            return self._parse_output(result.stdout)

        except subprocess.CalledProcessError as e:
            logger.error(f"Error running NetMHCIIpan: {e}")
            logger.error(f"Error output: {e.stderr}")
            return None
            
        except Exception as e:
            logger.error(f"Unexpected error in prediction: {e}")
            return None
            
        finally:
            if 'input_file' in locals():
                os.remove(input_file)

    def _parse_output(self, output: str) -> pd.DataFrame:
        """Parse NetMHCIIpan output into a DataFrame."""
        lines = output.split('\n')
        data = []
        current_allele = None
        
        for line in lines:
            if line.startswith('# Allele:'):
                current_allele = line.split(':')[1].strip()
                
            elif line.strip() and not line.startswith(('#', '-')):
                parts = line.split()
                if len(parts) >= 12 and parts[0].isdigit():
                    data.append({
                        'Allele': current_allele,
                        'Pos': int(parts[0]),
                        'MHC': parts[1],
                        'Peptide': parts[2],
                        'Of': int(parts[3]),
                        'Core': parts[4],
                        'Core_Rel': float(parts[5]),
                        'Inverted': int(parts[6]),
                        'Identity': parts[7],
                        'Score_EL': float(parts[8]),
                        '%Rank_EL': float(parts[9]),
                        'Exp_Bind': parts[10] if parts[10] != 'NA' else None,
                        'Score_BA': float(parts[11]),
                        '%Rank_BA': float(parts[12]),
                        'Affinity(nM)': float(parts[13]),
                        'BindLevel': parts[15] if len(parts) > 15 else None
                    })
        
        return pd.DataFrame(data)

    def process_batch(
        self,
        peptides: List[str],
        alleles: List[str],
        context: bool = False
    ) -> Optional[pd.DataFrame]:
        """Process a single batch of peptides."""
        return self.run_prediction(peptides, alleles, context)

    def process_samples(
        self,
        sample_data: Dict[str, Dict],
        min_batch_size: int = 1000,
        max_batch_size: int = 10000
    ) -> Dict[str, pd.DataFrame]:
        """
        Process multiple samples with optimized batch processing.
        
        Args:
            sample_data: Dictionary of sample information
            min_batch_size: Minimum batch size
            max_batch_size: Maximum batch size
            
        Returns:
            Dictionary mapping sample IDs to results DataFrames
        """
        total_peptides = sum(len(data['peptides']) for data in sample_data.values())
        
        # Use ResourceManager to optimize parameters
        batch_size, n_processes = ResourceManager.optimize_batch_parameters(
            total_items=total_peptides,
            memory_per_item_mb=20,  # Estimated memory per peptide
            min_batch_size=min_batch_size,
            max_batch_size=max_batch_size
        )
        
        # Initialize progress tracking
        progress = ProgressTracker(
            total=len(sample_data),
            desc="Processing samples"
        )
        
        results = {}
        for sample_id, data in sample_data.items():
            try:
                # Process peptides in batches
                sample_results = BatchProcessor.process_batches(
                    items=data['peptides'],
                    process_func=lambda batch: self.process_batch(
                        batch, data['alleles']
                    ),
                    batch_size=batch_size,
                    n_jobs=n_processes,
                    desc=f"Sample {sample_id}"
                )
                
                if sample_results:
                    # Combine results and add sample ID
                    combined_df = pd.concat(sample_results, ignore_index=True)
                    combined_df['sample_id'] = sample_id
                    
                    # Save results
                    self.save_sample_results(combined_df, sample_id)
                    results[sample_id] = combined_df
                    progress.update(True)
                else:
                    logger.warning(f"No results for sample {sample_id}")
                    progress.update(False)
                    
            except Exception as e:
                logger.error(f"Error processing sample {sample_id}: {e}")
                progress.update(False)
        
        progress.close()
        return results

    def save_sample_results(
        self,
        results: pd.DataFrame,
        sample_id: str,
        mode: str = 'csv'
    ) -> bool:
        """
        Save sample results using FileManager.
        
        Args:
            results: DataFrame containing prediction results
            sample_id: Sample identifier
            mode: File format mode
            
        Returns:
            bool indicating success
        """
        filename = f"{sample_id}_results.{mode}"
        return self.file_manager.safe_save(results, filename, mode)

    def load_sample_results(
        self,
        sample_id: str,
        mode: str = 'csv'
    ) -> Optional[pd.DataFrame]:
        """
        Load sample results using FileManager.
        
        Args:
            sample_id: Sample identifier
            mode: File format mode
            
        Returns:
            DataFrame containing results or None if loading fails
        """
        filename = f"{sample_id}_results.{mode}"
        return self.file_manager.safe_load(filename, mode)