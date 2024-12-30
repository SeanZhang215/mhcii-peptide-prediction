"""
Data processing utilities for MHC-peptide prediction pipeline.

This module provides classes and functions for:
- Peptide sequence transformations
- Sample data processing
- Allele name mapping
- Data validation and preparation
"""

import pandas as pd
import logging
from typing import List, Dict, Set, Optional, Tuple, Any
import itertools
from dataclasses import dataclass
import re
from pathlib import Path

from .utils import FileManager, BatchProcessor, ProgressTracker, ResourceManager

logger = logging.getLogger(__name__)

@dataclass
class ProcessedSample:
    """Container for processed sample data."""
    peptides: List[str]
    alleles: List[str]
    raw_counts: int
    processed_df: Optional[pd.DataFrame] = None

    def validate(self) -> Dict[str, Any]:
        """Validate sample data and return statistics."""
        return {
            'total_peptides': len(self.peptides),
            'mean_length': sum(len(p) for p in self.peptides) / len(self.peptides),
            'num_variants': len(self.processed_df) if self.processed_df is not None else 0,
            'num_alleles': len(self.alleles),
            'has_processed_df': self.processed_df is not None
        }

class PeptideProcessor:
    """Handles peptide sequence processing and transformations."""
    
    @staticmethod
    def flip_first_two(peptide: str) -> str:
        """Flip first two residues of peptide."""
        if len(peptide) < 2:
            return peptide
        return peptide[1] + peptide[0] + peptide[2:]

    @staticmethod
    def flip_last_two(peptide: str) -> str:
        """Flip last two residues of peptide."""
        if len(peptide) < 2:
            return peptide
        return peptide[:-2] + peptide[-1] + peptide[-2]

    @staticmethod
    def flip_both(peptide: str) -> str:
        """Flip both first and last two residues."""
        return PeptideProcessor.flip_last_two(PeptideProcessor.flip_first_two(peptide))
    
    @staticmethod
    def invert_peptide(peptide: str) -> str:
        """Reverse the peptide sequence."""
        return peptide[::-1]

    def process_single_peptide(self, peptide: str) -> Dict[str, Any]:
        """Process a single peptide to create variants."""
        variants = []
        for inversion in [0, 1]:
            base_peptide = self.invert_peptide(peptide) if inversion else peptide
            for flip_type in ['raw', 'first_two', 'last_two', 'both']:
                if flip_type == 'raw':
                    new_peptide = base_peptide
                elif flip_type == 'first_two':
                    new_peptide = self.flip_first_two(base_peptide)
                elif flip_type == 'last_two':
                    new_peptide = self.flip_last_two(base_peptide)
                else:  # 'both'
                    new_peptide = self.flip_both(base_peptide)
                
                variants.append({
                    'raw_peptides': peptide,
                    'inverted_manual': inversion,
                    'flipped': flip_type,
                    'new_peptides': new_peptide
                })
        return variants

    def create_peptide_variants(self, peptides: List[str]) -> pd.DataFrame:
        """
        Create all sequence variants for a list of peptides using batch processing.
        
        Args:
            peptides: List of peptide sequences
        
        Returns:
            DataFrame containing all variants with transformations applied
        """
        # Use BatchProcessor for efficient processing
        batch_size, n_jobs = ResourceManager.optimize_batch_parameters(
            total_items=len(peptides),
            memory_per_item_mb=1  
        )
        
        # Process peptides in batches
        variants = BatchProcessor.process_batches(
            items=peptides,
            process_func=lambda batch: [
                variant 
                for peptide in batch 
                for variant in self.process_single_peptide(peptide)
            ],
            batch_size=batch_size,
            n_jobs=n_jobs,
            desc="Creating peptide variants"
        )
        
        return pd.DataFrame(variants)

class AlleleMapper:
    """Handles MHC allele name mapping and validation."""
    
    def __init__(self, allele_list_path: str):
        """
        Initialize AlleleMapper with netMHCIIpan allele list.
        
        Args:
            allele_list_path: Path to netMHCIIpan allele.list file
        """
        self.file_manager = FileManager(str(Path(allele_list_path).parent))
        self.allele_list = self._load_allele_list(allele_list_path)
        
    def _load_allele_list(self, path: str) -> pd.DataFrame:
        """Load and validate allele list file using FileManager."""
        try:
            filename = Path(path).name
            df = self.file_manager.safe_load(filename)
            if df is None:
                raise FileNotFoundError(f"Could not load allele list from {path}")
            
            df.columns = ['Allele']
            return df
            
        except Exception as e:
            logger.error(f"Error loading allele list from {path}: {e}")
            raise
    
    def find_matching_alleles(self, allele: str) -> List[str]:
        """Find matching alleles in the netMHCIIpan allele list."""
        return self.allele_list[
            self.allele_list['Allele'].str.contains(allele, regex=False)
        ]['Allele'].tolist()
    
    def process_allele_groups(
        self,
        dq_alleles: Dict[str, List[str]],
        dp_alleles: Dict[str, List[str]],
        dr_alleles: List[str]
    ) -> Tuple[List[str], List[str]]:
        """Process allele groups to find valid combinations."""
        final_alleles = []
        unmatched_alleles = []
        
        # Process DR alleles
        for dr_allele in dr_alleles:
            matching = self.find_matching_alleles(dr_allele)
            if matching:
                final_alleles.extend(matching)
            else:
                unmatched_alleles.append(dr_allele)
        
        # Process DQ combinations
        if dq_alleles['A'] and dq_alleles['B']:
            for dqa, dqb in itertools.product(dq_alleles['A'], dq_alleles['B']):
                combined = f"HLA-{dqa}-{dqb}"
                matching = self.find_matching_alleles(combined)
                if matching:
                    final_alleles.extend(matching)
                else:
                    unmatched_alleles.append(combined)
        
        # Process DP combinations
        if dp_alleles['A'] and dp_alleles['B']:
            for dpa, dpb in itertools.product(dp_alleles['A'], dp_alleles['B']):
                combined = f"HLA-{dpa}-{dpb}"
                matching = self.find_matching_alleles(combined)
                if matching:
                    final_alleles.extend(matching)
                else:
                    unmatched_alleles.append(combined)
        
        return list(set(final_alleles)), list(set(unmatched_alleles))
    
    def process_hla_typing(self, hla_row: pd.Series) -> Tuple[List[str], List[str]]:
        """
        Process HLA typing data to get valid alleles.
        
        Args:
            hla_row: Row from HLA typing DataFrame
            
        Returns:
            Tuple of (valid_alleles, unmatched_alleles)
        """
        dq_alleles = {'A': [], 'B': []}
        dp_alleles = {'A': [], 'B': []}
        dr_alleles = []
        
        for col in hla_row.index[7:]:
            allele = hla_row[col]
            if pd.isna(allele):
                continue
                
            allele = (allele.replace('HLA-', '')
                     .replace('*', '_')
                     .replace(':', ''))
            
            # Sort alleles into respective groups
            if allele.startswith('DQA'):
                allele = allele.replace('_', '')
                dq_alleles['A'].append(allele)
            elif allele.startswith('DQB'):
                allele = allele.replace('_', '')
                dq_alleles['B'].append(allele)
            elif allele.startswith('DPA'):
                allele = allele.replace('_', '')
                dp_alleles['A'].append(allele)
            elif allele.startswith('DPB'):
                allele = allele.replace('_', '')
                dp_alleles['B'].append(allele)
            elif allele.startswith('DR'):
                dr_alleles.append(allele)
        
        return self.process_allele_groups(dq_alleles, dp_alleles, dr_alleles)

class SampleDataProcessor:
    """Processes and prepares sample data for prediction."""
    
    def __init__(self, min_peptide_length: int = 13):
        """
        Initialize processor with configuration.
        
        Args:
            min_peptide_length: Minimum peptide length to include
        """
        self.min_length = min_peptide_length
        self.peptide_processor = PeptideProcessor()
    
    def process_ms_data(self, df: pd.DataFrame) -> Dict[str, ProcessedSample]:
        """
        Process MS data and prepare samples for prediction.
        
        Args:
            df: DataFrame containing MS data
            
        Returns:
            Dictionary mapping sample IDs to processed data
        """
        processed_samples = {}
        progress = ProgressTracker(
            total=len(df['sampleId'].unique()),
            desc="Processing samples"
        )
        
        # Group by sample ID and get unique peptides
        sample_groups = df.groupby('sampleId')
        for sample_id, group in sample_groups:
            try:
                # Filter peptides by length
                peptides = [
                    peptide for peptide in group['cells.peptide'].unique()
                    if len(peptide) >= self.min_length
                ]
                
                if not peptides:
                    logger.warning(f"No valid peptides found for sample {sample_id}")
                    progress.update(False)
                    continue
                
                # Create variants for peptides
                processed_df = self.peptide_processor.create_peptide_variants(peptides)
                
                processed_samples[sample_id] = ProcessedSample(
                    peptides=peptides,
                    alleles=[], 
                    raw_counts=len(peptides),
                    processed_df=processed_df
                )
                progress.update(True)
                
            except Exception as e:
                logger.error(f"Error processing sample {sample_id}: {e}")
                progress.update(False)
        
        progress.close()
        return processed_samples
    
    def validate_processed_data(
        self,
        processed_samples: Dict[str, ProcessedSample]
    ) -> Dict[str, Dict[str, Any]]:
        """Validate processed sample data using sample validation."""
        return {
            sample_id: data.validate()
            for sample_id, data in processed_samples.items()
        }