"""
Visualization utilities for MHC peptide prediction analysis.

This module provides visualization functions for:
- Peptide-HLA binding distributions
- Sample comparisons
- Target analysis
- Gene expression analysis
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List, Dict, Optional, Tuple, Any
import numpy as np
from dataclasses import dataclass
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

@dataclass
class PlotConfig:
    """Configuration for plot styling."""
    figsize: Tuple[int, int] = (12, 8)
    palette: Optional[Dict[str, str]] = None
    style: str = "whitegrid"
    font_scale: float = 1.2
    context: str = "paper"

class BindingDistributionPlotter:
    """Handles visualization of peptide-HLA binding distributions."""

    def __init__(self, config: Optional[PlotConfig] = None):
        """Initialize plotter with configuration."""
        self.config = config or PlotConfig()
        self._setup_style()

    def _setup_style(self):
        """Set up plot styling."""
        sns.set_theme(
            style=self.config.style,
            font_scale=self.config.font_scale,
            context=self.config.context
        )
        if self.config.palette:
            sns.set_palette(self.config.palette)

    def plot_rank_distribution(
        self,
        data: pd.DataFrame,
        rank_col: str = '%Rank_EL',
        threshold: float = 5.0,
        title: Optional[str] = None
    ) -> plt.Figure:
        """
        Create violin plot of binding rank distribution.
        
        Args:
            data: DataFrame containing binding predictions
            rank_col: Name of rank column
            threshold: Binding threshold
            title: Optional plot title
            
        Returns:
            matplotlib Figure
        """
        fig, ax = plt.subplots(figsize=self.config.figsize)

        # Calculate statistics
        total_pairs = len(data)
        above_threshold = (data[rank_col] > threshold).sum() / total_pairs * 100
        below_threshold = (data[rank_col] <= threshold).sum() / total_pairs * 100

        # Create violin plot
        sns.violinplot(
            y=rank_col,
            data=data,
            ax=ax,
            color='skyblue',
            cut=0
        )

        # Add threshold line
        ax.axhline(
            threshold,
            color='red',
            linestyle='--',
            label=f'Threshold ({threshold})'
        )

        # Customize plot
        ax.set_yscale('log')
        ax.set_ylim(0.1, 100)
        ax.set_title(title or 'Binding Rank Distribution')
        ax.set_ylabel('Rank EL (%)')

        # Add statistics
        ax.text(
            0.02, 0.98,
            f'Above {threshold}%: {above_threshold:.1f}%\n'
            f'Below {threshold}%: {below_threshold:.1f}%',
            transform=ax.transAxes,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
        )

        plt.tight_layout()
        return fig

    def plot_comparison(
        self,
        data_dict: Dict[str, pd.DataFrame],
        rank_col: str = '%Rank_EL',
        threshold: float = 5.0,
        title: Optional[str] = None
    ) -> plt.Figure:
        """
        Create comparison violin plot for multiple datasets.
        
        Args:
            data_dict: Dictionary mapping names to DataFrames
            rank_col: Name of rank column
            threshold: Binding threshold
            title: Optional plot title
            
        Returns:
            matplotlib Figure
        """
        fig, ax = plt.subplots(figsize=self.config.figsize)

        # Prepare data for plotting
        plot_data = pd.concat([
            pd.DataFrame({
                'Group': name,
                'Rank': df[rank_col]
            })
            for name, df in data_dict.items()
        ])

        # Create violin plot
        sns.violinplot(
            x='Group',
            y='Rank',
            data=plot_data,
            ax=ax,
            cut=0
        )

        # Add threshold line
        ax.axhline(
            threshold,
            color='red',
            linestyle='--',
            label=f'Threshold ({threshold})'
        )

        # Customize plot
        ax.set_yscale('log')
        ax.set_ylim(0.1, 100)
        ax.set_title(title or 'Binding Rank Comparison')
        ax.set_ylabel('Rank EL (%)')

        # Add statistics for each group
        stats = []
        for name, df in data_dict.items():
            above = (df[rank_col] > threshold).mean() * 100
            below = (df[rank_col] <= threshold).mean() * 100
            stats.append(f'{name}:\nAbove: {above:.1f}%\nBelow: {below:.1f}%')

        ax.text(
            1.02, 0.98,
            '\n\n'.join(stats),
            transform=ax.transAxes,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
        )

        plt.tight_layout()
        return fig

class TargetAnalysisPlotter:
    """Handles visualization of target analysis results."""

    def __init__(self, config: Optional[PlotConfig] = None):
        """Initialize plotter with configuration."""
        self.config = config or PlotConfig()
        self._setup_style()

    def _setup_style(self):
        """Set up plot styling."""
        sns.set_theme(
            style=self.config.style,
            font_scale=self.config.font_scale,
            context=self.config.context
        )
        if self.config.palette:
            sns.set_palette(self.config.palette)

    def plot_target_overview(
        self,
        tumor_data: pd.DataFrame,
        nat_data: pd.DataFrame,
        genes_of_interest: List[str]
    ) -> Tuple[plt.Figure, plt.Figure, plt.Figure, plt.Figure]:
        """
        Create comprehensive target analysis plots.
        
        Args:
            tumor_data: Tumor sample data
            nat_data: NAT sample data
            genes_of_interest: List of target genes
            
        Returns:
            Tuple of figures (peptide_counts, ion_scores, peptide_lengths, expression)
        """
        # Create figures
        fig1 = self._plot_peptide_counts(tumor_data, nat_data, genes_of_interest)
        fig2 = self._plot_ion_scores(tumor_data, genes_of_interest)
        fig3 = self._plot_peptide_lengths(tumor_data, nat_data)
        fig4 = self._plot_gene_expression(tumor_data, nat_data, genes_of_interest)

        return fig1, fig2, fig3, fig4

    def _plot_peptide_counts(
        self,
        tumor_data: pd.DataFrame,
        nat_data: pd.DataFrame,
        genes: List[str]
    ) -> plt.Figure:
        """Plot peptide counts for target genes."""
        fig, ax = plt.subplots(figsize=self.config.figsize)

        # Calculate peptide counts
        tumor_counts = tumor_data[tumor_data['gene_name'].isin(genes)].groupby('gene_name')['Peptide'].nunique()
        nat_counts = nat_data[nat_data['gene_name'].isin(genes)].groupby('gene_name')['Peptide'].nunique()

        # Create bar plot
        x = np.arange(len(genes))
        width = 0.35

        ax.bar(x - width/2, [tumor_counts.get(g, 0) for g in genes], width, label='Tumor')
        ax.bar(x + width/2, [nat_counts.get(g, 0) for g in genes], width, label='NAT')

        ax.set_ylabel('Number of Unique Peptides')
        ax.set_title('Peptide Counts by Gene')
        ax.set_xticks(x)
        ax.set_xticklabels(genes, rotation=45)
        ax.legend()

        plt.tight_layout()
        return fig

    def _plot_ion_scores(self, tumor_data: pd.DataFrame, genes: List[str]) -> plt.Figure:
        """Plot ion score distributions for target genes."""
        fig, ax = plt.subplots(figsize=self.config.figsize)

        data = tumor_data[tumor_data['gene_name'].isin(genes)]
        sns.boxplot(
            x='gene_name',
            y='IonScore',
            data=data,
            ax=ax
        )

        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
        ax.set_title('Ion Score Distribution in Tumor Samples')
        ax.set_ylabel('Ion Score')

        plt.tight_layout()
        return fig

    def _plot_peptide_lengths(
        self,
        tumor_data: pd.DataFrame,
        nat_data: pd.DataFrame
    ) -> plt.Figure:
        """Plot peptide length distributions."""
        fig, ax = plt.subplots(figsize=self.config.figsize)

        tumor_data['length'] = tumor_data['Peptide'].str.len()
        nat_data['length'] = nat_data['Peptide'].str.len()

        sns.histplot(
            data=tumor_data,
            x='length',
            label='Tumor',
            alpha=0.5,
            ax=ax
        )
        sns.histplot(
            data=nat_data,
            x='length',
            label='NAT',
            alpha=0.5,
            ax=ax
        )

        ax.set_title('Peptide Length Distribution')
        ax.set_xlabel('Peptide Length')
        ax.set_ylabel('Count')
        ax.legend()

        plt.tight_layout()
        return fig

    def _plot_gene_expression(
        self,
        tumor_data: pd.DataFrame,
        nat_data: pd.DataFrame,
        genes: List[str]
    ) -> plt.Figure:
        """Plot gene expression comparison."""
        fig, ax = plt.subplots(figsize=self.config.figsize)

        data = pd.concat([
            tumor_data[tumor_data['gene_name'].isin(genes)].assign(group='Tumor'),
            nat_data[nat_data['gene_name'].isin(genes)].assign(group='NAT')
        ])

        sns.boxplot(
            x='gene_name',
            y='expression',
            hue='group',
            data=data,
            ax=ax
        )

        ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
        ax.set_title('Gene Expression Comparison')
        ax.set_ylabel('Expression Level (log2 TPM)')

        plt.tight_layout()
        return fig

def save_figures(
    figures: Dict[str, plt.Figure],
    output_dir: str,
    format: str = 'png',
    dpi: int = 300
) -> None:
    """
    Save multiple figures to files.
    
    Args:
        figures: Dictionary mapping names to figures
        output_dir: Directory to save figures
        format: File format for figures
        dpi: Resolution for saved figures
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    for name, fig in figures.items():
        try:
            filename = f"{name}.{format}"
            fig.savefig(
                output_path / filename,
                format=format,
                dpi=dpi,
                bbox_inches='tight'
            )
            logger.info(f"Saved figure {name} to {filename}")
        except Exception as e:
            logger.error(f"Error saving figure {name}: {e}")