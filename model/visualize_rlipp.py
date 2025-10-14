#!/usr/bin/env python3
"""
Visualize RLIPP scores from calculate_rlipp.py output

Creates plots to help interpret RLIPP results:
1. Distribution of RLIPP scores (histogram)
2. Top terms by RLIPP (bar plot)
3. RLIPP vs layer position (scatter plot)
4. Predictive power comparison (term vs children)
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_rlipp_distribution(df, output_prefix):
    """Plot histogram of RLIPP scores"""
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    ax.hist(df['rlipp'], bins=50, edgecolor='black', alpha=0.7)
    ax.axvline(x=0, color='red', linestyle='--', linewidth=2, label='Zero RLIPP')
    ax.set_xlabel('RLIPP Score', fontsize=12)
    ax.set_ylabel('Number of Terms', fontsize=12)
    ax.set_title('Distribution of RLIPP Scores', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Add statistics text
    stats_text = f"Mean: {df['rlipp'].mean():.3f}\n"
    stats_text += f"Median: {df['rlipp'].median():.3f}\n"
    stats_text += f"Positive: {(df['rlipp'] > 0).sum()}\n"
    stats_text += f"Negative: {(df['rlipp'] < 0).sum()}"
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_distribution.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_prefix}_distribution.png")
    plt.close()


def plot_top_terms(df, output_prefix, top_n=20):
    """Plot top N terms by RLIPP score"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    # Top positive RLIPP
    top_positive = df.nlargest(top_n, 'rlipp')
    colors_pos = ['green' if x > 0 else 'red' for x in top_positive['rlipp']]
    ax1.barh(range(len(top_positive)), top_positive['rlipp'], color=colors_pos, alpha=0.7)
    ax1.set_yticks(range(len(top_positive)))
    ax1.set_yticklabels(top_positive['term'], fontsize=8)
    ax1.set_xlabel('RLIPP Score', fontsize=12)
    ax1.set_title(f'Top {top_n} Terms by RLIPP Score (Positive)', fontsize=14, fontweight='bold')
    ax1.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    ax1.grid(True, alpha=0.3, axis='x')
    ax1.invert_yaxis()

    # Bottom (most negative) RLIPP
    bottom_negative = df.nsmallest(top_n, 'rlipp')
    colors_neg = ['green' if x > 0 else 'red' for x in bottom_negative['rlipp']]
    ax2.barh(range(len(bottom_negative)), bottom_negative['rlipp'], color=colors_neg, alpha=0.7)
    ax2.set_yticks(range(len(bottom_negative)))
    ax2.set_yticklabels(bottom_negative['term'], fontsize=8)
    ax2.set_xlabel('RLIPP Score', fontsize=12)
    ax2.set_title(f'Bottom {top_n} Terms by RLIPP Score (Negative)', fontsize=14, fontweight='bold')
    ax2.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    ax2.grid(True, alpha=0.3, axis='x')
    ax2.invert_yaxis()

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_top_terms.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_prefix}_top_terms.png")
    plt.close()


def plot_rlipp_by_layer(df, output_prefix):
    """Plot RLIPP scores by layer position"""
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    # Scatter plot with color coding
    scatter = ax.scatter(df['layer'], df['rlipp'],
                        c=df['rlipp'], cmap='RdYlGn',
                        alpha=0.6, s=50, edgecolors='black', linewidth=0.5)

    ax.axhline(y=0, color='black', linestyle='--', linewidth=1)
    ax.set_xlabel('Layer (0 = leaves, higher = closer to root)', fontsize=12)
    ax.set_ylabel('RLIPP Score', fontsize=12)
    ax.set_title('RLIPP Score by Hierarchy Layer', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('RLIPP Score', fontsize=10)

    # Add layer statistics
    layer_stats = df.groupby('layer')['rlipp'].agg(['mean', 'std', 'count'])
    for layer, row in layer_stats.iterrows():
        ax.plot(layer, row['mean'], 'ro', markersize=10, markeredgecolor='black',
                markeredgewidth=1.5, label='Mean' if layer == layer_stats.index[0] else '')

    ax.legend()
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_by_layer.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_prefix}_by_layer.png")
    plt.close()


def plot_power_comparison(df, output_prefix):
    """Plot term power vs children power"""
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    # Scatter plot
    colors = ['green' if x > 0 else 'red' for x in df['rlipp']]
    ax.scatter(df['children_power'], df['term_power'],
              c=df['rlipp'], cmap='RdYlGn', alpha=0.6, s=50,
              edgecolors='black', linewidth=0.5)

    # Add diagonal line (equal power)
    min_val = min(df['children_power'].min(), df['term_power'].min())
    max_val = max(df['children_power'].max(), df['term_power'].max())
    ax.plot([min_val, max_val], [min_val, max_val],
            'k--', linewidth=2, label='Equal Power', alpha=0.5)

    ax.set_xlabel('Children Predictive Power (Correlation)', fontsize=12)
    ax.set_ylabel('Term Predictive Power (Correlation)', fontsize=12)
    ax.set_title('Predictive Power: Term vs Children', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap='RdYlGn',
                                norm=plt.Normalize(vmin=df['rlipp'].min(),
                                                  vmax=df['rlipp'].max()))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('RLIPP Score', fontsize=10)

    # Add text indicating regions
    ax.text(0.05, 0.95, 'Term better\n(Positive RLIPP)',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))
    ax.text(0.95, 0.05, 'Children better\n(Negative RLIPP)',
            transform=ax.transAxes, fontsize=10, verticalalignment='bottom',
            horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.5))

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_power_comparison.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_prefix}_power_comparison.png")
    plt.close()


def plot_children_analysis(df, output_prefix):
    """Analyze RLIPP by number of children"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Box plot by number of children (grouped)
    df['children_group'] = pd.cut(df['num_children'], bins=[0, 2, 5, 10, 20, 100],
                                   labels=['1-2', '3-5', '6-10', '11-20', '20+'])

    df.boxplot(column='rlipp', by='children_group', ax=ax1)
    ax1.set_xlabel('Number of Children', fontsize=12)
    ax1.set_ylabel('RLIPP Score', fontsize=12)
    ax1.set_title('RLIPP Distribution by Number of Children', fontsize=12)
    ax1.axhline(y=0, color='red', linestyle='--', linewidth=1)
    plt.sca(ax1)
    plt.xticks(rotation=0)

    # Scatter: number of children vs RLIPP
    ax2.scatter(df['num_children'], df['rlipp'], alpha=0.5, s=30)
    ax2.set_xlabel('Number of Children', fontsize=12)
    ax2.set_ylabel('RLIPP Score', fontsize=12)
    ax2.set_title('RLIPP vs Number of Children', fontsize=12)
    ax2.axhline(y=0, color='red', linestyle='--', linewidth=1)
    ax2.grid(True, alpha=0.3)

    plt.suptitle('')  # Remove automatic title
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_children_analysis.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_prefix}_children_analysis.png")
    plt.close()


def generate_summary_report(df, output_file):
    """Generate text summary report"""
    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("RLIPP Score Analysis Summary Report\n")
        f.write("="*80 + "\n\n")

        f.write(f"Total number of terms: {len(df)}\n\n")

        f.write("RLIPP Score Statistics:\n")
        f.write("-"*40 + "\n")
        f.write(f"Mean:     {df['rlipp'].mean():>10.4f}\n")
        f.write(f"Median:   {df['rlipp'].median():>10.4f}\n")
        f.write(f"Std Dev:  {df['rlipp'].std():>10.4f}\n")
        f.write(f"Min:      {df['rlipp'].min():>10.4f}\n")
        f.write(f"Max:      {df['rlipp'].max():>10.4f}\n\n")

        f.write("RLIPP Distribution:\n")
        f.write("-"*40 + "\n")
        f.write(f"Positive RLIPP (>0):     {(df['rlipp'] > 0).sum():>6} ({100*(df['rlipp'] > 0).sum()/len(df):.1f}%)\n")
        f.write(f"Negative RLIPP (<0):     {(df['rlipp'] < 0).sum():>6} ({100*(df['rlipp'] < 0).sum()/len(df):.1f}%)\n")
        f.write(f"Near-zero RLIPP (±0.01): {(df['rlipp'].abs() < 0.01).sum():>6} ({100*(df['rlipp'].abs() < 0.01).sum()/len(df):.1f}%)\n")
        f.write(f"High RLIPP (>0.3):       {(df['rlipp'] > 0.3).sum():>6} ({100*(df['rlipp'] > 0.3).sum()/len(df):.1f}%)\n\n")

        f.write("Top 10 Terms by RLIPP Score:\n")
        f.write("-"*40 + "\n")
        for idx, row in df.nlargest(10, 'rlipp').iterrows():
            f.write(f"{row['term']:<30} {row['rlipp']:>8.4f} (Layer {row['layer']}, {row['num_children']} children)\n")

        f.write("\n" + "Bottom 10 Terms by RLIPP Score:\n")
        f.write("-"*40 + "\n")
        for idx, row in df.nsmallest(10, 'rlipp').iterrows():
            f.write(f"{row['term']:<30} {row['rlipp']:>8.4f} (Layer {row['layer']}, {row['num_children']} children)\n")

        f.write("\n" + "="*80 + "\n")

    print(f"Saved: {output_file}")


def main():
    parser = argparse.ArgumentParser(description='Visualize RLIPP scores')
    parser.add_argument('-input', help='Input RLIPP scores TSV file', type=str, required=True)
    parser.add_argument('-output', help='Output prefix for plots', type=str, default='rlipp_viz')
    parser.add_argument('-top_n', help='Number of top terms to show', type=int, default=20)

    args = parser.parse_args()

    # Read data
    print(f"Reading RLIPP scores from: {args.input}")
    df = pd.read_csv(args.input, sep='\t')
    print(f"Loaded {len(df)} terms\n")

    # Create visualizations
    print("Generating visualizations...")
    plot_rlipp_distribution(df, args.output)
    plot_top_terms(df, args.output, args.top_n)
    plot_rlipp_by_layer(df, args.output)
    plot_power_comparison(df, args.output)
    plot_children_analysis(df, args.output)

    # Generate text report
    print("\nGenerating summary report...")
    generate_summary_report(df, f'{args.output}_report.txt')

    print("\n" + "="*60)
    print("Visualization complete!")
    print(f"Created {5} plots and 1 summary report")
    print("="*60)


if __name__ == "__main__":
    main()
