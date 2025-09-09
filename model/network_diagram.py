#!/usr/bin/env python3
"""
Simple illustration of the DCellNN network architecture
Based on the output from test/out.txt
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch
import numpy as np

def create_network_diagram():
    fig, ax = plt.subplots(1, 1, figsize=(14, 10))
    
    # Colors
    gene_color = '#FF6B6B'      # Red for genes
    direct_color = '#4ECDC4'    # Teal for direct gene layers
    term_color = '#45B7D1'      # Blue for GO term processing
    dropout_color = '#96CEB4'   # Green for dropout layers
    root_color = '#FECA57'      # Yellow for root
    final_color = '#FF9FF3'     # Pink for final output
    
    # Gene Input Layer
    gene_box = FancyBboxPatch((1, 8), 2, 1, boxstyle="round,pad=0.1", 
                             facecolor=gene_color, edgecolor='black', linewidth=2)
    ax.add_patch(gene_box)
    ax.text(2, 8.5, 'Gene Input\n(2 genes)', ha='center', va='center', fontweight='bold', fontsize=10)
    
    # Direct Gene Layers (Sample)
    direct_positions = [(0.5, 6.5), (1.5, 6.5), (2.5, 6.5), (3.5, 6.5)]
    direct_labels = ['GO:0006805\ndirect_gene', 'GO:0008150\ndirect_gene', 'GO:0009987\ndirect_gene', '... 59 more']
    
    for i, (x, y) in enumerate(direct_positions[:4]):
        if i < 3:
            box = FancyBboxPatch((x, y), 1, 0.8, boxstyle="round,pad=0.05", 
                                facecolor=direct_color, edgecolor='black', linewidth=1)
            ax.add_patch(box)
            ax.text(x+0.5, y+0.4, direct_labels[i], ha='center', va='center', fontsize=8)
        else:
            ax.text(x+0.5, y+0.4, direct_labels[i], ha='center', va='center', fontsize=8, style='italic')
    
    # Term Processing Layers (Hierarchical)
    # Level 1 - Specific processes
    level1_pos = [(0.5, 5), (1.5, 5), (2.5, 5), (3.5, 5), (4.5, 5)]
    level1_labels = ['Metabolic\nProcesses', 'Transport\nProcesses', 'Regulation\nProcesses', 'Cellular\nProcesses', '...']
    
    for i, (x, y) in enumerate(level1_pos):
        if i < 4:
            # Main processing box
            box = FancyBboxPatch((x, y), 1, 1.2, boxstyle="round,pad=0.05", 
                                facecolor=term_color, edgecolor='black', linewidth=1)
            ax.add_patch(box)
            ax.text(x+0.5, y+0.9, level1_labels[i], ha='center', va='center', fontsize=8, fontweight='bold')
            ax.text(x+0.5, y+0.5, '50 hidden', ha='center', va='center', fontsize=7)
            ax.text(x+0.5, y+0.1, 'BatchNorm', ha='center', va='center', fontsize=7)
        else:
            ax.text(x+0.5, y+0.6, level1_labels[i], ha='center', va='center', fontsize=8, style='italic')
    
    # Level 2 - General processes (with dropout)
    level2_pos = [(1, 3), (2.5, 3), (4, 3)]
    level2_labels = ['General\nMetabolism', 'Transport &\nLocalization', 'Biological\nRegulation']
    
    for i, (x, y) in enumerate(level2_pos):
        # Dropout layer
        dropout_box = FancyBboxPatch((x, y+1.5), 1, 0.3, boxstyle="round,pad=0.02", 
                                    facecolor=dropout_color, edgecolor='black', linewidth=1)
        ax.add_patch(dropout_box)
        ax.text(x+0.5, y+1.65, 'Dropout 0.1', ha='center', va='center', fontsize=7)
        
        # Main processing box
        box = FancyBboxPatch((x, y), 1, 1.2, boxstyle="round,pad=0.05", 
                            facecolor=term_color, edgecolor='black', linewidth=1)
        ax.add_patch(box)
        ax.text(x+0.5, y+0.9, level2_labels[i], ha='center', va='center', fontsize=8, fontweight='bold')
        ax.text(x+0.5, y+0.5, '50 hidden', ha='center', va='center', fontsize=7)
        ax.text(x+0.5, y+0.1, 'BatchNorm', ha='center', va='center', fontsize=7)
    
    # ROOT Level
    root_box = FancyBboxPatch((2, 0.5), 1.5, 1.2, boxstyle="round,pad=0.1", 
                             facecolor=root_color, edgecolor='black', linewidth=2)
    ax.add_patch(root_box)
    ax.text(2.75, 1.3, 'ROOT', ha='center', va='center', fontsize=12, fontweight='bold')
    ax.text(2.75, 1.0, '1750 inputs', ha='center', va='center', fontsize=8)
    ax.text(2.75, 0.7, '50 hidden', ha='center', va='center', fontsize=8)
    
    # Final Output
    final_box = FancyBboxPatch((6, 0.5), 1.5, 1.2, boxstyle="round,pad=0.1", 
                              facecolor=final_color, edgecolor='black', linewidth=2)
    ax.add_patch(final_box)
    ax.text(6.75, 1.3, 'Final Output', ha='center', va='center', fontsize=12, fontweight='bold')
    ax.text(6.75, 1.0, '50 → 1 → 1', ha='center', va='center', fontsize=8)
    ax.text(6.75, 0.7, 'Prediction', ha='center', va='center', fontsize=8)
    
    # Arrows showing data flow
    # Gene input to direct layers
    for x, y in direct_positions[:3]:
        ax.annotate('', xy=(x+0.5, y+0.8), xytext=(2, 8),
                   arrowprops=dict(arrowstyle='->', lw=1, color='gray'))
    
    # Direct layers to level 1
    for i in range(3):
        ax.annotate('', xy=(level1_pos[i][0]+0.5, level1_pos[i][1]+1.2), 
                   xytext=(direct_positions[i][0]+0.5, direct_positions[i][1]),
                   arrowprops=dict(arrowstyle='->', lw=1, color='gray'))
    
    # Level 1 to level 2
    ax.annotate('', xy=(1.5, 4.8), xytext=(1, 5),
               arrowprops=dict(arrowstyle='->', lw=1.5, color='blue'))
    ax.annotate('', xy=(1.5, 4.8), xytext=(2, 5),
               arrowprops=dict(arrowstyle='->', lw=1.5, color='blue'))
    ax.annotate('', xy=(3, 4.8), xytext=(2.5, 5),
               arrowprops=dict(arrowstyle='->', lw=1.5, color='blue'))
    ax.annotate('', xy=(3, 4.8), xytext=(3.5, 5),
               arrowprops=dict(arrowstyle='->', lw=1.5, color='blue'))
    
    # Level 2 to ROOT
    for x, y in level2_pos:
        ax.annotate('', xy=(2.75, 1.7), xytext=(x+0.5, y),
                   arrowprops=dict(arrowstyle='->', lw=2, color='orange'))
    
    # ROOT to Final
    ax.annotate('', xy=(6, 1.1), xytext=(3.5, 1.1),
               arrowprops=dict(arrowstyle='->', lw=3, color='red'))
    
    # Auxiliary outputs (showing multi-level predictions)
    aux_y = [6, 4.5, 2.5]
    aux_labels = ['Term-level\nPredictions', 'Pathway-level\nPredictions', 'System-level\nPrediction']
    
    for i, (y, label) in enumerate(zip(aux_y, aux_labels)):
        aux_box = FancyBboxPatch((8, y), 1.5, 0.8, boxstyle="round,pad=0.05", 
                                facecolor='lightgray', edgecolor='black', linewidth=1, alpha=0.7)
        ax.add_patch(aux_box)
        ax.text(8.75, y+0.4, label, ha='center', va='center', fontsize=8)
        
        # Arrows from main pathway
        if i == 0:  # From level 1
            ax.annotate('', xy=(8, y+0.4), xytext=(5.5, 5.6),
                       arrowprops=dict(arrowstyle='->', lw=1, color='gray', linestyle='--'))
        elif i == 1:  # From level 2
            ax.annotate('', xy=(8, y+0.4), xytext=(5, 3.6),
                       arrowprops=dict(arrowstyle='->', lw=1, color='gray', linestyle='--'))
        else:  # From root
            ax.annotate('', xy=(8, y+0.4), xytext=(3.5, 1.1),
                       arrowprops=dict(arrowstyle='->', lw=1, color='gray', linestyle='--'))
    
    # Add title and labels
    ax.text(5, 9.5, 'DCellNN Architecture Overview', ha='center', va='center', 
            fontsize=16, fontweight='bold')
    
    # Legend
    legend_elements = [
        patches.Patch(color=gene_color, label='Gene Input (2 genes)'),
        patches.Patch(color=direct_color, label='Direct Gene Layers (62 terms)'),
        patches.Patch(color=term_color, label='GO Term Processing'),
        patches.Patch(color=dropout_color, label='Dropout Regularization'),
        patches.Patch(color=root_color, label='ROOT Integration'),
        patches.Patch(color=final_color, label='Final Prediction')
    ]
    
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(0.02, 0.98))
    
    # Add annotations
    ax.text(5, -0.5, 'Key Features:\n• Hierarchical biological pathway modeling\n• Multi-level predictions (auxiliary outputs)\n• Sparse connectivity based on Gene Ontology\n• 50 hidden units per term, ~3,100 core parameters', 
            ha='center', va='top', fontsize=10, 
            bbox=dict(boxstyle="round,pad=0.5", facecolor='lightyellow', alpha=0.8))
    
    # Set axis properties
    ax.set_xlim(-0.5, 10)
    ax.set_ylim(-1.5, 10)
    ax.set_aspect('equal')
    ax.axis('off')
    
    plt.tight_layout()
    return fig

if __name__ == "__main__":
    fig = create_network_diagram()
    
    # Save the diagram
    plt.savefig('/work/users/m/j/mjn15/bct_variants/dcell_architecture_diagram.png', 
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('/work/users/m/j/mjn15/bct_variants/dcell_architecture_diagram.pdf', 
                bbox_inches='tight', facecolor='white')
    
    print("Network diagram saved as:")
    print("- dcell_architecture_diagram.png")
    print("- dcell_architecture_diagram.pdf")
    
    plt.show()
