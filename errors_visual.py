import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def visualize_errors(file_path='demux_results/benchmark_trimmed_sequences.xlsx', output_dir='demux_results'):
    """
    Reads benchmark Excel file and visualizes mean Total_Penalty, Mismatches_Perc, and Gaps_Perc per SampleID.
    """
    if not os.path.exists(file_path):
        print(f"Error: File {file_path} not found.")
        return

    # Load data
    try:
        df = pd.read_excel(file_path)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return

    # Columns to analyze
    metrics = ['Total_Penalty', 'Mismatches_Perc', 'Gaps_Perc']
    
    # Check if columns exist
    for metric in metrics:
        if metric not in df.columns:
            print(f"Error: Column '{metric}' not found in {file_path}")
            return
    if 'SampleID' not in df.columns:
        print(f"Error: Column 'SampleID' not found in {file_path}")
        return

    # Set visualization style
    sns.set_theme(style="whitegrid")
    
    # Create a figure with 3 subplots
    fig, axes = plt.subplots(3, 1, figsize=(12, 18))
    
    titles = [
        'Средний штраф (Total Penalty) по группам',
        'Средний % несовпадений (Mismatches Perc) по группам',
        'Средний % гэпов (Gaps Perc) по группам'
    ]

    for i, metric in enumerate(metrics):
        # Calculate mean per SampleID and sort
        group_means = df.groupby('SampleID')[metric].mean().sort_values(ascending=False).reset_index()
        
        sns.barplot(
            data=group_means, 
            x='SampleID', 
            y=metric, 
            ax=axes[i], 
            palette='magma',
            hue='SampleID',
            legend=False
        )
        
        axes[i].set_title(titles[i], fontsize=14)
        # Use set_xticks to avoid warnings in newer matplotlib/seaborn
        axes[i].set_xticks(range(len(group_means['SampleID'])))
        axes[i].set_xticklabels(group_means['SampleID'], rotation=45, ha='right')
        axes[i].set_xlabel('Sample ID')
        axes[i].set_ylabel(f'Mean {metric}')

    plt.tight_layout()
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    output_path = os.path.join(output_dir, 'errors_visualization.png')
    plt.savefig(output_path, dpi=300)
    print(f"Visualization saved to {output_path}")

if __name__ == "__main__":
    visualize_errors()
