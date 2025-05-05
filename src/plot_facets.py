#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(
        description="Plot Pt facet metrics from CSV"
    )
    parser.add_argument(
        "csv_file",
        help="Path to CSV file containing facet metrics"
    )
    parser.add_argument(
        "-o", "--outdir",
        default=".",
        help="Directory to save plots (will be created if it doesn't exist)"
    )
    args = parser.parse_args()

    # Create output directory if needed
    os.makedirs(args.outdir, exist_ok=True)

    # Increase global font size for readability
    plt.rcParams.update({'font.size': 14})

    # Read the data
    df = pd.read_csv(args.csv_file)

    # Define metrics: column name → (axis label, plot file name)
    metrics = {
        "tau_facet_ns": ("Mean Residence Time τ (ns)", "tau_vs_frame.png"),
        "K_D_facet": ("Dissociation Constant K_D", "kd_vs_frame.png"),
        "DeltaG_kJ_per_mol": ("Free Energy ΔG (kJ/mol)", "dg_vs_frame.png"),
    }

    # Generate individual plots
    for col, (ylabel, filename) in metrics.items():
        fig, ax = plt.subplots(figsize=(8, 6))
        for facet in df["facet"].unique():
            subset = df[df["facet"] == facet]
            ax.plot(
                subset["end_frame"],
                subset[col],
                marker="o",
                label=facet
            )
        ax.set_xlabel("End Frame")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{ylabel} vs End Frame")
        ax.legend(fontsize=12)
        fig.tight_layout()
        fig.savefig(os.path.join(args.outdir, filename))
        plt.close(fig)

    # Generate combined plot
    fig, axes = plt.subplots(3, 1, figsize=(8, 18), sharex=True)
    for ax, (col, (ylabel, _)) in zip(axes, metrics.items()):
        for facet in df["facet"].unique():
            subset = df[df["facet"] == facet]
            ax.plot(
                subset["end_frame"],
                subset[col],
                marker="o",
                label=facet
            )
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel)
        ax.legend(fontsize=12)
    axes[-1].set_xlabel("End Frame")
    fig.tight_layout()
    combined_path = os.path.join(args.outdir, "all_metrics.png")
    fig.savefig(combined_path)
    plt.close(fig)

if __name__ == "__main__":
    main()

