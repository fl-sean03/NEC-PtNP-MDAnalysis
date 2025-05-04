#!/usr/bin/env python3
"""
Compute binding metrics for NEC fragments on the Pt nanoparticle surface.

This script calculates per-molecule and per-facet binding statistics,
including mean residence time, dissociation constant (K_D), and free energy (ΔG).

It requires the following input files:
  - CSV of residence events (from residence_event_analysis.py)
  - CSV of filtered surface analysis data (from filter_surface_fragments.py)
  - CSV of simulation properties (from simulation_info.py)

Outputs:
  - CSV with per-molecule binding metrics
  - CSV with per-facet binding metrics

Usage:
  python binding_metrics.py \\
      --residence_events_csv /path/to/residence_events.csv \\
      --filtered_surface_csv /path/to/surface_filtered.csv \\
      --simulation_properties_csv /path/to/simulation_info.csv \\
      [--output_molecule_csv molecule_results.csv] \\
      [--output_facet_csv facet_results.csv]
"""
import argparse
import pandas as pd
import numpy as np
import os

def main():
    parser = argparse.ArgumentParser(
        description="Compute binding metrics from event and filtered surface analysis CSVs"
    )
    parser.add_argument(
        "--residence_events_csv", required=True, help="Path to residence events CSV (from residence_event_analysis.py)"
    )
    parser.add_argument(
        "--filtered_surface_csv", required=True, help="Path to filtered surface analysis CSV (from filter_surface_fragments.py)"
    )
    parser.add_argument(
        "--simulation_properties_csv", required=True, help="Path to simulation properties CSV (from simulation_info.py)"
    )
    parser.add_argument(
        "--output_molecule_csv", default="molecule_results.csv",
        help="Output path for per-molecule CSV"
    )
    parser.add_argument(
        "--output_facet_csv", default="facet_results.csv",
        help="Output path for per-facet CSV"
    )
    args = parser.parse_args()

    print("--- Running Binding Metrics Calculation ---")
    print(f"Loading data from: {args.residence_events_csv}, {args.filtered_surface_csv}, {args.simulation_properties_csv}")

    # Read CSVs
    events = pd.read_csv(args.residence_events_csv)
    filtered_surface_data = pd.read_csv(args.filtered_surface_csv)
    props = pd.read_csv(args.simulation_properties_csv)

    print("Data loaded. Performing calculations...")

    # Ensure fragment_id is treated as the same type for merging
    events['fragment_id'] = events['fragment_id'].astype(str)
    filtered_surface_data['fragment_id'] = filtered_surface_data['fragment_id'].astype(str)

    # Determine the most frequent pt_classification for each fragment_id from the filtered surface data
    if not filtered_surface_data.empty:
        fragment_facet_mapping = filtered_surface_data.groupby('fragment_id')['pt_classification'].agg(lambda x: x.mode()[0] if not x.mode().empty else None).reset_index(name='dominant_facet')
    else:
        fragment_facet_mapping = pd.DataFrame(columns=['fragment_id', 'dominant_facet'])
        print("Filtered surface data is empty. Cannot determine fragment-facet mapping.")


    # Merge events with the dominant facet classification
    # We need to merge events with the dominant facet for each fragment.
    # The original script seemed to imply a single facet per fragment for the 'facet' column.
    # Let's merge events with the dominant facet mapping.
    df = events.merge(
        fragment_facet_mapping,
        on='fragment_id',
        how='left'
    )

    # Rename dominant_facet to facet for consistency with original script's output structure
    df.rename(columns={'dominant_facet': 'facet'}, inplace=True)


    # Simulation time (ns)
    T_sim_row = props.loc[props['property'] == 'total_time_ns', 'value']
    if T_sim_row.empty:
        raise ValueError("Could not find 'total_time_ns' in properties CSV.")
    T_sim = float(T_sim_row.iloc[0])
    print(f"Total simulation time: {T_sim:.3f} ns")

    # Per-molecule stats
    # Filter out fragments that didn't merge with classification data if necessary
    mol = (
        df.dropna(subset=['facet']) # Ensure we only process fragments with a determined facet
          .groupby('fragment_id')
          .agg(
              T_on_ns=('duration_ns', 'sum'),
              N_events=('duration_ns', 'count'),
              facet=('facet', 'first') # 'first' is safe after grouping by fragment_id and dropping NaNs
          )
          .reset_index()
    )

    if mol.empty:
        print("No fragments with facet classification found in events data. Skipping per-molecule and per-facet stats.")
        # Still create empty output files to avoid pipeline errors
        os.makedirs(os.path.dirname(args.output_molecule_csv) or '.', exist_ok=True)
        os.makedirs(os.path.dirname(args.output_facet_csv) or '.', exist_ok=True)
        pd.DataFrame(columns=['fragment_id', 'T_on_ns', 'N_events', 'facet', 'tau_mean_ns', 'T_off_ns', 'K_D']).to_csv(args.output_molecule_csv, index=False)
        pd.DataFrame(columns=['facet', 'tau_facet_ns', 'K_D_facet', 'DeltaG_kJ_per_mol']).to_csv(args.output_facet_csv, index=False)
        print(f"Created empty output files: {args.output_molecule_csv}, {args.output_facet_csv}")
        return


    mol['tau_mean_ns'] = mol['T_on_ns'] / mol['N_events']
    mol['T_off_ns']   = T_sim - mol['T_on_ns']
    # Handle cases where T_on_ns is zero to avoid division by zero for K_D
    mol['K_D'] = mol.apply(lambda row: row['T_off_ns'] / row['T_on_ns'] if row['T_on_ns'] > 0 else np.inf, axis=1)

    print(f"Computed per-molecule stats for {len(mol)} fragments.")

    # Per-facet stats
    facet = (
        mol.groupby('facet')
           .agg(
               tau_facet_ns=('tau_mean_ns', 'mean'),
               K_D_facet=('K_D', lambda x: np.mean(x[np.isfinite(x)])) # Exclude inf K_D for mean
           )
           .reset_index()
    )

    # Compute free energy ΔG = R T ln(K_D)
    R = 8.314    # J/(mol·K)
    T = 453.0    # K (Hardcoded temperature from original script)
    # Compute DeltaG only for facets with finite K_D_facet
    facet['DeltaG_kJ_per_mol'] = facet.apply(
        lambda row: (R * T * np.log(row['K_D_facet'])) / 1000 if np.isfinite(row['K_D_facet']) and row['K_D_facet'] > 0 else np.nan,
        axis=1
    )

    print(f"Computed per-facet stats for {len(facet)} facets.")

    # Export results
    os.makedirs(os.path.dirname(args.output_molecule_csv) or '.', exist_ok=True)
    os.makedirs(os.path.dirname(args.output_facet_csv) or '.', exist_ok=True)

    mol.to_csv(args.output_molecule_csv, index=False)
    facet.to_csv(args.output_facet_csv, index=False)

    print(f"Per-molecule results saved to {args.output_molecule_csv}")
    print(f"Per-facet results saved to {args.output_facet_csv}")
    print("--- Binding Metrics Calculation Finished ---")


if __name__ == "__main__":
    main()
