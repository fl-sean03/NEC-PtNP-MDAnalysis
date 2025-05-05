#!/usr/bin/env python3
"""
Compute binding metrics for NEC fragments on the Pt nanoparticle surface.

This script calculates per-molecule and per-facet binding statistics,
including mean residence time, dissociation constant (K_D), and free energy (ΔG).

It can analyze the full processed range or a specific timeframe within it.

Requires:
  - CSV of residence events (from residence_event_analysis.py)
  - CSV of filtered surface analysis data (from filter_surface_fragments.py)
  - JSON of simulation properties (from simulation_info.py)

Outputs:
  - CSV with per-molecule binding metrics
  - CSV with per-facet binding metrics

Usage (Full Range):
  python binding_metrics.py \\
      --residence_events_csv /path/to/residence_events.csv \\
      --filtered_surface_csv /path/to/surface_filtered.csv \\
      --sim-info-json /path/to/simulation_info.json \\
      [--output_molecule_csv molecule_results.csv] \\
      [--output_facet_csv facet_results.csv]

Usage (Specific Timeframe):
  python binding_metrics.py \\
      --residence_events_csv /path/to/residence_events.csv \\
      --filtered_surface_csv /path/to/surface_filtered.csv \\
      --sim-info-json /path/to/simulation_info.json \\
      --analysis-start-frame <start_frame> \\
      --analysis-end-frame <end_frame> \\
      --min-res-time <min_time_ns> \\
      [--output_molecule_csv molecule_results.csv] \\
      [--output_facet_csv facet_results.csv]
"""
import argparse
import pandas as pd
import numpy as np
import os
import json

def main():
    parser = argparse.ArgumentParser(
        description="Compute binding metrics from event and filtered surface analysis data."
    )
    # Input Files
    parser.add_argument(
        "--residence_events_csv", required=True, help="Path to residence events CSV (from residence_event_analysis.py)"
    )
    parser.add_argument(
        "--filtered_surface_csv", required=True, help="Path to filtered surface analysis CSV (from filter_surface_fragments.py)"
    )
    parser.add_argument(
        "--sim-info-json", required=True, help="Path to simulation_info.json with time info (ns)."
    )
    # Output Files
    parser.add_argument(
        "--output_molecule_csv", default="molecule_results.csv",
        help="Output path for per-molecule CSV"
    )
    parser.add_argument(
        "--output_facet_csv", default="facet_results.csv",
        help="Output path for per-facet CSV"
    )
    # Timeframe Analysis Arguments
    parser.add_argument(
        "--analysis-start-frame", type=int, default=None,
        help="Optional start frame index (inclusive) for the analysis window."
    )
    parser.add_argument(
        "--analysis-end-frame", type=int, default=None,
        help="Optional end frame index (exclusive) for the analysis window."
    )
    parser.add_argument(
        "--min-res-time", type=float, default=None,
        help="Minimum residence time (in ns) threshold used in residence_event_analysis. Required if timeframe filtering is active."
    )
    # Legacy argument (not used but kept for potential compatibility if needed)
    parser.add_argument(
        "--simulation_properties_csv", help="Legacy path to simulation properties CSV (use --sim-info-json instead)."
    )

    args = parser.parse_args()

    print("--- Running Binding Metrics Calculation ---")
    print(f"Loading data from: {args.residence_events_csv}, {args.filtered_surface_csv}, {args.sim_info_json}")

    # --- 1. Load Data ---
    events = pd.read_csv(args.residence_events_csv)
    filtered_surface_data = pd.read_csv(args.filtered_surface_csv)
    with open(args.sim_info_json, 'r') as f:
        sim_info_data = json.load(f)

    print("Data loaded. Performing calculations...")

    # --- 2. Determine Timeframe and T_sim ---
    timeframe_analysis_active = args.analysis_start_frame is not None and args.analysis_end_frame is not None

    if timeframe_analysis_active:
        # --- 2a. Timeframe Analysis ---
        start_frame = args.analysis_start_frame
        end_frame = args.analysis_end_frame # Exclusive

        if args.min_res_time is None:
            raise ValueError("--min-res-time is required when --analysis-start-frame and --analysis-end-frame are provided.")
        if start_frame >= end_frame:
             raise ValueError("Analysis end frame must be greater than start frame.")

        print(f"Analyzing timeframe: frames [{start_frame}, {end_frame})")

        # Get time step from sim_info_json
        if "full_simulation" not in sim_info_data or "time_between_frames_ns" not in sim_info_data["full_simulation"]:
            raise KeyError("simulation_info.json missing 'full_simulation' or 'time_between_frames_ns'.")
        dt_ns = float(sim_info_data["full_simulation"]["time_between_frames_ns"])

        # Calculate effective simulation time (T_sim) for the analysis window
        analysis_window_frames = end_frame - start_frame
        T_sim = analysis_window_frames * dt_ns
        print(f"Effective simulation time for analysis window: {T_sim:.3f} ns")

        # Filter events to include only those overlapping the analysis timeframe
        print(f"Events before timeframe filtering: {len(events)}")
        events_in_range = events[
            (events['start_frame'] < end_frame) &
            (events['end_frame'] >= start_frame) # Use >= for end_frame as it's inclusive in events data
        ].copy()
        print(f"Events after timeframe filtering: {len(events_in_range)}")

        # Adjust event durations to be strictly within the analysis timeframe
        events_in_range['clipped_start_frame'] = events_in_range['start_frame'].clip(lower=start_frame)
        # Clip end_frame to be one less than the exclusive analysis end frame
        events_in_range['clipped_end_frame'] = events_in_range['end_frame'].clip(upper=end_frame - 1)

        # Calculate clipped duration in frames (inclusive start, inclusive end) and convert to ns
        # Add 1 because both clipped_start_frame and clipped_end_frame are inclusive indices
        events_in_range['clipped_duration_frames'] = (events_in_range['clipped_end_frame'] - events_in_range['clipped_start_frame']) + 1
        events_in_range['clipped_duration_ns'] = events_in_range['clipped_duration_frames'] * dt_ns

        # Filter out events with zero or negative duration after clipping
        events_in_range = events_in_range[events_in_range['clipped_duration_ns'] > 1e-9].copy() # Use small tolerance
        print(f"Events after clipping and removing zero duration: {len(events_in_range)}")

        # Filter out events whose clipped duration falls below the minimum residence time
        events_in_range = events_in_range[events_in_range['clipped_duration_ns'] >= args.min_res_time].copy()
        print(f"Events after filtering by min_res_time ({args.min_res_time} ns): {len(events_in_range)}")

        # Use clipped duration for subsequent calculations
        events_in_range['duration_ns'] = events_in_range['clipped_duration_ns']
        events = events_in_range[['fragment_id', 'event_id', 'start_frame', 'end_frame', 'start_time_ns', 'end_time_ns', 'duration_ns']].copy() # Select relevant columns

    else:
        # --- 2b. Full Range Analysis ---
        if "processed_range" not in sim_info_data or "total_time_ns" not in sim_info_data["processed_range"]:
             raise KeyError("simulation_info.json missing 'processed_range' or 'total_time_ns'.")
        sim_range = sim_info_data["processed_range"]
        T_sim = float(sim_range.get("total_time_ns"))
        print(f"Using total simulation time from processed range: {T_sim:.3f} ns")

    # --- 3. Facet Mapping ---
    # Ensure fragment_id is treated as the same type for merging
    events['fragment_id'] = events['fragment_id'].astype(str)
    filtered_surface_data['fragment_id'] = filtered_surface_data['fragment_id'].astype(str)

    # Determine the most frequent pt_classification for each fragment_id
    # This mapping uses the full filtered_surface_data regardless of the analysis timeframe
    if not filtered_surface_data.empty:
        fragment_facet_mapping = filtered_surface_data.groupby('fragment_id')['pt_classification'].agg(
            lambda x: x.mode()[0] if not x.mode().empty else None
        ).reset_index(name='dominant_facet')
    else:
        fragment_facet_mapping = pd.DataFrame(columns=['fragment_id', 'dominant_facet'])
        print("Warning: Filtered surface data is empty. Cannot determine fragment-facet mapping.")

    # --- 4. Merge Events with Facets ---
    if events.empty:
         print("No residence events found within the specified timeframe or overall. Skipping calculations.")
         df = pd.DataFrame(columns=['fragment_id', 'event_id', 'duration_ns', 'facet']) # Empty df with needed columns
    else:
        df = events.merge(
            fragment_facet_mapping,
            on='fragment_id',
            how='left'
        )
        df.rename(columns={'dominant_facet': 'facet'}, inplace=True)

    # --- 5. Per-Molecule Statistics ---
    if df.empty or df['facet'].isnull().all():
         print("No fragments with facet classification found in the analyzed events. Skipping per-molecule and per-facet stats.")
         mol = pd.DataFrame(columns=['fragment_id', 'T_on_ns', 'N_events', 'facet', 'tau_mean_ns', 'T_off_ns', 'K_D'])
    else:
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
            print("Aggregation resulted in empty per-molecule data. Skipping stats.")
        else:
            mol['tau_mean_ns'] = mol['T_on_ns'] / mol['N_events']
            mol['T_off_ns'] = T_sim - mol['T_on_ns']
            # Handle cases where T_on_ns is zero or T_sim is less than T_on_ns (due to potential float issues)
            mol['K_D'] = mol.apply(
                lambda row: row['T_off_ns'] / row['T_on_ns'] if row['T_on_ns'] > 1e-9 and row['T_off_ns'] >= 0 else np.inf,
                axis=1
            )
            # Set K_D to 0 if T_off is effectively zero (fully bound within the window)
            mol.loc[mol['T_off_ns'] < 1e-9, 'K_D'] = 0.0

            print(f"Computed per-molecule stats for {len(mol)} fragments.")

    # --- 6. Per-Facet Statistics ---
    if mol.empty:
        print("Skipping per-facet stats as per-molecule data is empty.")
        facet = pd.DataFrame(columns=['facet', 'tau_facet_ns', 'K_D_facet', 'DeltaG_kJ_per_mol'])
    else:
        facet = (
            mol.groupby('facet')
               .agg(
                   tau_facet_ns=('tau_mean_ns', 'mean'),
                   # Calculate mean K_D excluding infinite values
                   K_D_facet=('K_D', lambda x: np.mean(x[np.isfinite(x)]))
               )
               .reset_index()
        )

        # Compute free energy ΔG = R T ln(K_D)
        R = 8.314    # J/(mol·K)
        T = 453.0    # K (Hardcoded temperature - consider making this an input)
        # Compute DeltaG only for facets with finite and positive K_D_facet
        facet['DeltaG_kJ_per_mol'] = facet.apply(
            lambda row: (R * T * np.log(row['K_D_facet'])) / 1000 if pd.notna(row['K_D_facet']) and np.isfinite(row['K_D_facet']) and row['K_D_facet'] > 0 else np.nan,
            axis=1
        )
        print(f"Computed per-facet stats for {len(facet)} facets.")


    # --- 7. Export Results ---
    # Ensure output directories exist
    os.makedirs(os.path.dirname(args.output_molecule_csv) or '.', exist_ok=True)
    os.makedirs(os.path.dirname(args.output_facet_csv) or '.', exist_ok=True)

    # Save molecule results
    mol_cols = ['fragment_id', 'T_on_ns', 'N_events', 'facet', 'tau_mean_ns', 'T_off_ns', 'K_D']
    mol_final = pd.DataFrame(mol, columns=mol_cols) # Ensure columns exist even if empty
    mol_final.to_csv(args.output_molecule_csv, index=False, float_format='%.6g')
    print(f"Per-molecule results saved to {args.output_molecule_csv}")

    # Save facet results
    facet_cols = ['facet', 'tau_facet_ns', 'K_D_facet', 'DeltaG_kJ_per_mol']
    facet_final = pd.DataFrame(facet, columns=facet_cols) # Ensure columns exist even if empty
    facet_final.to_csv(args.output_facet_csv, index=False, float_format='%.6g')
    print(f"Per-facet results saved to {args.output_facet_csv}")

    print("--- Binding Metrics Calculation Finished ---")


if __name__ == "__main__":
    main()