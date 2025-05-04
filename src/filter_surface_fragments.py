#!/usr/bin/env python3
"""
Filter surface fragment analysis by a cutoff distance.
Reads a combined surface analysis CSV produced by fragment_surface_analysis.py,
filters per-fragment-per-frame pairs where the minimum atom distance to Pt surface ≤ cutoff,
and from those, selects per-atom entries within cutoff.

Produces combined outputs in the output directory:
    - surface_filtered.csv: rows for atoms within cutoff in qualifying fragment-frame pairs
    - surface_filtered_summary.json: list of summary dicts for each qualifying fragment-frame
    - surface_filtered_summary.md: Markdown report of filtered fragment-frame results

Usage:
    python filter_surface_fragments.py \
        --input-csv /path/to/surface_analysis.csv \
        --output-dir /path/to/output_dir \
        [--cutoff 3.0]
"""
import os
import argparse
import json
import pandas as pd
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

def main():
    parser = argparse.ArgumentParser(
        description="Filter surface fragment analysis by cutoff distance per frame."
    )
    parser.add_argument(
        "--input-csv", required=True,
        help="Combined CSV from fragment_surface_analysis.py"
    )
    parser.add_argument(
        "--output-dir", default=".",
        help="Directory to save filtered outputs."
    )
    parser.add_argument(
        "--cutoff", type=float, default=3.0,
        help="Distance cutoff (Å) to select fragment-frame pairs and atoms."
    )
    args = parser.parse_args()

    # Load data
    df = pd.read_csv(args.input_csv)
    required = {"fragment_id", "frame_index", "pt_min_distance", "pt_classification"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {missing}")

    # Identify fragment-frame minima
    minima = (
        df.groupby(["fragment_id", "frame_index"])["pt_min_distance"]
          .min()
          .reset_index(name="min_distance")
    )
    # Keep only pairs under cutoff
    keeps = minima[minima["min_distance"] <= args.cutoff].assign(keep=True)
    if keeps.empty:
        print(f"No fragment/frame pairs with min_distance ≤ {args.cutoff} Å; no output generated.")
        return

    # Merge and filter
    df2 = df.merge(
        keeps[['fragment_id', 'frame_index', 'keep']],
        on=["fragment_id", "frame_index"], how="left"
    )
    df2['keep'] = df2['keep'].fillna(False)
    df_filtered = (
        df2.loc[df2['keep'] & (df2['pt_min_distance'] <= args.cutoff)]
        .drop(columns=["keep"])
        .copy()
    )

    # Summaries per fragment-frame
    summaries = []
    for _, row in keeps.iterrows():
        frag = int(row['fragment_id'])
        frame = int(row['frame_index'])
        min_dist = float(row['min_distance'])
        sub_all = df[(df['fragment_id']==frag) & (df['frame_index']==frame)]
        sub_within = df_filtered[(df_filtered['fragment_id']==frag) & (df_filtered['frame_index']==frame)]

        total_atoms = int(len(sub_all))
        n_within = int(len(sub_within))
        class_counts = (
            sub_within['pt_classification']
            .value_counts()
            .to_dict()
        )
        class_counts = {str(k): int(v) for k, v in class_counts.items()}

        summaries.append({
            'fragment_id': frag,
            'frame_index': frame,
            'cutoff': args.cutoff,
            'min_distance': min_dist,
            'total_atoms': total_atoms,
            'atoms_within_cutoff': n_within,
            'pt_classification_within_counts': class_counts
        })

    # Ensure output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Write outputs
    csv_out = os.path.join(args.output_dir, "surface_filtered.csv")
    df_filtered.to_csv(csv_out, index=False)
    print(f"Saved filtered CSV: {csv_out}")

    json_out = os.path.join(args.output_dir, "surface_filtered_summary.json")
    with open(json_out, 'w') as fj:
        json.dump(summaries, fj, indent=2)
    print(f"Saved JSON summary: {json_out}")

    md_out = os.path.join(args.output_dir, "surface_filtered_summary.md")
    with open(md_out, 'w') as fm:
        fm.write(f"# Filtered Surface Analysis (cutoff = {args.cutoff} Å)\n\n")
        for s in summaries:
            fm.write(f"## Fragment {s['fragment_id']} (Frame {s['frame_index']})\n")
            fm.write(f"- Min distance: **{s['min_distance']:.2f} Å**\n")
            fm.write(f"- Total NEC atoms: **{s['total_atoms']}**\n")
            fm.write(f"- Atoms within cutoff: **{s['atoms_within_cutoff']}**\n")
            fm.write(f"- pt_classification counts within cutoff:\n")
            for cls, cnt in s['pt_classification_within_counts'].items():
                fm.write(f"  - {cls}: {cnt}\n")
            fm.write("\n")
    print(f"Saved Markdown report: {md_out}")

if __name__ == '__main__':
    main()
