#!/usr/bin/env python3
"""
Script to analyze whether NEC fragments lie on the Pt nanoparticle surface,
with frame-level checkpointing to resume from the last processed frame.
Produces:
    - surface_analysis.csv
    - surface_summary.json
    - surface_summary.md

Usage:
    python fragment_surface_analysis.py \
        --psf /path/to/system.psf \
        --dcd /path/to/trajectory.dcd \
        --pt-csv /path/to/coordination_numbers.csv \
        --output-dir /path/to/output_dir \
        [--fragment-ids 0 2 5] \
        [--frame-index 0] \
        [--resume] \
        [-v|--verbose] \
        [--all-atoms]
"""
import os
import argparse
import json
import numpy as np
import pandas as pd
import MDAnalysis as mda
from scipy.spatial import cKDTree
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

def main():
    parser = argparse.ArgumentParser(
        description="Surface analysis for NEC fragments with checkpointing."
    )
    parser.add_argument("--psf", required=True, help="Path to the PSF file.")
    parser.add_argument("--dcd", required=True, help="Path to the DCD trajectory file.")
    parser.add_argument(
        "--pt-csv", required=True,
        help="CSV of Pt atom classifications (coordination numbers + is_surface)."
    )
    parser.add_argument(
        "--output-dir", default=".",
        help="Directory to save outputs (surface_analysis.csv, .json, .md)."
    )
    parser.add_argument(
        "--fragment-ids", type=int, nargs="*", default=None,
        help="IDs of fragments to analyze; if omitted, analyze all."
    )
    parser.add_argument(
        "--frame-index", type=int,
        help="Frame index to analyze; if set, only that frame will run (useful for testing)."
    )
    parser.add_argument(
        "--resume", action="store_true",
        help="Resume from the last processed frame in surface_analysis.csv."
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable verbose output to terminal."
    )
    parser.add_argument(
        "--all-atoms", action="store_true",
        help="Save data for all NEC atoms; default is only the closest per fragment."
    )
    args = parser.parse_args()

    def vprint(msg):
        if args.verbose:
            print(msg)

    # Prepare output paths
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    analysis_csv = os.path.join(output_dir, "surface_analysis.csv")
    summary_json = os.path.join(output_dir, "surface_summary.json")
    summary_md   = os.path.join(output_dir, "surface_summary.md")

    # Debug prints for resume logic
    print(f"Checking for resume file: {analysis_csv}")
    print(f"File exists: {os.path.exists(analysis_csv)}")
    if os.path.exists(analysis_csv):
        print(f"File size: {os.path.getsize(analysis_csv)}")

    # Determine starting frame based on checkpoint in analysis_csv
    if os.path.exists(analysis_csv) and os.path.getsize(analysis_csv) > 0:
        with open(analysis_csv, 'r') as f:
            f.readline()  # skip header
            last_line = None
            for line in f:
                last_line = line
        if last_line:
            last_frame = int(last_line.split(',')[1])
            start_frame = last_frame + 1
            write_header = False
            mode = 'a'
            vprint(f"Resuming from frame {start_frame} (last recorded: {last_frame})")
        else:
            start_frame = 0
            write_header = True
            mode = 'w'
            vprint("Checkpoint file exists but empty, starting at frame 0")
    else:
        start_frame = 0
        write_header = True
        mode = 'w'
        vprint("Starting fresh run from frame 0")

    # Load Pt classification
    df_pt = pd.read_csv(args.pt_csv)
    pt_surface = df_pt[df_pt["is_surface"]]
    vprint(f"Loaded Pt classification: {len(pt_surface)} surface atoms")

    # Load Universe
    u = mda.Universe(args.psf, args.dcd)
    n_frames = len(u.trajectory)
    vprint(f"Loaded trajectory: {n_frames} frames total")

    # Determine fragments
    nec_atoms = u.select_atoms("not name PT*")
    fragments = nec_atoms.fragments
    total_frags = len(fragments)
    fragment_ids = list(range(total_frags)) if args.fragment_ids is None else [f for f in args.fragment_ids if 0 <= f < total_frags]
    vprint(f"Fragments to analyze: {fragment_ids}")

    # Static Pt selections
    all_pt = u.select_atoms("name PT*")
    surf_idx = pt_surface["atom_index"].values
    pt_classes = pt_surface["classification"].values

    # Prepare summary list
    summaries = []

    # Open analysis CSV for writing/appending
    f_csv = open(analysis_csv, mode)
    if write_header:
        f_csv.write("fragment_id,frame_index,atom_index,pt_min_distance,pt_classification\n")
        f_csv.flush()

    # Frame iteration
    if args.frame_index is not None:
        if not (0 <= args.frame_index < n_frames):
            raise ValueError(f"frame-index {args.frame_index} out of range (0–{n_frames-1})")
        frames = [u.trajectory[args.frame_index]]
        vprint(f"Isolated single frame {args.frame-index}")
    else:
        frames = u.trajectory[start_frame:]

    for ts in frames:
        frame = ts.frame
        vprint(f"Processing frame {frame}")

        tree = cKDTree(all_pt.positions[surf_idx])
        frame_rows = []

        for frag_id in fragment_ids:
            frag = fragments[frag_id]
            dist_arr, nn_idx = tree.query(frag.positions, k=1)
            nearest_cls = pt_classes[nn_idx]

            if args.all_atoms:
                for idx, d, cls in zip(frag.indices, dist_arr, nearest_cls):
                    frame_rows.append(f"{frag_id},{frame},{idx},{d},{cls}\n")
            else:
                i0 = np.argmin(dist_arr)
                frame_rows.append(
                    f"{frag_id},{frame},{frag.indices[i0]},{dist_arr[i0]},{nearest_cls[i0]}\n"
                )

            class_counts = pd.Series(nearest_cls).value_counts().to_dict()
            summaries.append({
                "fragment_id": frag_id,
                "frame_index": frame,
                "n_atoms": int(len(dist_arr)),
                "min_distance": float(dist_arr.min()),
                "max_distance": float(dist_arr.max()),
                "mean_distance": float(dist_arr.mean()),
                "median_distance": float(np.median(dist_arr)),
                "std_distance": float(dist_arr.std()),
                "pt_classification_counts": {str(k): int(v) for k, v in class_counts.items()}
            })

        for row in frame_rows:
            f_csv.write(row)
        f_csv.flush()

    f_csv.close()

    # Write summary JSON and Markdown
    with open(summary_json, "w") as fj:
        json.dump(summaries, fj, indent=2)

    with open(summary_md, "w") as fm:
        fm.write("# Surface Analysis Report\n\n")
        for s in summaries:
            fm.write(f"## Fragment {s['fragment_id']} (Frame {s['frame_index']})\n")
            fm.write(f"- NEC atoms: **{s['n_atoms']}**\n")
            fm.write(f"- Min distance: **{s['min_distance']:.2f} Å**\n")
            fm.write(f"- Max distance: **{s['max_distance']:.2f} Å**\n")
            fm.write(f"- Mean distance: **{s['mean_distance']:.2f} Å**\n")
            fm.write(f"- Median distance: **{s['median_distance']:.2f} Å**\n")
            fm.write(f"- Std deviation: **{s['std_distance']:.2f} Å**\n")
            fm.write("- Pt classification counts:\n")
            for cls, cnt in s['pt_classification_counts'].items():
                fm.write(f"  - {cls}: {cnt}\n")
            fm.write("\n")

if __name__ == "__main__":
    main()