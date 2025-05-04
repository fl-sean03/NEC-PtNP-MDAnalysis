#!/usr/bin/env python3
"""
Script to analyze whether NEC fragments lie on the Pt nanoparticle surface.
Computes, for each NEC fragment in the specified frames, either the single closest NEC atom or all NEC atoms,
and records the minimum distance to any Pt surface atom and the classification of that nearest Pt atom.

Can process a single frame (with --frame-index) or all frames (leave --frame-index unset).
By default, only the closest NEC atom per fragment per frame is written to CSV for reduced output;
use --all-atoms to save all per-atom distances.

Produces combined outputs (in the output directory):
    - surface_analysis.csv: CSV with columns fragment_id, frame_index, atom_index, pt_min_distance, pt_classification
    - surface_summary.json: JSON list of per-fragment-per-frame summary dicts
    - surface_summary.md: Markdown report with a section per fragment-frame

Usage:
    python fragment_surface_analysis.py \
        --psf /path/to/system.psf \
        --dcd /path/to/trajectory.dcd \
        --pt-csv /path/to/coordination_numbers.csv \
        [--fragment-ids 0 2 5] \
        [--frame-index 0] \
        --output-dir /path/to/output_dir \
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
        description="Surface analysis for NEC fragments over one or multiple frames."
    )
    parser.add_argument("--psf", required=True, help="Path to the PSF file.")
    parser.add_argument("--dcd", required=True, help="Path to the DCD trajectory file.")
    parser.add_argument(
        "--pt-csv", required=True,
        help="CSV of Pt atom classifications (coordination numbers + is_surface)."
    )
    parser.add_argument(
        "--fragment-ids", type=int, nargs="*", default=None,
        help="IDs of fragments to analyze; if omitted, analyze all."
    )
    parser.add_argument(
        "--frame-index", type=int,
        help="Frame index to analyze; if omitted, analyze all frames."
    )
    parser.add_argument(
        "--output-dir", default=".", help="Directory to save outputs.")
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

    # Load Pt classification
    df_pt = pd.read_csv(args.pt_csv)
    pt_surface = df_pt[df_pt["is_surface"]]
    vprint(f"Loaded Pt classification: {len(pt_surface)} surface atoms")

    # Load Universe
    u = mda.Universe(args.psf, args.dcd)
    n_frames = len(u.trajectory)
    vprint(f"Loaded trajectory: {n_frames} frames")

    # Determine fragments
    nec_atoms = u.select_atoms("not name PT*")
    fragments = nec_atoms.fragments
    total_frags = len(fragments)
    fragment_ids = list(range(total_frags)) if args.fragment_ids is None else [f for f in args.fragment_ids if 0 <= f < total_frags]
    vprint(f"Analyzing fragments: {fragment_ids}")

    # Prepare output structures
    all_rows = []
    summaries = []

    # Static Pt selections
    all_pt = u.select_atoms("name PT*")
    surf_idx = pt_surface["atom_index"].values
    pt_classes = pt_surface["classification"].values

    # Frame iterator
    if args.frame_index is None:
        frames = u.trajectory
    else:
        if not (0 <= args.frame_index < n_frames):
            raise ValueError(f"frame-index {args.frame_index} out of range (0–{n_frames-1})")
        frames = [u.trajectory[args.frame_index]]
        vprint(f"Isolated frame {args.frame_index}")

    for ts in frames:
        frame = ts.frame
        vprint(f"Processing frame {frame}")

        # Build KDTree for this frame
        tree = cKDTree(all_pt.positions[surf_idx])

        for frag_id in fragment_ids:
            frag = fragments[frag_id]
            dist_arr, nn_idx = tree.query(frag.positions, k=1)
            nearest_cls = pt_classes[nn_idx]

            # Determine rows to save
            if args.all_atoms:
                for idx, d, cls in zip(frag.indices, dist_arr, nearest_cls):
                    all_rows.append({
                        "fragment_id": frag_id,
                        "frame_index": frame,
                        "atom_index": int(idx),
                        "pt_min_distance": float(d),
                        "pt_classification": cls
                    })
            else:
                i0 = np.argmin(dist_arr)
                all_rows.append({
                    "fragment_id": frag_id,
                    "frame_index": frame,
                    "atom_index": int(frag.indices[i0]),
                    "pt_min_distance": float(dist_arr[i0]),
                    "pt_classification": nearest_cls[i0]
                })

            # Summary always includes all atoms
            class_counts = pd.Series(nearest_cls).value_counts().to_dict()
            summary = {
                "fragment_id": frag_id,
                "frame_index": frame,
                "n_atoms": int(len(dist_arr)),
                "min_distance": float(dist_arr.min()),
                "max_distance": float(dist_arr.max()),
                "mean_distance": float(dist_arr.mean()),
                "median_distance": float(np.median(dist_arr)),
                "std_distance": float(dist_arr.std()),
                "pt_classification_counts": {str(k): int(v) for k, v in class_counts.items()}
            }
            summaries.append(summary)
            vprint(f"  Fragment {frag_id}: wrote {1 if not args.all_atoms else len(dist_arr)} rows, mean dist {summary['mean_distance']:.2f}")

    # Write outputs
    os.makedirs(args.output_dir, exist_ok=True)
    pd.DataFrame(all_rows).to_csv(os.path.join(args.output_dir, "surface_analysis.csv"), index=False)
    with open(os.path.join(args.output_dir, "surface_summary.json"), "w") as fj:
        json.dump(summaries, fj, indent=2)
    with open(os.path.join(args.output_dir, "surface_summary.md"), "w") as fm:
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
