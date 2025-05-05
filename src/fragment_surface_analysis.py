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
    parser.add_argument(
        "--start-frame", type=int, default=None,
        help="Optional start frame for analysis (0-based index, None for beginning)."
    )
    parser.add_argument(
        "--end-frame", type=int, default=None,
        help="Optional end frame for analysis (0-based index, None for end)."
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

    # Load Universe first to get total number of frames
    u = mda.Universe(args.psf, args.dcd)
    n_frames_full = len(u.trajectory)
    vprint(f"Loaded trajectory: {n_frames_full} frames total")

    # Determine the effective start and end frames for this run
    # These are the frames within the FULL trajectory that we will process.
    effective_start_frame = args.start_frame if args.start_frame is not None else 0
    effective_end_frame = args.end_frame if args.end_frame is not None else n_frames_full - 1

    if effective_start_frame < 0:
        vprint(f"Warning: Specified start frame ({args.start_frame}) is negative. Using 0.")
        effective_start_frame = 0

    if effective_end_frame >= n_frames_full:
         vprint(f"Warning: Specified end frame ({args.end_frame}) is beyond the last frame ({n_frames_full - 1}). Using the last frame.")
         effective_end_frame = n_frames_full - 1

    if effective_start_frame > effective_end_frame:
        print(f"Error: Effective start frame ({effective_start_frame}) is after effective end frame ({effective_end_frame}). Exiting.")
        return # Or raise an error

    vprint(f"Processing frames from {effective_start_frame} to {effective_end_frame} (inclusive).")


    # Determine starting frame for iteration based on checkpoint and specified range
    resume_from_frame = effective_start_frame # Default start frame for iteration

    if args.resume and os.path.exists(analysis_csv) and os.path.getsize(analysis_csv) > 0:
        try:
            df_checkpoint = pd.read_csv(analysis_csv)
            if not df_checkpoint.empty:
                last_frame_checkpoint = df_checkpoint['frame_index'].max()
                vprint(f"Found checkpoint file. Last processed frame in checkpoint: {last_frame_checkpoint}")

                # If the last frame in the checkpoint is within the specified range, resume from the next frame
                if effective_start_frame <= last_frame_checkpoint <= effective_end_frame:
                     resume_from_frame = last_frame_checkpoint + 1
                     vprint(f"Resuming iteration from frame {resume_from_frame} based on checkpoint.")
                elif last_frame_checkpoint < effective_start_frame:
                     # Checkpoint is before the start of the new range, start from the beginning of the new range
                     resume_from_frame = effective_start_frame
                     vprint(f"Checkpoint ({last_frame_checkpoint}) is before the specified start frame ({effective_start_frame}). Starting iteration from {resume_from_frame}.")
                else: # last_frame_checkpoint > effective_end_frame
                     # Checkpoint is after the end of the new range, the range is already processed
                     print(f"Specified frame range ({effective_start_frame}-{effective_end_frame}) has already been processed based on checkpoint ({last_frame_checkpoint}). Exiting.")
                     return # Exit as the range is already done

            else:
                vprint("Checkpoint file is empty. Starting iteration from the beginning of the specified range.")
                resume_from_frame = effective_start_frame # Start from the beginning of the specified range

        except Exception as e:
            vprint(f"Error reading checkpoint file {analysis_csv}: {e}. Starting iteration from the beginning of the specified range.")
            resume_from_frame = effective_start_frame # Start from the beginning of the specified range

    # Determine if we need to write the header
    write_header = not (args.resume and os.path.exists(analysis_csv) and os.path.getsize(analysis_csv) > 0 and resume_from_frame > effective_start_frame)
    mode = 'w' if write_header else 'a'

    if resume_from_frame > effective_end_frame:
         print(f"Resume frame ({resume_from_frame}) is beyond the effective end frame ({effective_end_frame}). Analysis for this range is complete. Exiting.")
         return # Exit if the range is already fully processed

    vprint(f"Actual iteration will start from frame {resume_from_frame}.")

    # Load Pt classification
    df_pt = pd.read_csv(args.pt_csv)
    pt_surface = df_pt[df_pt["is_surface"]]
    vprint(f"Loaded Pt classification: {len(pt_surface)} surface atoms")

    # Load Universe - now with slicing based on effective range for iteration
    # MDAnalysis slice is [start:stop:step], stop is exclusive.
    # We iterate from resume_from_frame up to and including effective_end_frame.
    # So the slice should be [resume_from_frame : effective_end_frame + 1]
    u_range = u.trajectory[resume_from_frame : effective_end_frame + 1]
    n_frames_range = len(u_range)
    vprint(f"Loaded trajectory slice for processing: {n_frames_range} frames")

    if n_frames_range == 0:
        print(f"No frames to process in the specified range ({effective_start_frame}-{effective_end_frame}) starting from resume frame ({resume_from_frame}). Exiting.")
        return # Exit if the sliced trajectory is empty

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

    # Prepare summary list (will only contain summaries for frames processed in this run)
    summaries = []

    # Open analysis CSV for writing/appending
    f_csv = open(analysis_csv, mode)
    if write_header:
        f_csv.write("fragment_id,frame_index,atom_index,pt_min_distance,pt_classification\n")
        f_csv.flush()

    # Frame iteration
    # The u_range trajectory object already represents the sliced trajectory
    # Iterating through u_range will give TimeSteps with original frame indices
    frames_to_process = u_range
    if args.frame_index is not None:
        # If a single frame is specified, override the range and process only that frame
        # Need to ensure the specified frame is within the full trajectory bounds
        if not (0 <= args.frame_index < n_frames_full):
             raise ValueError(f"frame-index {args.frame_index} out of range (0–{n_frames_full-1})")
        # Also check if the single frame is within the requested range if range is specified
        if args.start_frame is not None and args.frame_index < args.start_frame:
             print(f"Warning: Specified single frame ({args.frame_index}) is before the requested start frame ({args.start_frame}). Processing the single frame anyway.")
        if args.end_frame is not None and args.frame_index > args.end_frame:
             print(f"Warning: Specified single frame ({args.frame_index}) is after the requested end frame ({args.end_frame}). Processing the single frame anyway.")

        frames_to_process = u.trajectory[args.frame_index : args.frame_index + 1] # Slice for a single frame
        vprint(f"Isolated single frame {args.frame_index} for processing.")
        # If processing a single frame, we should probably overwrite the CSV, not append
        # unless the single frame is the *next* frame after a checkpoint.
        # Let's simplify: if frame_index is used, always overwrite for clarity.
        # This might conflict with resume, need to clarify behavior.
        # Based on the original script, frame_index overrides the normal frame iteration.
        # Let's keep that behavior but ensure it respects the resume logic if resume is also used.
        # If resume is used AND frame_index is used, and the frame_index is <= last_frame_checkpoint,
        # it means this single frame was already processed.
        # Let's refine: if frame_index is used, ignore resume and process only that frame,
        # overwriting the output file. This seems the most straightforward interpretation
        # of the original script's behavior.
        vprint("Single frame analysis requested. Ignoring resume and range parameters for frame iteration.")
        f_csv.close() # Close the potentially appended file
        mode = 'w' # Switch to write mode
        f_csv = open(analysis_csv, mode)
        f_csv.write("fragment_id,frame_index,atom_index,pt_min_distance,pt_classification\n")
        f_csv.flush()
        summaries = [] # Clear summaries for single frame run


    for ts in frames_to_process:
        frame = ts.frame # This is the original frame index
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
    # These summaries are only for the frames processed in this specific run.
    # If resuming, the summary files will only contain data from the resumed frames.
    # This might not be ideal; ideally, the summary should cover the *entire* processed range
    # if resuming. However, modifying the summary logic to read existing summaries and append
    # is more complex and not explicitly in the plan. Let's stick to the plan's implication
    # that summaries are generated from the frames processed in the current run.
    # A note in the documentation update should clarify this behavior.
    with open(summary_json, "w") as fj:
        json.dump(summaries, fj, indent=2)

    with open(summary_md, "w") as fm:
        fm.write("# Surface Analysis Report\n\n")
        # Add info about the processed range
        fm.write(f"## Analysis Range: Frames {effective_start_frame} to {effective_end_frame} (inclusive)\n\n")
        if args.resume and resume_from_frame > effective_start_frame:
             fm.write(f"*(Resumed from frame {resume_from_frame})*\n\n")

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
        # Iterate through the effective frame range
        # MDAnalysis slice is [start:stop:step], stop is exclusive.
        # We iterate from effective_start_frame up to and including effective_end_frame.
        # So the slice should be [effective_start_frame : effective_end_frame + 1]
        frames = u.trajectory[effective_start_frame : effective_end_frame + 1]
        vprint(f"Iterating through frames {effective_start_frame} to {effective_end_frame} (inclusive).")

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