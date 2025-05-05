#!/usr/bin/env python3
"""
Extract metadata from a PSF/DCD molecular dynamics simulation,
including counts of Pt atoms, NEC fragments, and timing information.

Integration timestep is reported in femtoseconds (fs); all other time metrics
are converted to nanoseconds (ns).

Outputs:
    - simulation_info.csv: table of properties and values
    - simulation_info.json: structured JSON of metadata
    - simulation_info.md: human-readable Markdown summary

Usage:
    python simulation_info.py \
        --psf /path/to/system.psf \
        --dcd /path/to/trajectory.dcd \
        --steps-per-frame 1000 \
        --output-dir /path/to/output_dir
"""
import os
import argparse
import json
import MDAnalysis as mda
import numpy as np
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

def main():
    parser = argparse.ArgumentParser(
        description="Extract simulation metadata with mixed time units."
    )
    parser.add_argument("--psf", required=True, help="Path to PSF file.")
    parser.add_argument("--dcd", required=True, help="Path to DCD trajectory.")
    parser.add_argument(
        "--steps-per-frame", type=int,
        help="Number of MD integration steps between saved frames (e.g. 1000).",
        required=True
    )
    parser.add_argument(
        "--output-dir", default=".",
        help="Directory to save outputs."
    )
    parser.add_argument(
        "--start-frame", type=int, default=None,
        help="Optional start frame for analysis (0-based index, None for beginning)."
    )
    parser.add_argument(
        "--end-frame", type=int, default=None,
        help="Optional end frame for analysis (0-based index, None for end)."
    )
    parser.add_argument(
        "--timestep-fs", type=float, default=None,
        help="Optional override for integration timestep in femtoseconds."
    )
    args = parser.parse_args()

    # Load Universe
    u = mda.Universe(args.psf, args.dcd)

    # --- Full Simulation Info ---
    # Atom counts
    n_total_atoms_full = int(u.atoms.n_atoms)
    n_pt_atoms_full = int(len(u.select_atoms("name PT*")))
    n_nec_fragments_full = int(len(u.select_atoms("not name PT*").fragments))

    # Frames
    n_frames_full = int(len(u.trajectory))

    # Determine integration timestep and source
    if args.timestep_fs is not None:
        dt_fs = args.timestep_fs
        timestep_source = "User argument"
    else:
        # Integration timestep from DCD header (AKMA units == fs)
        dt_fs = float(u.trajectory.dt)
        timestep_source = "DCD header"

    # Steps per frame from user
    spf = args.steps_per_frame

    # Time between frames in fs, then convert to ns
    dt_frame_fs = dt_fs * spf
    dt_frame_ns = dt_frame_fs / 1e6

    # Timeline in ns (start at zero) for FULL simulation
    start_ns_full = 0.0
    end_ns_full = dt_frame_ns * (n_frames_full - 1)
    total_ns_full = end_ns_full - start_ns_full

    # Box dimensions (Å) from first frame of FULL simulation
    try:
        dims_raw_full = u.trajectory.ts.dimensions[:3]
        box_lengths_full = [float(x) for x in dims_raw_full]
    except Exception:
        box_lengths_full = None

    # Atom type breakdown for FULL simulation
    unique_full, counts_full = np.unique(u.atoms.names, return_counts=True)
    atom_type_counts_full = {str(name): int(cnt) for name, cnt in zip(unique_full.tolist(), counts_full.tolist())}

    # --- Processed Range Info ---
    start_frame = args.start_frame
    end_frame = args.end_frame

    # Slice trajectory based on start and end frames
    # MDAnalysis slicing is [start:stop:step], stop is exclusive.
    # If end_frame is None, slice goes to the end.
    # If start_frame is None, slice starts from the beginning.
    # If both are None, it's the full trajectory.
    # We need to handle the case where end_frame is specified but is the last frame.
    # MDAnalysis slice [start:end] includes 'start' but excludes 'end'.
    # To include the frame specified by end_frame, we need to slice up to end_frame + 1.
    # If end_frame is not None:
    #     # Ensure end_frame is within bounds of the full trajectory
    #     if end_frame >= n_frames_full:
    #         print(f"Warning: Specified end frame ({end_frame}) is beyond the last frame of the trajectory ({n_frames_full - 1}). Using the last frame.")
    #         end_frame_slice = n_frames_full
    #     else:
    #         end_frame_slice = end_frame + 1
    # else:
    #     end_frame_slice = None # Slice to the end

    # if start_frame is not None and start_frame < 0:
    #      print(f"Warning: Specified start frame ({start_frame}) is negative. Using the first frame (0).")
    #      start_frame_slice = 0
    # else:
    #      start_frame_slice = start_frame # None or positive index

    # Apply slicing - MDAnalysis handles None for start/end correctly
    # and handles end index being inclusive when slicing trajectory objects directly.
    # We need to be careful with the end_frame + 1 logic if we were to use
    # u.trajectory[start:end+1] but the plan suggests using the arguments directly
    # with MDAnalysis slicing which is [start:stop:step] where stop is exclusive.
    # Let's re-read the MDAnalysis slicing documentation or test this.
    # For now, let's assume the original slicing logic was correct for the desired range.
    # Reverting to the original slicing logic based on the previous code state.

    if end_frame is not None:
        # Ensure end_frame is within bounds of the full trajectory
        if end_frame >= n_frames_full:
            print(f"Warning: Specified end frame ({end_frame}) is beyond the last frame of the trajectory ({n_frames_full - 1}). Using the last frame.")
            end_frame_slice = n_frames_full
        else:
            end_frame_slice = end_frame + 1
    else:
        end_frame_slice = None # Slice to the end

    if start_frame is not None and start_frame < 0:
         print(f"Warning: Specified start frame ({start_frame}) is negative. Using the first frame (0).")
         start_frame_slice = 0
    else:
         start_frame_slice = start_frame # None or positive index


    u_range = u.trajectory[start_frame_slice:end_frame_slice]


    # Frames in processed range
    n_frames_range = int(len(u_range))

    # Timeline in ns for PROCESSED RANGE
    # Need to calculate the actual start and end frame indices within the full trajectory
    # that correspond to the start and end of the sliced trajectory.
    # If start_frame was None, the actual start index is 0. Otherwise, it's start_frame.
    actual_start_frame_index = start_frame if start_frame is not None else 0
    # If end_frame was None, the actual end index is n_frames_full - 1. Otherwise, it's end_frame.
    # If end_frame was specified, the actual end index is end_frame.
    actual_end_frame_index = end_frame if end_frame is not None and end_frame < n_frames_full else n_frames_full - 1


    start_ns_range = dt_frame_ns * actual_start_frame_index
    end_ns_range = dt_frame_ns * actual_end_frame_index
    total_ns_range = end_ns_range - start_ns_range

    # --- Assemble Output Data ---

    # Data for JSON output
    info_full = {
        "psf_file": args.psf,
        "dcd_file": args.dcd,
        "n_total_atoms": n_total_atoms_full,
        "n_pt_atoms": n_pt_atoms_full,
        "n_nec_fragments": n_nec_fragments_full,
        "n_frames": n_frames_full,
        "integration_timestep_fs": dt_fs,
        "timestep_source": timestep_source, # Source of the timestep
        "steps_per_frame": spf,
        "time_between_frames_ns": dt_frame_ns,
        "start_time_ns": start_ns_full,
        "end_time_ns": end_ns_full,
        "total_time_ns": total_ns_full,
        "box_lengths": box_lengths_full,
        "atom_type_counts": atom_type_counts_full
    }

    info_range = {
        "start_frame_index": start_frame, # Original requested start frame
        "end_frame_index": end_frame,     # Original requested end frame
        "n_frames": n_frames_range,
        "start_time_ns": start_ns_range,
        "end_time_ns": end_ns_range,
        "total_time_ns": total_ns_range,
        # Include time_between_frames_ns in range info for convenience
        "time_between_frames_ns": dt_frame_ns,
        "integration_timestep_fs": dt_fs, # Include timestep in range info
        "timestep_source": timestep_source, # Include timestep source in range info
        "steps_per_frame": spf, # Include steps per frame in range info
    }

    # Combined JSON structure
    info_combined = {
        "full_simulation": info_full,
        "processed_range": info_range
    }

    # Data for CSV output
    import pandas as pd

    # Full simulation data
    rows_full = [
        ("scope", "full"),
        ("psf_file", args.psf),
        ("dcd_file", args.dcd),
        ("n_total_atoms", n_total_atoms_full),
        ("n_pt_atoms", n_pt_atoms_full),
        ("n_nec_fragments", n_nec_fragments_full),
        ("n_frames", n_frames_full),
        ("integration_timestep_fs", dt_fs),
        ("timestep_source", timestep_source), # Add source to CSV
        ("steps_per_frame", spf),
        ("time_between_frames_ns", dt_frame_ns),
        ("start_time_ns", start_ns_full),
        ("end_time_ns", end_ns_full),
        ("total_time_ns", total_ns_full),
        ("box_lengths", box_lengths_full),
        ("atom_type_counts", atom_type_counts_full)
    ]

    # Processed range data
    rows_range = [
        ("scope", "range"),
        ("start_frame_index", start_frame),
        ("end_frame_index", end_frame),
        ("n_frames", n_frames_range),
        ("time_between_frames_ns", dt_frame_ns), # Include for consistency
        ("start_time_ns", start_ns_range),
        ("end_time_ns", end_ns_range),
        ("total_time_ns", total_ns_range),
        ("integration_timestep_fs", dt_fs), # Add timestep to range CSV
        ("timestep_source", timestep_source), # Add source to range CSV
        ("steps_per_frame", spf), # Add steps per frame to range CSV
        # Other metrics like atom counts, box lengths are for the full simulation
    ]

    # Combine rows and create DataFrame
    rows_combined = rows_full + rows_range
    df = pd.DataFrame(rows_combined, columns=["property", "value"])

    # Ensure output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Write simulation_info.csv
    csv_path = os.path.join(args.output_dir, "simulation_info.csv")
    df.to_csv(csv_path, index=False)
    print(f"Saved CSV: {csv_path}")

    # Ensure output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Write simulation_info.csv
    csv_path = os.path.join(args.output_dir, "simulation_info.csv")
    df.to_csv(csv_path, index=False)
    print(f"Saved CSV: {csv_path}")


    # Write simulation_info.json
    json_path = os.path.join(args.output_dir, "simulation_info.json")
    with open(json_path, 'w') as jf:
        json.dump(info_combined, jf, indent=2)
    print(f"Saved JSON: {json_path}")

    # Write simulation_info.md
    md_path = os.path.join(args.output_dir, "simulation_info.md")
    with open(md_path, 'w') as mf:
        mf.write("# Simulation Metadata Summary\n\n")

        mf.write("## Full Simulation Information\n\n")
        mf.write(f"- PSF file: `{args.psf}`\n")
        mf.write(f"- DCD file: `{args.dcd}`\n")
        mf.write(f"- Total atoms: **{info_full['n_total_atoms']}**\n")
        mf.write(f"- Pt atoms: **{info_full['n_pt_atoms']}**\n")
        mf.write(f"- NEC fragments: **{info_full['n_nec_fragments']}**\n")
        mf.write(f"- Number of frames: **{info_full['n_frames']}**\n")
        mf.write(f"- Integration timestep: **{info_full['integration_timestep_fs']:.6f}** fs (Source: {info_full['timestep_source']})\n")
        mf.write(f"- Steps per frame: **{info_full['steps_per_frame']}**\n")
        mf.write(f"- Time between frames: **{info_full['time_between_frames_ns']:.6f}** ns\n")
        mf.write(f"- Start time: **{info_full['start_time_ns']:.6f}** ns\n")
        mf.write(f"- End time: **{info_full['end_time_ns']:.6f}** ns\n")
        mf.write(f"- Total simulation time: **{info_full['total_time_ns']:.6f}** ns\n")
        if info_full['box_lengths'] is not None:
            mf.write(f"- Box lengths (Å): **{info_full['box_lengths'][0]:.2f}, {info_full['box_lengths'][1]:.2f}, {info_full['box_lengths'][2]:.2f}**\n")
        mf.write("- Atom counts by type:\n")
        for name, cnt in info_full['atom_type_counts'].items():
            mf.write(f"  - {name}: {cnt}\n")

        mf.write("\n## Processed Frame Range Information\n\n")
        mf.write(f"- Requested start frame (0-based): **{info_range['start_frame_index']}**\n")
        mf.write(f"- Requested end frame (0-based): **{info_range['end_frame_index']}**\n")
        mf.write(f"- Number of frames in range: **{info_range['n_frames']}**\n")
        mf.write(f"- Time between frames: **{info_range['time_between_frames_ns']:.6f}** ns\n")
        mf.write(f"- Start time of range: **{info_range['start_time_ns']:.6f}** ns\n")
        mf.write(f"- End time of range: **{info_range['end_time_ns']:.6f}** ns\n")
        mf.write(f"- Total time in range: **{info_range['total_time_ns']:.6f}** ns\n")
        mf.write(f"- Integration timestep: **{info_range['integration_timestep_fs']:.6f}** fs (Source: {info_range['timestep_source']})\n") # Added source info for range
        mf.write(f"- Steps per frame: **{info_range['steps_per_frame']}**\n") # Added steps per frame for range


    print(f"Saved Markdown: {md_path}")

if __name__ == "__main__":
    main()
