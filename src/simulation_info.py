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
    args = parser.parse_args()

    # Load Universe
    u = mda.Universe(args.psf, args.dcd)

    # Atom counts
    n_total_atoms = int(u.atoms.n_atoms)
    n_pt_atoms = int(len(u.select_atoms("name PT*")))
    n_nec_fragments = int(len(u.select_atoms("not name PT*").fragments))

    # Frames
    n_frames = int(len(u.trajectory))

    # Integration timestep from DCD header (AKMA units == fs)
    dt_fs = float(u.trajectory.dt)

    # Steps per frame from user
    spf = args.steps_per_frame

    # Time between frames in fs, then convert to ns
    dt_frame_fs = dt_fs * spf
    dt_frame_ns = dt_frame_fs / 1e6

    # Timeline in ns (start at zero)
    start_ns = 0.0
    end_ns = dt_frame_ns * (n_frames - 1)
    total_ns = end_ns - start_ns

    # Box dimensions (Å) from first frame
    try:
        dims_raw = u.trajectory.ts.dimensions[:3]
        box_lengths = [float(x) for x in dims_raw]
    except Exception:
        box_lengths = None

    # Atom type breakdown
    unique, counts = np.unique(u.atoms.names, return_counts=True)
    atom_type_counts = {str(name): int(cnt) for name, cnt in zip(unique.tolist(), counts.tolist())}

    # Assemble CSV rows
    rows = [
        ("psf_file", args.psf),
        ("dcd_file", args.dcd),
        ("n_total_atoms", n_total_atoms),
        ("n_pt_atoms", n_pt_atoms),
        ("n_nec_fragments", n_nec_fragments),
        ("n_frames", n_frames),
        ("integration_timestep_fs", dt_fs),
        ("steps_per_frame", spf),
        ("time_between_frames_ns", dt_frame_ns),
        ("start_time_ns", start_ns),
        ("end_time_ns", end_ns),
        ("total_time_ns", total_ns),
        ("box_lengths", box_lengths),
        ("atom_type_counts", atom_type_counts)
    ]
    import pandas as pd
    df = pd.DataFrame(rows, columns=["property", "value"])

    # Ensure output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Write simulation_info.csv
    csv_path = os.path.join(args.output_dir, "simulation_info.csv")
    df.to_csv(csv_path, index=False)
    print(f"Saved CSV: {csv_path}")

    # JSON
    info = {
        "psf_file": args.psf,
        "dcd_file": args.dcd,
        "n_total_atoms": n_total_atoms,
        "n_pt_atoms": n_pt_atoms,
        "n_nec_fragments": n_nec_fragments,
        "n_frames": n_frames,
        "integration_timestep_fs": dt_fs,
        "steps_per_frame": spf,
        "time_between_frames_ns": dt_frame_ns,
        "start_time_ns": start_ns,
        "end_time_ns": end_ns,
        "total_time_ns": total_ns,
        "box_lengths": box_lengths,
        "atom_type_counts": atom_type_counts
    }
    json_path = os.path.join(args.output_dir, "simulation_info.json")
    with open(json_path, 'w') as jf:
        json.dump(info, jf, indent=2)
    print(f"Saved JSON: {json_path}")

    # Markdown
    md_path = os.path.join(args.output_dir, "simulation_info.md")
    with open(md_path, 'w') as mf:
        mf.write("# Simulation Metadata Summary\n\n")
        mf.write(f"- PSF file: `{args.psf}`\n")
        mf.write(f"- DCD file: `{args.dcd}`\n")
        mf.write(f"- Total atoms: **{n_total_atoms}**\n")
        mf.write(f"- Pt atoms: **{n_pt_atoms}**\n")
        mf.write(f"- NEC fragments: **{n_nec_fragments}**\n")
        mf.write(f"- Number of frames: **{n_frames}**\n")
        mf.write(f"- Integration timestep: **{dt_fs:.6f}** fs\n")
        mf.write(f"- Steps per frame: **{spf}**\n")
        mf.write(f"- Time between frames: **{dt_frame_ns:.6f}** ns\n")
        mf.write(f"- Start time: **{start_ns:.6f}** ns\n")
        mf.write(f"- End time: **{end_ns:.6f}** ns\n")
        mf.write(f"- Total simulation time: **{total_ns:.6f}** ns\n")
        if box_lengths is not None:
            mf.write(f"- Box lengths (Å): **{box_lengths[0]:.2f}, {box_lengths[1]:.2f}, {box_lengths[2]:.2f}**\n")
        mf.write("- Atom counts by type:\n")
        for name, cnt in atom_type_counts.items():
            mf.write(f"  - {name}: {cnt}\n")
    print(f"Saved Markdown: {md_path}")

if __name__ == "__main__":
    main()
