#!/usr/bin/env python3
"""
Classify Pt atoms based on coordination numbers derived from RDF analysis or a user-supplied cutoff.

This script computes an average Pt–Pt RDF on a specified frame (or range),
finds the first minimum after the main peak to define a cutoff,
then counts neighbors within that cutoff via a KDTree to assign coordination numbers.
Atoms are classified into categories like "bulk", "facet", "edge", and "vertex".

Usage:
    python classify_pt_atoms.py \\
        --psf /path/to/system.psf \\
        --dcd /path/to/trajectory.dcd \\
        --output-dir /path/to/output_dir \\
        [--prefix system_prefix] \\
        [--frame-index INDEX] \\
        [--use-cutoff CUT_OFF] \\
        [--r-max 10.0] [--nbins 200] \\
        [--min-frame INDEX] [--max-frame INDEX] \\
        [--plot-rdf]

Outputs:
    - {prefix}_coordination_numbers.csv: atom_index, coordination_number, classification, is_surface
    - {prefix}_coordination_numbers.txt: groupings by classification
    - {prefix}_classification_summary.md: Markdown summary of classification metrics
    - (optional) {prefix}_rdf_curve.png: visualization of g(r) and chosen cutoff
"""
import os
import argparse
import json
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF
from scipy.signal import argrelextrema
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

# Default thresholds for classification
def classify_atom(cn):
    if cn >= 11:
        return 'bulk'
    elif cn == 9:
        return '111 facet'
    elif cn == 8:
        return '100 facet'
    elif 6 <= cn < 8:
        return 'edge'
    else:
        return 'vertex'


def main():
    parser = argparse.ArgumentParser(description='Classify Pt atoms by coordination number.')
    parser.add_argument('--psf', required=True, help='PSF file')
    parser.add_argument('--dcd', required=True, help='DCD file')
    parser.add_argument('--output-dir', required=True, help='Directory for output files')
    parser.add_argument('--prefix', default='system', help='Prefix for output names')
    parser.add_argument('--frame-index', type=int,
                        help='Single frame to analyze (default: middle frame)')
    parser.add_argument('--min-frame', type=int,
                        help='Start frame for RDF averaging')
    parser.add_argument('--max-frame', type=int,
                        help='End frame for RDF averaging')
    parser.add_argument('--use-cutoff', type=float,
                        help='Bypass RDF: use this cutoff (Å) directly')
    parser.add_argument('--r-max', type=float, default=10.0,
                        help='Max distance (Å) for RDF')
    parser.add_argument('--nbins', type=int, default=200,
                        help='Number of bins for RDF')
    parser.add_argument('--plot-rdf', action='store_true',
                        help='Save RDF curve plot with cutoff')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load Universe and determine frames
    u = mda.Universe(args.psf, args.dcd)
    n_frames = len(u.trajectory)
    if args.min_frame is not None or args.max_frame is not None:
        start = args.min_frame or 0
        stop = (args.max_frame + 1) if args.max_frame is not None else n_frames
        frames = list(range(start, min(stop, n_frames)))
    else:
        idx = args.frame_index if args.frame_index is not None else n_frames // 2
        frames = [idx]

    # Select Pt atoms
    pt = u.select_atoms('name PT*')
    n_pt = len(pt)
    print(f'Classifying {n_pt} Pt atoms over frames {frames}')

    # Determine cutoff
    if args.use_cutoff is not None:
        dist_cut = args.use_cutoff
        print(f'Using user-supplied cutoff: {dist_cut:.2f} Å')
    else:
        # Compute RDF over selected frames
        rdf = InterRDF(pt, pt, nbins=args.nbins, range=(0.0, args.r_max))
        rdf.run(start=frames[0], stop=frames[-1]+1)
        g_r = rdf.rdf
        bins = rdf.bins
        # locate first peak then first minimum
        peaks = argrelextrema(g_r, np.greater)[0]
        if len(peaks) == 0:
            raise RuntimeError('No RDF peak found; cannot auto-select cutoff')
        peak = peaks[0]
        minima = argrelextrema(g_r, np.less)[0]
        dist_cut = None
        for mi in minima:
            if mi > peak:
                dist_cut = bins[mi]
                break
        if dist_cut is None:
            raise RuntimeError('No RDF minimum after first peak; increase r-max or nbins')
        print(f'Auto-detected cutoff: {dist_cut:.2f} Å')
        if args.plot_rdf:
            fig, ax = plt.subplots()
            ax.plot(bins, g_r, label='g(r)')
            ax.axvline(dist_cut, color='k', linestyle='--', label=f'cutoff={dist_cut:.2f} Å')
            ax.set_xlabel('r (Å)')
            ax.set_ylabel('g(r)')
            ax.legend()
            fig.savefig(os.path.join(args.output_dir, f'{args.prefix}_rdf_curve.png'), dpi=200)
            plt.close(fig)

    # KDTree and coordination counting
    coords = pt.positions
    tree = cKDTree(coords)
    neighbor_lists = tree.query_ball_point(coords, dist_cut)
    cn = np.array([len(lst)-1 for lst in neighbor_lists], dtype=int)

    # Classification
    classes = [classify_atom(x) for x in cn]
    df = pd.DataFrame({
        'atom_index': np.arange(n_pt),
        'coordination_number': cn,
        'classification': classes
    })
    df['is_surface'] = df['classification'] != 'bulk'

    # --- Generate Classification Summary ---
    total_pt_atoms = n_pt
    surface_atoms = df[df['is_surface']].shape[0]
    classification_counts = df['classification'].value_counts().to_dict()

    summary_md = f"""# Platinum Atom Classification Summary

*   Total Pt Atoms: **{total_pt_atoms}**
*   Total Surface Atoms: **{surface_atoms}**
*   Atoms per Category:
"""
    for category, count in classification_counts.items():
        summary_md += f"    *   {category}: **{count}**\n"

    summary_md_file = os.path.join(args.output_dir, f'{args.prefix}_classification_summary.md')
    with open(summary_md_file, 'w') as f:
        f.write(summary_md)
    print(f'Wrote Markdown summary: {summary_md_file}')

    # Save CSV
    csv_file = os.path.join(args.output_dir, f'{args.prefix}_coordination_numbers.csv')
    df.to_csv(csv_file, index=False)
    print(f'Wrote CSV: {csv_file}')

    # Save text grouping
    txt_file = os.path.join(args.output_dir, f'{args.prefix}_coordination_numbers.txt')
    with open(txt_file, 'w') as f:
        for grp, sub in df.groupby('classification'):
            f.write(f'{grp}\n')
            f.write(' '.join(map(str, sub['atom_index'].tolist())) + '\n\n')
        f.write('surface\n')
        f.write(' '.join(map(str, df[df['is_surface']]['atom_index'].tolist())))
    print(f'Wrote text file: {txt_file}')


if __name__ == '__main__':
    main()
