#!/usr/bin/env python3
"""
Quantify residence events of NEC fragments on the Pt nanoparticle surface.

A "residence event" is defined as:
  1) A contiguous period where a fragment is within the spatial cutoff (as in surface_filtered.csv)
     for at least `min_res_time` (ns).
  2) Gaps of off-surface shorter than or equal to `max_off_time` (ns) merge two events.

Reads:
  - surface_filtered.csv: per-atom, per-frame entries for fragments within cutoff
  - simulation_info.json: JSON with simulation metadata including time_between_frames_ns, start_time_ns, n_frames

Writes to the output directory:
  - residence_events.csv: one row per event (fragment_id, event_id, start_frame, end_frame, start_time_ns, end_time_ns, duration_ns)
  - residence_events_summary.json: summary per fragment with event list and stats
  - residence_events_summary.md: Markdown report of events and summaries

Usage:
  python residence_event_analysis.py \
      --filtered-csv /path/to/surface_filtered.csv \
      --sim-info-json /path/to/simulation_info.json \
      [--min-res-time 0.1] \
      [--max-off-time 0.1] \
      --output-dir /path/to/output_dir
"""
import os
import argparse
import json
import pandas as pd
import numpy as np
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

def main():
    parser = argparse.ArgumentParser(
        description="Quantify residence events per fragment from surface_filtered data.")
    parser.add_argument(
        "--filtered-csv", required=True,
        help="Path to surface_filtered.csv (from filter_surface_fragments.py)"
    )
    parser.add_argument(
        "--sim-info-json", required=True,
        help="Path to simulation_info.json with time info (ns)."
    )
    parser.add_argument(
        "--min-res-time", type=float, default=0.1,
        help="Minimum contiguous residence time (ns) to count as an event."
    )
    parser.add_argument(
        "--max-off-time", type=float, default=0.1,
        help="Maximum off-surface gap (ns) to merge adjacent events."
    )
    parser.add_argument(
        "--output-dir", default=".", help="Directory to save outputs."
    )
    args = parser.parse_args()

    # Load filtered data
    df = pd.read_csv(args.filtered_csv)
    # Load sim info JSON
    with open(args.sim_info_json) as f:
        sim = json.load(f)
    # Build frame->time mapping in ns
    n_frames = int(sim.get("n_frames"))
    start_ns = float(sim.get("start_time_ns"))
    dt_ns = float(sim.get("time_between_frames_ns"))
    times = start_ns + np.arange(n_frames) * dt_ns

    # Identify fragments (cast to Python ints to avoid JSON serialization errors)
    fragment_ids = sorted(df['fragment_id'].unique())
    fragment_ids = [int(x) for x in fragment_ids]

    # Container for all events
    all_events = []
    summaries = []

    # Iterate per fragment
    for frag in fragment_ids:
        # Boolean on/off series over all frames
        on_frames = set(df.loc[df['fragment_id']==frag, 'frame_index'])
        on_series = np.array([i in on_frames for i in range(n_frames)], dtype=bool)

        # Detect raw runs
        runs = []  # list of (start_idx, end_idx)
        in_run = False
        for i, on in enumerate(on_series):
            if on and not in_run:
                start = i
                in_run = True
            elif not on and in_run:
                end = i-1
                runs.append((start, end))
                in_run = False
        # handle run to final frame
        if in_run:
            runs.append((start, n_frames-1))

        # Compute run durations and filter by min_res_time
        filtered_runs = []
        for (s, e) in runs:
            t0 = times[s]
            t1 = times[e]
            dur = t1 - t0
            if dur >= args.min_res_time:
                filtered_runs.append({'start': s, 'end': e, 'start_ns': t0, 'end_ns': t1, 'dur_ns': dur})

        # Merge runs separated by small off-gaps
        merged = []
        for run in filtered_runs:
            if not merged:
                merged.append(run.copy())
            else:
                prev = merged[-1]
                gap = run['start_ns'] - prev['end_ns']
                if gap <= args.max_off_time:
                    # extend previous event
                    prev['end'] = run['end']
                    prev['end_ns'] = run['end_ns']
                    prev['dur_ns'] = prev['end_ns'] - prev['start_ns']
                else:
                    merged.append(run.copy())

        # Record events for this fragment
        frag_events = []
        for eid, run in enumerate(merged):
            evt = {
                'fragment_id': frag,
                'event_id': eid,
                'start_frame': int(run['start']),
                'end_frame': int(run['end']),
                'start_time_ns': float(run['start_ns']),
                'end_time_ns': float(run['end_ns']),
                'duration_ns': float(run['dur_ns'])
            }
            all_events.append(evt)
            frag_events.append(evt)

        # Summarize for this fragment
        durations = [e['duration_ns'] for e in frag_events]
        if durations:
            summaries.append({
                'fragment_id': frag,
                'n_events': len(durations),
                'total_residence_ns': float(sum(durations)),
                'mean_event_ns': float(np.mean(durations)),
                'events': [
                    {'event_id': e['event_id'],
                     'start_time_ns': e['start_time_ns'],
                     'end_time_ns': e['end_time_ns'],
                     'duration_ns': e['duration_ns']}
                    for e in frag_events
                ]
            })
        else:
            summaries.append({
                'fragment_id': frag,
                'n_events': 0,
                'total_residence_ns': 0.0,
                'mean_event_ns': 0.0,
                'events': []
            })

    # Ensure output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Write CSV of events
    df_evt = pd.DataFrame(all_events)
    csv_out = os.path.join(args.output_dir, 'residence_events.csv')
    df_evt.to_csv(csv_out, index=False)
    print(f"Saved events CSV: {csv_out}")

    # Write JSON summary
    json_out = os.path.join(args.output_dir, 'residence_events_summary.json')
    with open(json_out, 'w') as fj:
        json.dump(summaries, fj, indent=2)
    print(f"Saved summary JSON: {json_out}")

    # Write Markdown report
    md_out = os.path.join(args.output_dir, 'residence_events_summary.md')
    with open(md_out, 'w') as fm:
        fm.write("# Residence Event Analysis\n\n")
        fm.write(f"Min residence time: **{args.min_res_time} ns**, ")
        fm.write(f"Max off-gap: **{args.max_off_time} ns**\n\n")
        for s in summaries:
            fm.write(f"## Fragment {s['fragment_id']}\n")
            fm.write(f"- Number of events: **{s['n_events']}**\n")
            fm.write(f"- Total residence: **{s['total_residence_ns']:.3f} ns**\n")
            fm.write(f"- Mean event duration: **{s['mean_event_ns']:.3f} ns**\n")
            if s['events']:
                fm.write("### Events:\n")
                for e in s['events']:
                    fm.write(
                        f"- Event {e['event_id']}: {e['duration_ns']:.3f} ns ("
                        f"{e['start_time_ns']:.3f} → {e['end_time_ns']:.3f} ns)\n"
                    )
            fm.write("\n")
    print(f"Saved Markdown report: {md_out}")

if __name__ == '__main__':
    main()
