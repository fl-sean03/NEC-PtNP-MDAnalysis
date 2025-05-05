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
import sys # Import sys for stderr printing

warnings.filterwarnings("ignore", category=DeprecationWarning)

def save_checkpoint(file_path, unprocessed_fragment_ids, processed_events, processed_summaries):
    """Saves the current state to a checkpoint file."""
    checkpoint_data = {
        "unprocessed_fragment_ids": unprocessed_fragment_ids,
        "processed_events": processed_events,
        "processed_summaries": processed_summaries
    }
    try:
        # Ensure directory exists before saving
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        with open(file_path, 'w') as f:
            json.dump(checkpoint_data, f, indent=2)
    except Exception as e:
        print(f"Error saving checkpoint file {file_path}: {e}", file=sys.stderr)
        # In a more robust system, consider logging or alternative saving

def load_checkpoint(file_path):
    """Loads state from a checkpoint file."""
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return None # No checkpoint or empty file

    try:
        with open(file_path, 'r') as f:
            checkpoint_data = json.load(f)
            # Basic validation
            if all(k in checkpoint_data for k in ["unprocessed_fragment_ids", "processed_events", "processed_summaries"]):
                return checkpoint_data
            else:
                print(f"Checkpoint file {file_path} is missing required keys. Starting from scratch.", file=sys.stderr)
                return None
    except json.JSONDecodeError:
        print(f"Checkpoint file {file_path} is corrupted (invalid JSON). Starting from scratch.", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error loading checkpoint file {file_path}: {e}. Starting from scratch.", file=sys.stderr)
        return None

# Define checkpoint save interval (e.g., every 10 fragments)
SAVE_INTERVAL = 10

def process_fragment(frag, df, times, args, n_frames):
    """Processes a single fragment to find residence events."""
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
        frag_events.append(evt)

    # Summarize for this fragment
    durations = [e['duration_ns'] for e in frag_events]
    if durations:
        frag_summary = {
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
        }
    else:
        frag_summary = {
            'fragment_id': frag,
            'n_events': 0,
            'total_residence_ns': 0.0,
            'mean_event_ns': 0.0,
            'events': []
        }

    return frag_events, frag_summary


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
    parser.add_argument(
        "--checkpoint-file", default=None,
        help="Path to the checkpoint file. Defaults to [output_dir]/residence_events_checkpoint.json"
    )
    args = parser.parse_args()

    # Define checkpoint file path, using default if not provided
    checkpoint_file = args.checkpoint_file if args.checkpoint_file else os.path.join(args.output_dir, "residence_events_checkpoint.json")

    # Initialize state variables
    unprocessed_fragment_ids = []
    processed_events = []
    processed_summaries = []

    # Attempt to load checkpoint
    checkpoint_data = load_checkpoint(checkpoint_file)

    if checkpoint_data:
        unprocessed_fragment_ids = checkpoint_data["unprocessed_fragment_ids"]
        processed_events = checkpoint_data["processed_events"]
        processed_summaries = checkpoint_data["processed_summaries"]
        print(f"Resuming from checkpoint: {checkpoint_file}")
    else:
        print("No valid checkpoint found. Starting from scratch.")
        # Load filtered data to get all fragment IDs if starting from scratch
        df_initial = pd.read_csv(args.filtered_csv)
        original_fragment_ids = sorted(df_initial['fragment_id'].unique())
        unprocessed_fragment_ids = [int(x) for x in original_fragment_ids]
        processed_events = []
        processed_summaries = []

    # Load simulation info (needed regardless of checkpoint)
    with open(args.sim_info_json) as f:
        sim_info_data = json.load(f)

    # Get time information for the processed range
    # The simulation_info.json now contains 'full_simulation' and 'processed_range' keys
    if "processed_range" not in sim_info_data:
        raise KeyError("simulation_info.json does not contain the 'processed_range' key. Ensure simulation_info.py was run with a recent version.")

    sim_range = sim_info_data["processed_range"]

    # Build frame->time mapping in ns for the PROCESSED RANGE
    # Note: n_frames here refers to the number of frames *in the processed range*
    # The frame indices in surface_filtered.csv are still the original indices.
    # We need to adjust the time mapping to align with the original frame indices
    # but use the time information specific to the range.

    n_frames_range = int(sim_range.get("n_frames"))
    start_ns_range = float(sim_range.get("start_time_ns"))
    dt_ns = float(sim_range.get("time_between_frames_ns"))
    start_frame_index_range = int(sim_range.get("start_frame_index", 0)) # Get the original start frame index of the range

    # The times array should map original frame indices within the range to their times.
    # The frames in df['frame_index'] are the original frame indices.
    # We need a time array that covers the full original range of frames,
    # but with times calculated based on the range's start time and dt.
    # Let's get the full simulation info to know the total number of frames for the time array size.
    if "full_simulation" not in sim_info_data:
         raise KeyError("simulation_info.json does not contain the 'full_simulation' key.")
    sim_full = sim_info_data["full_simulation"]
    n_frames_full = int(sim_full.get("n_frames"))

    # Create a time array for the full simulation length, but calculate times based on the range's dt
    # The time for frame `i` (original index) is `start_ns_range + (i - start_frame_index_range) * dt_ns`
    # This assumes the time_between_frames_ns is constant across the full trajectory.
    # This approach ensures that the times array aligns with the original frame indices in surface_filtered.csv.
    times = np.zeros(n_frames_full)
    # Only calculate times for the frames within the processed range to avoid issues with indices outside the range
    # The frames in df['frame_index'] will be within the original range [start_frame_index_range, start_frame_index_range + n_frames_range - 1]
    # Let's create the times array only for the frames present in the filtered CSV to be safe.
    # Or, create the full time array but ensure calculations use the correct start time and dt.

    # Let's create the time array for the full simulation length, but calculate times relative to the start of the range.
    # This way, the index `i` in the `times` array corresponds to the original frame index `i`.
    # The time for original frame index `i` is `start_ns_range + (i - start_frame_index_range) * dt_ns`.
    # This is only valid for frames within the processed range.
    # We need to ensure that when process_fragment is called, it only considers frames within the processed range.
    # The df passed to process_fragment should already be filtered by the pipeline based on the range.
    # However, the `on_series` is created based on `range(n_frames)`, where `n_frames` is currently the full number of frames.
    # This needs to be `n_frames_full` to align with the original frame indices in `df`.
    # The `times` array should also be indexed by original frame index.

    # Let's recalculate the `times` array to cover the full range of original frame indices,
    # but with times relative to the start of the processed range.
    # The time for original frame index `i` is `start_ns_range + (i - start_frame_index_range) * dt_ns`.
    # This is only meaningful for frames `i` where `start_frame_index_range <= i < start_frame_index_range + n_frames_range`.
    # For frames outside this range, the time calculation might be incorrect, but those frames shouldn't be in df anyway.

    # Let's create the times array for the full simulation length, indexed by original frame index.
    times = np.array([start_ns_range + (i - start_frame_index_range) * dt_ns for i in range(n_frames_full)])

    # The process_fragment function uses `range(n_frames)` to create `on_series`.
    # This `n_frames` should be `n_frames_full` to match the original frame indices in `df`.
    # The `times` array is indexed by original frame index.
    # The `df` contains original frame indices.
    # This seems consistent. The key is that `times` is indexed by original frame index,
    # and the times themselves are calculated based on the range's start time and dt.
    # The `on_series` should be based on the full range of original frame indices to correctly
    # represent the presence/absence of the fragment across the *entire* original trajectory,
    # even though `df` only contains data for frames within the processed range.
    # This is because the gap merging logic needs to consider frames outside the processed range
    # if they are within the `max_off_time` gap.

    # Let's keep n_frames as n_frames_full for the on_series creation.
    n_frames = n_frames_full # Use full number of frames for on_series indexing

    # The times array is correctly indexed by original frame index.
    # The times themselves are relative to the start of the processed range.
    # This seems correct for calculating durations and event times within the processed range.

    # Load filtered data (needed regardless of checkpoint for processing)
    df = pd.read_csv(args.filtered_csv)

    # Determine fragments to process in this run
    # Create a copy as we will modify unprocessed_fragment_ids during iteration
    fragments_to_process_this_run = list(unprocessed_fragment_ids)

    # Filter the dataframe to only include fragments that need processing in this run
    df_to_process = df[df['fragment_id'].isin(fragments_to_process_this_run)].copy()

    # Ensure the frame indices in the filtered dataframe are within the expected range
    # This is a sanity check, the pipeline should handle this filtering.
    # df_to_process = df_to_process[(df_to_process['frame_index'] >= start_frame_index_range) &
    #                               (df_to_process['frame_index'] < start_frame_index_range + n_frames_range)].copy()
    # No, the df should contain original frame indices within the processed range.
    # The filtering should happen upstream in filter_surface_fragments.py or pipeline.py.
    # Let's assume df contains the correct frames based on the range.

    fragments_processed_this_run = 0

    fragments_processed_this_run = 0

    # Processing Loop
    for frag in fragments_to_process_this_run:
        print(f"Processing fragment {frag}...")

        # Process the fragment using the filtered dataframe and the full times array
        frag_events, frag_summary = process_fragment(frag, df_to_process, times, args, n_frames)

        # Append results to processed lists
        processed_events.extend(frag_events)
        if frag_events: # Only append summary if there are events
            processed_summaries.append(frag_summary)

        # Update unprocessed_fragment_ids
        if frag in unprocessed_fragment_ids: # Ensure it's still in the list before removing
             unprocessed_fragment_ids.remove(frag)

        fragments_processed_this_run += 1

        # Save Checkpoint Periodically
        if fragments_processed_this_run % SAVE_INTERVAL == 0:
             save_checkpoint(checkpoint_file, unprocessed_fragment_ids, processed_events, processed_summaries)
             print(f"Checkpoint saved after processing {fragments_processed_this_run} fragments in this run.")


    # Final Output (only if all fragments processed)
    if not unprocessed_fragment_ids: # Check if the list is empty
        print("All fragments processed. Writing final outputs.")
        # Ensure output directory
        os.makedirs(args.output_dir, exist_ok=True)

        # Write CSV of events using processed_events
        df_evt = pd.DataFrame(processed_events)
        csv_out = os.path.join(args.output_dir, 'residence_events.csv')
        df_evt.to_csv(csv_out, index=False)
        print(f"Saved events CSV: {csv_out}")

        # Write JSON summary using processed_summaries
        json_out = os.path.join(args.output_dir, 'residence_events_summary.json')
        with open(json_out, 'w') as fj:
            json.dump(processed_summaries, fj, indent=2)
        print(f"Saved summary JSON: {json_out}")

        # Write Markdown report using processed_summaries
        md_out = os.path.join(args.output_dir, 'residence_events_summary.md')
        with open(md_out, 'w') as fm:
            fm.write("# Residence Event Analysis\n\n")
            fm.write(f"Min residence time: **{args.min_res_time} ns**, ")
            fm.write(f"Max off-gap: **{args.max_off_time} ns**\n\n")
            for s in processed_summaries:
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
            fm.write("\n")
        print(f"Saved Markdown report: {md_out}")

        # Clean up checkpoint file on successful completion
        if os.path.exists(checkpoint_file):
            try:
                os.remove(checkpoint_file)
                print(f"Checkpoint file {checkpoint_file} removed on successful completion.")
            except Exception as e:
                print(f"Error removing checkpoint file {checkpoint_file}: {e}", file=sys.stderr)
    else:
        print(f"Processing interrupted. {len(unprocessed_fragment_ids)} fragments remaining. Checkpoint saved to {checkpoint_file}")


if __name__ == '__main__':
    main()
