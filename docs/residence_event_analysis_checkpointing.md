# `residence_event_analysis.py` Checkpointing Implementation Plan

This document outlines the implementation plan for adding checkpointing capabilities to the `residence_event_analysis.py` script. This will allow the script to save its progress periodically and resume from the last saved state if interrupted.

## Objective

Enable `residence_event_analysis.py` to save and load its processing state to allow for resuming analysis from an interruption point.

## Requirements

1.  The script must save its current state (unprocessed fragment IDs, accumulated events, accumulated summaries) to a checkpoint file.
2.  On startup, the script must check for and load an existing checkpoint file if found.
3.  If a checkpoint is loaded, the script must resume processing from the point indicated by the checkpointed state.
4.  The checkpointing mechanism should be robust to file errors during saving/loading.
5.  The checkpoint file path should be configurable or follow a predictable pattern within the output directory.

## State to Checkpoint

The following data needs to be saved to the checkpoint file:

-   `unprocessed_fragment_ids`: A list of integer fragment IDs that have not yet been fully processed.
-   `processed_events`: A list of dictionaries, where each dictionary represents a residence event found for fragments processed so far. This corresponds to the `all_events` list in the current script.
-   `processed_summaries`: A list of dictionaries, where each dictionary represents the summary statistics for fragments processed so far. This corresponds to the `summaries` list in the current script.

## Checkpoint File

-   The checkpoint file will be a JSON file.
-   A suggested default path is `[output_dir]/residence_events_checkpoint.json`.
-   The script should handle cases where the file does not exist or is corrupted.

## Pseudocode

```pseudocode
FUNCTION main():
    # Argument Parsing
    args = parse_arguments() # Add argument for checkpoint file path

    # Define checkpoint file path
    checkpoint_file = os.path.join(args.output_dir, "residence_events_checkpoint.json")

    # Initialize state variables
    unprocessed_fragment_ids = []
    processed_events = []
    processed_summaries = []
    start_fragment_index = 0 # Index in the original fragment_ids list

    # Load Checkpoint
    IF checkpoint_file exists AND is not empty:
        TRY:
            checkpoint_data = load_json(checkpoint_file)
            unprocessed_fragment_ids = checkpoint_data["unprocessed_fragment_ids"]
            processed_events = checkpoint_data["processed_events"]
            processed_summaries = checkpoint_data["processed_summaries"]
            PRINT "Resuming from checkpoint."
        EXCEPT Exception as e:
            PRINT f"Error loading checkpoint file: {e}. Starting from scratch."
            # Fallback to starting from scratch
            df = pd.read_csv(args.filtered_csv)
            original_fragment_ids = sorted(df['fragment_id'].unique())
            unprocessed_fragment_ids = [int(x) for x in original_fragment_ids]
            processed_events = []
            processed_summaries = []
    ELSE:
        PRINT "No checkpoint found. Starting from scratch."
        df = pd.read_csv(args.filtered_csv)
        original_fragment_ids = sorted(df['fragment_id'].unique())
        unprocessed_fragment_ids = [int(x) for x in original_fragment_ids]
        processed_events = []
        processed_summaries = []

    # Load simulation info (needed regardless of checkpoint)
    sim = load_json(args.sim_info_json)
    n_frames = int(sim.get("n_frames"))
    start_ns = float(sim.get("start_time_ns"))
    dt_ns = float(sim.get("time_between_frames_ns"))
    times = start_ns + np.arange(n_frames) * dt_ns

    # Load filtered data (needed regardless of checkpoint)
    df = pd.read_csv(args.filtered_csv)

    # Determine the starting point in the fragment processing loop
    # Find the index of the first unprocessed fragment in the original list
    # This is needed if we want to iterate through the original list and skip processed ones
    # Alternatively, iterate directly through unprocessed_fragment_ids

    # Let's iterate directly through unprocessed_fragment_ids for simplicity
    fragments_to_process = list(unprocessed_fragment_ids) # Create a copy to modify unprocessed_fragment_ids during iteration

    # Processing Loop
    FOR frag IN fragments_to_process:
        PRINT f"Processing fragment {frag}..."

        # --- Existing processing logic for a single fragment ---
        # (Lines 81-139 from original script, adapted to use 'frag')
        on_frames = set(df.loc[df['fragment_id']==frag, 'frame_index'])
        on_series = np.array([i in on_frames for i in range(n_frames)], dtype=bool)

        # Detect raw runs
        runs = []  # list of (start_idx, end_idx)
        in_run = False
        FOR i, on IN enumerate(on_series):
            IF on AND NOT in_run:
                start = i
                in_run = TRUE
            ELIF NOT on AND in_run:
                end = i-1
                runs.append((start, end))
                in_run = FALSE
        IF in_run:
            runs.append((start, n_frames-1))

        # Compute run durations and filter by min_res_time
        filtered_runs = []
        FOR (s, e) IN runs:
            t0 = times[s]
            t1 = times[e]
            dur = t1 - t0
            IF dur >= args.min_res_time:
                filtered_runs.append({'start': s, 'end': e, 'start_ns': t0, 'end_ns': t1, 'dur_ns': dur})

        # Merge runs separated by small off-gaps
        merged = []
        FOR run IN filtered_runs:
            IF NOT merged:
                merged.append(run.copy())
            ELSE:
                prev = merged[-1]
                gap = run['start_ns'] - prev['end_ns']
                IF gap <= args.max_off_time:
                    prev['end'] = run['end']
                    prev['end_ns'] = run['end_ns']
                    prev['dur_ns'] = prev['end_ns'] - prev['start_ns']
                ELSE:
                    merged.append(run.copy())

        # Record events for this fragment and append to processed_events
        frag_events = []
        FOR eid, run IN enumerate(merged):
            evt = {
                'fragment_id': frag,
                'event_id': eid,
                'start_frame': int(run['start']),
                'end_frame': int(run['end']),
                'start_time_ns': float(run['start_ns']),
                'end_time_ns': float(run['end_ns']),
                'duration_ns': float(run['dur_ns'])
            }
            processed_events.append(evt)
            frag_events.append(evt)

        # Summarize for this fragment and append to processed_summaries
        durations = [e['duration_ns'] for e in frag_events]
        IF durations:
            processed_summaries.append({
                'fragment_id': frag,
                'n_events': len(durations),
                'total_residence_ns': float(sum(durations)),
                'mean_event_ns': float(np.mean(durations)),
                'events': [
                    {'event_id': e['event_id'],
                     'start_time_ns': e['start_time_ns'],
                     'end_time_ns': e['end_time_ns'],
                     'duration_ns': e['duration_ns']}
                    FOR e IN frag_events
                ]
            })
        ELSE:
            processed_summaries.append({
                'fragment_id': frag,
                'n_events': 0,
                'total_residence_ns': 0.0,
                'mean_event_ns': 0.0,
                'events': []
            })
        # --- End existing processing logic ---

        # Update unprocessed_fragment_ids
        unprocessed_fragment_ids.remove(frag)

        # Save Checkpoint Periodically (e.g., after every N fragments)
        # Define N (e.g., 10 or 100)
        # IF (fragments_processed_count % N == 0): # Need a counter for processed fragments in this run
        #     save_checkpoint(checkpoint_file, unprocessed_fragment_ids, processed_events, processed_summaries)
        # For simplicity in pseudocode, let's save after each fragment for now.
        # A more robust implementation would save periodically or on signal.
        save_checkpoint(checkpoint_file, unprocessed_fragment_ids, processed_events, processed_summaries)
        PRINT f"Checkpoint saved after processing fragment {frag}."


    # Final Output (only if all fragments processed)
    IF NOT unprocessed_fragment_ids: # Check if the list is empty
        PRINT "All fragments processed. Writing final outputs."
        # Ensure output directory
        os.makedirs(args.output_dir, exist_ok=True)

        # Write CSV of events using processed_events
        df_evt = pd.DataFrame(processed_events)
        csv_out = os.path.join(args.output_dir, 'residence_events.csv')
        df_evt.to_csv(csv_out, index=False)
        PRINT f"Saved events CSV: {csv_out}"

        # Write JSON summary using processed_summaries
        json_out = os.path.join(args.output_dir, 'residence_events_summary.json')
        WITH open(json_out, 'w') as fj:
            json.dump(processed_summaries, fj, indent=2)
        PRINT f"Saved summary JSON: {json_out}"

        # Write Markdown report using processed_summaries
        md_out = os.path.join(args.output_dir, 'residence_events_summary.md')
        WITH open(md_out, 'w') as fm:
            fm.write("# Residence Event Analysis\n\n")
            fm.write(f"Min residence time: **{args.min_res_time} ns**, ")
            fm.write(f"Max off-gap: **{args.max_off_time} ns**\n\n")
            FOR s IN processed_summaries:
                fm.write(f"## Fragment {s['fragment_id']}\n")
                fm.write(f"- Number of events: **{s['n_events']}**\n")
                fm.write(f"- Total residence: **{s['total_residence_ns']:.3f} ns**\n")
                fm.write(f"- Mean event duration: **{s['mean_event_ns']:.3f} ns**\n")
                IF s['events']:
                    fm.write("### Events:\n")
                    FOR e IN s['events']:
                        fm.write(
                            f"- Event {e['event_id']}: {e['duration_ns']:.3f} ns ("
                            f"{e['start_time_ns']:.3f} → {e['end_time_ns']:.3f} ns)\n"
                        )
                fm.write("\n")
        PRINT f"Saved Markdown report: {md_out}")

        # Optional: Clean up checkpoint file on successful completion
        IF checkpoint_file exists:
            os.remove(checkpoint_file)
            PRINT "Checkpoint file removed on successful completion."

FUNCTION save_checkpoint(file_path, unprocessed_fragment_ids, processed_events, processed_summaries):
    checkpoint_data = {
        "unprocessed_fragment_ids": unprocessed_fragment_ids,
        "processed_events": processed_events,
        "processed_summaries": processed_summaries
    }
    TRY:
        WITH open(file_path, 'w') as f:
            json.dump(checkpoint_data, f, indent=2)
    EXCEPT Exception as e:
        PRINT f"Error saving checkpoint file: {e}"
        # Consider more robust error handling, e.g., saving to a temporary file first

FUNCTION load_json(file_path):
    WITH open(file_path, 'r') as f:
        RETURN json.load(f)

```

## TDD Anchors

1.  **Test Case: No Checkpoint File**
    *   **Input:** Run the script with no existing checkpoint file.
    *   **Expected Output:** Script starts processing from the first fragment, `unprocessed_fragment_ids` initially contains all fragment IDs, `processed_events` and `processed_summaries` are empty. A checkpoint file is created during execution.
2.  **Test Case: Valid Checkpoint File**
    *   **Input:** Run the script with a valid checkpoint file containing a partial state (some fragments processed).
    *   **Expected Output:** Script loads the state, `unprocessed_fragment_ids` contains the remaining fragment IDs, `processed_events` and `processed_summaries` contain the data from the checkpoint. Processing resumes from the first unprocessed fragment. The final output files are correct and complete.
3.  **Test Case: Corrupted Checkpoint File**
    *   **Input:** Run the script with a checkpoint file that contains invalid JSON.
    *   **Expected Output:** Script catches the loading error, prints a warning, and starts processing from scratch as if no checkpoint was found.
4.  **Test Case: Interruption and Resume**
    *   **Input:** Run the script, interrupt it after some fragments have been processed (before completion). Then run the script again.
    *   **Expected Output:** The second run loads the checkpoint and resumes processing from where it left off. The final output files are identical to a successful run without interruption.
5.  **Test Case: Checkpoint File Content**
    *   **Input:** Run the script and let it process a few fragments.
    *   **Expected Output:** Verify the content of the checkpoint file matches the expected state (`unprocessed_fragment_ids`, `processed_events`, `processed_summaries`) after those fragments have been processed.

## Integration with `pipeline.py`

The `pipeline.py` script will need to be updated to pass the appropriate output directory to `residence_event_analysis.py` so the script can determine the default checkpoint file path. The pipeline-level `.done` file mechanism will still be used to track whether the `residence_event_analysis` step has completed successfully (i.e., all fragments processed and final outputs written). The internal checkpointing allows the script to be restarted by the pipeline and resume its work within that step.