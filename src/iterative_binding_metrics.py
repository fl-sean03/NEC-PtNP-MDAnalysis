import subprocess
import pandas as pd
import json
import argparse
import os
import math # Import math module

def run_iterative_binding_metrics(events_csv, classification_csv, properties_csv, sim_info_json, min_res_time, output_csv, verbose=False):
    """
    Reruns binding_metrics.py iteratively over increasing frame ranges
    and collects facet-specific binding metrics.
    """
    if verbose:
        print("Starting iterative binding metrics analysis...")
    # Read simulation info to get total number of frames
    with open(sim_info_json, 'r') as f:
        sim_info = json.load(f)
    # Access n_frames from the nested 'full_simulation' key
    total_frames = sim_info['full_simulation']['n_frames']
    time_between_frames_ns = sim_info['full_simulation']['time_between_frames_ns'] # Read time between frames

    # Calculate the minimum number of frames required to meet min_res_time
    min_frames_for_res_time = math.ceil(min_res_time / time_between_frames_ns)
    # Ensure the starting frame is at least 1
    start_frame_range = max(1, min_frames_for_res_time)

    # Initialize output CSV
    output_df = pd.DataFrame(columns=['end_frame', 'facet', 'tau_facet_ns', 'K_D_facet', 'DeltaG_kJ_per_mol'])
    
    # Ensure the output directory exists
    output_dir = os.path.dirname(output_csv)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_df.to_csv(output_csv, index=False)

    temp_output_csv_name = "temp_facet_results.csv"
    temp_output_csv = os.path.join(output_dir, temp_output_csv_name)

    # Iterate over increasing end frames with a specified step, starting from min_frames_for_res_time
    frame_step = args.frame_step # Get frame_step from parsed arguments
    for end_frame in range(start_frame_range, total_frames + 1, frame_step):
        # Ensure the last frame is included even if total_frames is not a multiple of frame_step
        current_end_frame = min(end_frame, total_frames)

        # Check if this is the last iteration covering the full range
        is_full_range_iteration = (current_end_frame == total_frames)

        if verbose:
            if is_full_range_iteration:
                print(f"Running binding metrics for the full simulation range (frames [0, {total_frames}))")
            else:
                print(f"Running binding metrics for frames [0, {current_end_frame})")

        # Construct the command to run binding_metrics.py
        command = [
            "python", "src/binding_metrics.py",
            "--residence_events_csv", events_csv,
            "--filtered_surface_csv", os.path.join(os.path.dirname(os.path.dirname(classification_csv)), "surface_filtered", "surface_filtered.csv"), # Dynamically construct path
            "--simulation_properties_csv", properties_csv,
            "--sim-info-json", sim_info_json, # Add sim-info-json argument
            "--output_facet_csv", temp_output_csv # Always output to temporary file
        ]

        # Add timeframe arguments only if it's not the full range iteration
        if not is_full_range_iteration:
            command.extend([
                "--analysis-start-frame", "0",
                "--analysis-end-frame", str(current_end_frame),
                "--min-res-time", str(min_res_time)
            ])

        if verbose:
            print(f"Executing command: {' '.join(command)}")

        # Execute the command
        try:
            subprocess.run(command, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            if verbose:
                print(f"Error running binding_metrics.py for frame range [0, {end_frame}):")
                print(e.stdout)
                print(e.stderr)
            continue # Continue to the next frame range even if one fails

        # Read the temporary output and append to the main output CSV
        if os.path.exists(temp_output_csv):
            try:
                temp_df = pd.read_csv(temp_output_csv)
                temp_df['end_frame'] = current_end_frame # Use current_end_frame
                temp_df = temp_df[['end_frame', 'facet', 'tau_facet_ns', 'K_D_facet', 'DeltaG_kJ_per_mol']]
                temp_df.to_csv(output_csv, mode='a', header=False, index=False)
            except pd.errors.EmptyDataError:
                if verbose:
                    print(f"Warning: Temporary output file {temp_output_csv} is empty for frame range [0, {current_end_frame}). Skipping.")
            except Exception as e:
                if verbose:
                    print(f"Error processing temporary output file {temp_output_csv} for frame range [0, {current_end_frame}): {e}")
            finally:
                # Always attempt to remove the temporary file
                if os.path.exists(temp_output_csv):
                    os.remove(temp_output_csv)
        else:
            if verbose:
                print(f"Warning: {temp_output_csv} not found after running for frame range [0, {current_end_frame})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run binding_metrics.py iteratively over increasing frame ranges.")
    parser.add_argument("--events-csv", required=True, help="Path to residence_events.csv")
    parser.add_argument("--classification-csv", required=True, help="Path to coordination_numbers.csv")
    parser.add_argument("--properties-csv", required=True, help="Path to simulation_info.csv")
    parser.add_argument("--sim-info-json", required=True, help="Path to simulation_info.json")
    parser.add_argument("--min-res-time", type=float, required=True, help="Minimum residence time (in ns) used in residence_event_analysis.")
    parser.add_argument("--output-csv", required=True, help="Path for the output CSV file to save iterative results.")
    parser.add_argument("--frame-step", type=int, default=1, help="Number of frames to increment the analysis window by in each iteration.")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output.")

    args = parser.parse_args()

    run_iterative_binding_metrics(
        args.events_csv,
        args.classification_csv,
        args.properties_csv,
        args.sim_info_json,
        args.min_res_time,
        args.output_csv,
        args.verbose # Pass the verbose flag
    )