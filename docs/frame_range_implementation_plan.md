# Implementation Plan: Adding Frame Range Selection to the Pipeline

This document details the plan to add functionality for specifying a subset of frames from the DCD trajectory for processing in the analysis pipeline.

## Objective

Modify the analysis pipeline to allow users to specify a start and end frame for processing DCD trajectories, enabling analysis of subsets of large simulations.

## Files to Edit/Modify

1.  `docs/scripts/pipeline.md`: Update documentation to describe the new feature.
2.  `config.py`: Add optional `start_frame` and `end_frame` parameters to relevant analysis steps.
3.  `src/simulation_info.py`: Modify to accept optional `start_frame` and `end_frame` arguments and process only the specified range, while also reporting full simulation info.
4.  `src/fragment_surface_analysis.py`: Modify to accept optional `start_frame` and `end_frame` arguments and process only the specified range, while also integrating with existing checkpointing.
5.  `src/fragment_surface_analysis_test_restart.py`: Modify to accept optional `start_frame` and `end_frame` arguments to allow testing of checkpointing within a specific range.
6.  `pipeline.py`: Modify to read optional `start_frame` and `end_frame` from `config.py` and pass them as arguments to `simulation_info.py` and `fragment_surface_analysis.py`. Also, pass `start_frame` to `classify_pt_atoms.py` as `--frame-index` when a range is specified.

## Detailed Edits

### 1. `docs/scripts/pipeline.md`

*   Add a new section titled "Specifying Frame Range and Timestep".
*   Explain that users can define `start_frame`, `end_frame`, and `timestep_fs` parameters in `config.py` for the `simulation_info` step.
*   Clarify that `start_frame` and `end_frame` are 0-based indices.
*   Explain that `timestep_fs` allows overriding the timestep read from the DCD header.
*   Describe how `pipeline.py` will handle these parameters and pass them to `simulation_info.py`.
*   Mention the impact on `fragment_surface_analysis.py` (it will process the specified frame range) and checkpointing (it will resume within the specified frame range).

### 2. `config.py`

*   Modify the `config.py` file to include optional `start_frame`, `end_frame`, and `timestep_fs` parameters within the configuration dictionary for the `simulation_info` step.
*   Also include optional `start_frame` and `end_frame` parameters for the `fragment_surface_analysis` step.
*   These parameters should be optional integer or float values. If they are not provided or set to a specific indicator (like `None` or -1), the scripts should default to processing the entire trajectory and reading the timestep from the DCD header.

### 3. `src/simulation_info.py`

*   Add command-line arguments `--start-frame`, `--end-frame`, and `--timestep-fs` using `argparse`. These should be optional integers (`start-frame`, `end-frame`) and an optional float (`timestep-fs`).
*   Modify the script to check if `--timestep-fs` is provided. If yes, use this value for `dt_fs`. If not, read `dt_fs` from the DCD header (`u.trajectory.dt`).
*   Modify the DCD reading logic (likely using MDAnalysis) to load the trajectory with the specified `start_frame` and `end_frame`. MDAnalysis `Universe.select_atoms(...).trajectory[start:end:step]` slicing can be used.
*   Ensure that the calculation of `n_frames`, `start_time_ns`, `end_time_ns`, and `total_time_ns` correctly reflects the specified frame range, using the determined `dt_fs` and `steps_per_frame`.
*   **New Requirement:** Modify the script to also calculate and save simulation information (total frames, total time, etc.) for the *full* DCD trajectory, in addition to the information for the specified frame range.
*   **Output Structure for Downstream Compatibility:** The output files (`simulation_info.csv`, `simulation_info.json`, `simulation_info.md`) must be structured to provide both the full simulation details and the details for the processed frame range in a way that is easily accessible and clearly distinguishable.
    *   **`simulation_info.json`**: This file should contain a top-level structure (e.g., a dictionary) with two main keys: one for "full_simulation" and one for "processed_range". Each key should contain a dictionary with the relevant metrics (`n_frames`, `start_time_ns`, `end_time_ns`, `total_time_ns`, `integration_timestep_fs`, `timestep_source`, etc.) for that scope. Downstream scripts like `residence_event_analysis.py` will need to be updated to read the "processed_range" section from this JSON.
    *   **`simulation_info.csv`**: This file could include columns that clearly indicate whether a row pertains to the "full_simulation" or the "processed_range". Alternatively, it could have separate sections or be structured to prioritize the "processed_range" info for easy reading by `binding_metrics.py`, while still including full simulation info for reference. A simple approach is to have columns like `scope` ('full' or 'range'), `n_frames`, `total_time_ns`, `integration_timestep_fs`, `timestep_source`, etc. `binding_metrics.py` will need to be updated to read the `total_time_ns` value specifically for the 'processed_range' scope.
    *   **`simulation_info.md`**: The Markdown report should clearly label sections for "Full Simulation Information" and "Processed Frame Range Information", presenting the relevant metrics under each heading, including the timestep and its source.

### 4. `src/fragment_surface_analysis.py`

*   Add command-line arguments `--start-frame` and `--end-frame` using `argparse`. These should be optional integers.
*   Modify the DCD reading logic to load the trajectory with the specified `start_frame` and `end_frame`.
*   **Crucially**, integrate the existing checkpointing logic with the new frame range functionality.
    *   When loading a checkpoint, the script should determine the last processed frame *within the specified range*.
    *   The script should then resume processing from the frame immediately following the last processed frame *within that range*.
    *   If the loaded checkpoint's last frame is outside the new specified range, the script should likely ignore the checkpoint and start from the beginning of the new range, or handle this edge case appropriately (e.g., raise a warning/error).
    *   The checkpoint file should potentially store the `start_frame` and `end_frame` it was created under to help with resuming within the correct range.
*   Ensure that frame indices and time calculations in the output files (`surface_analysis.csv`, `surface_summary.json`, `surface_summary.md`) are consistent with the processed frame range.

### 5. `src/fragment_surface_analysis_test_restart.py`

*   Add command-line arguments `--start-frame` and `--end-frame` using `argparse`. These should be optional integers.
*   Modify the DCD reading logic to load the trajectory with the specified `start_frame` and `end_frame` for testing purposes.

### 6. `pipeline.py`

*   Import the updated `config.py`.
*   When constructing the `subprocess` commands for `simulation_info.py` and `fragment_surface_analysis.py`:
    *   Read the `start_frame`, `end_frame`, and `timestep_fs` values from the `config.py` dictionary for the respective steps.
    *   If `start_frame`, `end_frame`, or `timestep_fs` values are present (not `None`), include the corresponding `--start-frame <value>`, `--end-frame <value>`, and `--timestep-fs <value>` in the command-line arguments passed to `simulation_info.py`.
    *   If `start_frame` or `end_frame` values are present (not `None`), include the corresponding `--start-frame <value>` and `--end-frame <value>` in the command-line arguments passed to `fragment_surface_analysis.py`.
    *   If `start_frame` is present (not `None`), pass it to `classify_pt_atoms.py` as `--frame-index <start_frame>`.

## Impact on Downstream Scripts (`residence_event_analysis.py`, `binding_metrics.py`)

`residence_event_analysis.py` and `binding_metrics.py` do not directly read the DCD trajectory. They rely on output files generated by preceding steps, specifically:

*   `residence_event_analysis.py` uses `simulation_info.json` (for time information) and `surface_filtered.csv` (for filtered frame data).
*   `binding_metrics.py` uses `simulation_info.csv` (for total simulation time) and outputs from `residence_event_analysis.py` and `classify_pt_atoms.py`.

The plan is to modify `src/simulation_info.py` to calculate and output time metrics (like `n_frames`, `start_time_ns`, `end_time_ns`, `total_time_ns`, `integration_timestep_fs`, `timestep_source`) that are specific to the *selected frame range*, in addition to the full simulation info.

By providing range-specific time information in the `simulation_info` outputs, the downstream scripts (`residence_event_analysis.py` and `binding_metrics.py`) will process and report results that are consistent with the specified frame range.

**Important Note:** While the core logic of these downstream scripts for calculating events and metrics does not need to change, they *will* require minor modifications to correctly parse the time information from the new, structured output of `simulation_info.py`. Specifically:

*   `residence_event_analysis.py` will need to be updated to read the time metrics (`n_frames`, `start_time_ns`, `time_between_frames_ns`) from the `"processed_range"` key in the `simulation_info.json` file.
*   `binding_metrics.py` will need to be updated to read the `total_time_ns` value specifically for the `"processed_range"` scope from the `simulation_info.csv` file (based on the chosen CSV structure, e.g., filtering by a 'scope' column).

These are minor parsing updates, not changes to the core analysis algorithms, but they are necessary to ensure the scripts use the correct time context for the processed frame range.

## Subtask Delegation

Once this plan is documented, a subtask can be created (e.g., using `new_task` in `sparc` mode) to perform the actual code modifications in `src/simulation_info.py`, `src/fragment_surface_analysis.py`, `config.py`, and `pipeline.py`, and to update `docs/scripts/pipeline.md`.