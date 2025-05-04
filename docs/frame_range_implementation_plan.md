# Implementation Plan: Adding Frame Range Selection to the Pipeline

This document details the plan to add functionality for specifying a subset of frames from the DCD trajectory for processing in the analysis pipeline.

## Objective

Modify the analysis pipeline to allow users to specify a start and end frame for processing DCD trajectories, enabling analysis of subsets of large simulations.

## Files to Edit/Modify

1.  `docs/scripts/pipeline.md`: Update documentation to describe the new feature.
2.  `config.py`: Add optional `start_frame` and `end_frame` parameters to relevant analysis steps.
3.  `src/simulation_info.py`: Modify to accept `start_frame` and `end_frame` arguments and process only the specified range.
4.  `src/fragment_surface_analysis.py`: Modify to accept `start_frame` and `end_frame` arguments and process only the specified range, while also integrating with existing checkpointing.
5.  `pipeline.py`: Modify to read `start_frame` and `end_frame` from `config.py` and pass them as arguments to `simulation_info.py` and `fragment_surface_analysis.py`.

## Detailed Edits

### 1. `docs/scripts/pipeline.md`

*   Add a new section titled "Specifying Frame Range for Analysis".
*   Explain that users can define `start_frame` and `end_frame` parameters in `config.py` for `simulation_info` and `fragment_surface_analysis` steps.
*   Clarify that `start_frame` and `end_frame` are 0-based indices.
*   Describe how `pipeline.py` will handle these parameters.
*   Mention the impact on checkpointing in `fragment_surface_analysis.py` (it will resume within the specified frame range).

### 2. `config.py`

*   Modify the `config.py` file to include optional `start_frame` and `end_frame` parameters within the configuration dictionaries for the `simulation_info` and `fragment_surface_analysis` steps.
*   These parameters should be optional integer values. If they are not provided or set to a specific indicator (like `None` or -1), the scripts should default to processing the entire trajectory.

### 3. `src/simulation_info.py`

*   Add command-line arguments `--start-frame` and `--end-frame` using `argparse`. These should be optional integers.
*   Modify the DCD reading logic (likely using MDAnalysis) to load the trajectory with the specified `start_frame` and `end_frame`. MDAnalysis `Universe.select_atoms(...).trajectory[start:end:step]` slicing can be used.
*   Ensure that the calculation of `n_frames`, `start_time_ns`, `end_time_ns`, and `total_time_ns` correctly reflects the specified frame range, not the entire DCD file.
*   **New Requirement:** Modify the script to also calculate and save simulation information (total frames, total time, etc.) for the *full* DCD trajectory, in addition to the information for the specified frame range. This might involve reading the full trajectory metadata separately or storing it before applying the frame range slice. The output format should clearly distinguish between full simulation info and frame-range-specific info (e.g., separate sections in JSON/Markdown, or distinct columns/rows in CSV).

### 4. `src/fragment_surface_analysis.py`

*   Add command-line arguments `--start-frame` and `--end-frame` using `argparse`. These should be optional integers.
*   Modify the DCD reading logic to load the trajectory with the specified `start_frame` and `end_frame`.
*   **Crucially**, integrate the existing checkpointing logic with the new frame range functionality.
    *   When loading a checkpoint, the script should determine the last processed frame *within the specified range*.
    *   The script should then resume processing from the frame immediately following the last processed frame *within that range*.
    *   If the loaded checkpoint's last frame is outside the new specified range, the script should likely ignore the checkpoint and start from the beginning of the new range, or handle this edge case appropriately (e.g., raise a warning/error).
    *   The checkpoint file should potentially store the `start_frame` and `end_frame` it was created under to help with resuming within the correct range.
*   Ensure that frame indices and time calculations in the output files (`surface_analysis.csv`, `surface_summary.json`, `surface_summary.md`) are consistent with the processed frame range.

### 5. `pipeline.py`

*   Import the updated `config.py`.
*   When constructing the `subprocess` commands for `simulation_info.py` and `fragment_surface_analysis.py`:
    *   Read the `start_frame` and `end_frame` values from the `config.py` dictionary for the respective steps.
    *   If these values are present (not `None`), include `--start-frame <value>` and `--end-frame <value>` in the command-line arguments passed to the scripts.

## Subtask Delegation

Once this plan is documented, a subtask can be created (e.g., using `new_task` in `sparc` mode) to perform the actual code modifications in `src/simulation_info.py`, `src/fragment_surface_analysis.py`, `config.py`, and `pipeline.py`, and to update `docs/scripts/pipeline.md`.