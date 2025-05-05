# Specification: Timeframe-Specific Analysis in binding_metrics.py

**Version:** 1.0
**Date:** 2025-05-05

## 1. Introduction & Motivation

Currently, the `binding_metrics.py` script calculates binding metrics based on the entire dataset provided in the input `residence_events.csv` file, using the total simulation time derived from `simulation_info.csv`.

There is a need to analyze binding behavior within specific time windows of the simulation without re-running the entire preceding pipeline (especially `fragment_surface_analysis.py` and `residence_event_analysis.py`) for just that window.

This update will add functionality to `binding_metrics.py` to accept optional start and end time arguments. When provided, the script will filter the input residence events to include only those occurring within the specified timeframe, adjust their durations accordingly, re-apply the minimum residence time threshold, and calculate metrics based only on this time window.

## 2. Overview of Changes

*   **`src/binding_metrics.py`:**
    *   Add new command-line arguments: `--analysis-start-time-ns`, `--analysis-end-time-ns`, `--min-res-time`.
    *   Implement logic to filter residence events based on the specified time window.
    *   Implement logic to adjust the duration of events that overlap the time window boundaries.
    *   Implement logic to filter out events whose *adjusted* duration falls below the provided `--min-res-time`.
    *   Calculate and use the duration of the specified time window as the effective simulation time (`T_sim`) for metric calculations (`T_off_ns`, `K_D`).
    *   Ensure the script retains its original behavior if the new time arguments are not provided.
*   **`pipeline.py`:**
    *   Modify the execution step for `binding_metrics.py` to retrieve the `min_res_time` value (used in the `residence_event_analysis` step) from `config.py` and pass it via the new `--min-res-time` argument.
*   **`docs/scripts/binding_metrics.md`:**
    *   Update documentation to reflect the new arguments and functionality.

## 3. Detailed Implementation Steps for `src/binding_metrics.py`

These steps should be implemented iteratively.

**Step 3.1: Add Command-Line Arguments**

*   Modify the `argparse` setup within the `main()` function.
*   Add three new optional arguments:
    *   `--analysis-start-time-ns`: Type `float`, default `None`. Help text: "Optional start time (in ns) for the analysis window."
    *   `--analysis-end-time-ns`: Type `float`, default `None`. Help text: "Optional end time (in ns) for the analysis window."
    *   `--min-res-time`: Type `float`, default `None`. Help text: "Minimum residence time (in ns) threshold used in residence_event_analysis. Required if timeframe filtering is active."
*   Ensure these arguments are parsed and stored (e.g., in `args.analysis_start_time_ns`, `args.analysis_end_time_ns`, `args.min_res_time`).

**Step 3.2: Add Conditional Logic for Timeframe Analysis**

*   Implement a check early in `main()` to see if both `args.analysis_start_time_ns` and `args.analysis_end_time_ns` are provided and are valid numbers.
    *   If they are provided, also check if `args.min_res_time` is provided. If not, raise an error or print a warning and exit, as it's required for correct filtering.
    *   If they are *not* provided, the script should proceed with its original logic using the full `events` DataFrame and `T_sim` from `simulation_info.csv`.
*   All subsequent filtering and adjustment steps (3.3 - 3.6) should be placed within the `if` block that executes when timeframe arguments are provided.

**Step 3.3: Implement Timeframe Filtering**

*   **Input:** The `events` DataFrame loaded from `residence_events.csv`.
*   **Condition:** Keep events where the event overlaps with the analysis window `[start_time, end_time]`. The condition is:
    `(events['start_time_ns'] < args.analysis_end_time_ns) & (events['end_time_ns'] > args.analysis_start_time_ns)`
*   **Action:** Create a new DataFrame (e.g., `filtered_events`) containing only the rows that satisfy the condition.
*   **Logging:** Add a verbose print statement indicating that timeframe filtering is active and reporting the number of events before and after filtering.

**Step 3.4: Implement Duration Adjustment**

*   **Input:** The `filtered_events` DataFrame from Step 3.3.
*   **Logic:** Calculate the adjusted start and end times clipped to the analysis window, then find the difference.
    *   `clipped_start_ns = events['start_time_ns'].clip(lower=args.analysis_start_time_ns)`
    *   `clipped_end_ns = events['end_time_ns'].clip(upper=args.analysis_end_time_ns)`
    *   `adjusted_duration_ns = clipped_end_ns - clipped_start_ns`
*   **Action:** Add a new column `adjusted_duration_ns` to the `filtered_events` DataFrame. Replace the original `duration_ns` column with these adjusted values for subsequent calculations within this timeframe analysis block.
*   **Note:** Ensure calculations handle potential floating-point inaccuracies if necessary.

**Step 3.5: Implement `min_res_time` Filtering**

*   **Input:** The `filtered_events` DataFrame with the `adjusted_duration_ns` column (or the updated `duration_ns` column).
*   **Condition:** Keep events where the adjusted duration meets the threshold:
    `filtered_events['duration_ns'] >= args.min_res_time` (assuming `duration_ns` was overwritten in Step 3.4, otherwise use `adjusted_duration_ns`).
*   **Action:** Further filter the `filtered_events` DataFrame, keeping only rows that satisfy this condition.
*   **Logging:** Add a verbose print statement reporting the number of events removed because their adjusted duration was below `min_res_time`.

**Step 3.6: Implement Effective Simulation Time Calculation**

*   **Input:** `args.analysis_start_time_ns`, `args.analysis_end_time_ns`.
*   **Logic:** Calculate the duration of the analysis window.
    *   `effective_T_sim = args.analysis_end_time_ns - args.analysis_start_time_ns`
*   **Action:** Store this value in a variable (e.g., `T_sim_eff`). This variable will be used instead of the `T_sim` read from `simulation_info.csv` when calculating `T_off_ns` and `K_D` *within the timeframe analysis block*.
*   **Logging:** Add a verbose print statement indicating the effective simulation time being used for calculations.

**Step 3.7: Integrate with Existing Calculations**

*   Ensure that all subsequent calculations (grouping, aggregation, `T_on_ns`, `N_events`, `tau_mean_ns`, `T_off_ns`, `K_D`, `DeltaG`) use the final `filtered_events` DataFrame (with adjusted durations) and the `T_sim_eff` calculated in Step 3.6.
*   The logic for determining the `dominant_facet` from `filtered_surface_data` likely does not need modification, as facet classification is usually considered static for a fragment. However, ensure the merge operation (`df = events.merge(...)`) uses the filtered `events` DataFrame.

**Step 3.8: Add Logging/Verbose Output**

*   Throughout the new logic (Steps 3.2-3.6), add informative print statements (ideally controlled by the existing `--verbose` flag) to show:
    *   That timeframe analysis is active.
    *   The specified start/end times and min_res_time.
    *   Number of events before/after timeframe filtering.
    *   Number of events before/after min_res_time filtering.
    *   The calculated effective simulation time.

## 4. Detailed Implementation Steps for `pipeline.py`

**Step 4.1: Read `min_res_time` from Config**

*   Locate the section in `pipeline.py` where parameters for the `residence_event_analysis` step are read from the configuration object (imported from `config.py`).
*   Extract the `min_res_time` value used for that step. Store it in a variable.

**Step 4.2: Pass `min_res_time` Argument**

*   Locate the `subprocess` call that executes `src/binding_metrics.py`.
*   Modify the command list (`cmd`) being constructed for the subprocess.
*   Add the `--min-res-time` argument followed by the `min_res_time` value retrieved in Step 4.1.

**Step 4.3: Consider Future Timeframe Arguments (Optional)**

*   Currently, `--analysis-start-time-ns` and `--analysis-end-time-ns` are intended as direct arguments to `binding_metrics.py`.
*   If pipeline-level control of this timeframe is desired later, `pipeline.py` would need to be further modified to:
    *   Read these optional start/end times from a new section in `config.py`.
    *   Pass them as arguments in the `subprocess` call if they are present in the config.
*   *This is not part of the immediate requirement but should be kept in mind for future extensibility.*

## 5. Testing Considerations

Create test cases covering various scenarios:

*   **Baseline:** Run without timeframe arguments; verify output matches original behavior.
*   **Filtering:**
    *   Provide a timeframe that excludes all events.
    *   Provide a timeframe that includes all events.
    *   Provide a timeframe that includes a subset of events.
*   **Boundary Conditions:**
    *   Event starts before, ends within the timeframe.
    *   Event starts within, ends after the timeframe.
    *   Event starts before, ends after (fully encompasses) the timeframe.
    *   Event start/end matches timeframe start/end exactly.
*   **`min_res_time` Interaction:**
    *   An event's original duration is > `min_res_time`, but its truncated duration is < `min_res_time` (should be excluded).
    *   An event's original duration is > `min_res_time`, and its truncated duration is also > `min_res_time` (should be included with adjusted duration).
*   **Input Validation:**
    *   Provide timeframe arguments but omit `--min-res-time`.
    *   Provide invalid start/end times (e.g., start > end).

## 6. Documentation Updates

*   Update `docs/scripts/binding_metrics.md`:
    *   Add the new arguments (`--analysis-start-time-ns`, `--analysis-end-time-ns`, `--min-res-time`) to the "Arguments" section, explaining their purpose, type, and default values.
    *   Add a new subsection explaining the timeframe analysis feature, detailing how events are filtered, durations adjusted, the `min_res_time` threshold is re-applied, and how the effective simulation time is determined for calculations when these arguments are used.
    *   Clarify that `--min-res-time` is required when providing timeframe arguments.