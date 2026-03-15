# MissileAlerts Analysis

This repository includes `Alert_distribution.jl`, a Julia script that:

1. Loads alerts from `alerts.jld2`
2. Filters by city and optional date range
3. Computes inter-alert intervals
4. Fits a Poisson-process baseline (exponential inter-arrival model)
5. Runs goodness-of-fit and refractory-period analyses
6. Generates summary plots and a refractory scan table

## Requirements

Install packages in Julia:

```julia
using Pkg
Pkg.add(["JLD2", "DataFrames", "CSV", "Plots", "TimeZones"])
```

(`Dates`, `Statistics`, and `Random` are part of Julia standard libraries.)

To generate `alerts.jld2` with `Alerts.jl`, also install:

```julia
using Pkg
Pkg.add(["HTTP", "JSON"])
```

## Generate `alerts.jld2`

Run:

```bash
julia Alerts.jl
```

What this does:

- Downloads cities/areas/polygons metadata from the Tzeva Adom API.
- Downloads historical alerts.
- Writes everything to `alerts.jld2` (including the `alerts` table used by `Alert_distribution.jl`).

Note: `Alerts.jl` currently rebuilds the full `alerts.jld2` file from the remote API and does not take a date-range argument.

Credit: `Alerts.jl` was written by **Yotam Ohad** (GitHub user: **yohad**).

## How to run

### Full command format

```bash
julia Alert_distribution.jl [city_pattern] [input_jld2] [output_png] [start_date] [end_date] [refractory_start_min] [refractory_end_min] [scan_step_min] [scan_window_width_min]
```

### Argument meanings

- `city_pattern` (default: `"תל אביב"`): substring used to match city names
- `input_jld2` (default: `"alerts.jld2"`): input data file
- `output_png` (default: `"tel_aviv_interval_distribution.png"`): main output plot file
- `start_date` (optional): `YYYY-MM-DD`
- `end_date` (optional): `YYYY-MM-DD`
- `refractory_start_min` (default: `30`): refractory window start in minutes
- `refractory_end_min` (default: `60`): refractory window end in minutes
- `scan_step_min` (default: `10`): step size for automatic scan over window starts
- `scan_window_width_min` (default: `30`): fixed width of scan windows in minutes

## What is generated

For `output_png = NAME.png`, the script generates:

1. `NAME.png`
   - Left panel: histogram of `log10` inter-alert intervals
   - Right panel: empirical CDF vs fitted exponential CDF + KS summary

2. `NAME_refractory_scan.csv`
   - Table with one row per scanned window
   - Includes: window bounds, counts, empirical conditional probability, Poisson expectation, refractory ratio, lower-tail p-value

3. `NAME_refractory_scan.png`
   - Top panel: refractory ratio (`empirical / Poisson`) vs window start
   - Bottom panel: one-sided lower-tail p-value vs window start

## Example runs

```bash
julia Alert_distribution.jl "תל אביב" "alerts.jld2" "Tel-Aviv-later.png" 2026-03-04 2026-03-15 30 60 10 60
julia Alert_distribution.jl "תל אביב" "alerts.jld2" "Tel-Aviv-earlier.png" 2026-02-28 2026-03-03 30 60 10 60
```

These commands compare two periods (an earlier and a later phase) for Tel-Aviv alerts, and generate both the main distribution plot and automatic refractory scan outputs for each period.

## Interpreting key outputs

- **KS p-value**: tests whether inter-arrival times are consistent with the fitted exponential model.
- **Refractory ratio < 1**: fewer events than Poisson expectation in that window.
- **Lower-tail p-value** (scan and selected window): evidence strength for "lower-than-Poisson" in that window.

Use these together (effect size + p-value) when comparing periods.
