using JLD2
using DataFrames
using CSV
using Dates
using TimeZones
using Statistics
using Random
using Plots
using Plots.PlotMeasures

"""
Plot alert interval distribution for a city pattern from alerts.jld2.

Usage:
	julia Alert_distribution.jl [city_pattern] [input_jld2] [output_png] [start_date] [end_date] [refractory_start_min] [refractory_end_min] [scan_step_min] [scan_window_width_min]

Defaults:
	city_pattern = "תל אביב"
	input_jld2   = "alerts.jld2"
	output_png   = "tel_aviv_interval_distribution.png"
	start_date   = (none)
	end_date     = (none)
	refractory_start_min = 30
	refractory_end_min   = 60
	scan_step_min        = 10
	scan_window_width_min = 30

Date format:
	YYYY-MM-DD (e.g. 2026-02-28 2026-03-03)
"""

function ks_statistic_exponential(samples::Vector{Float64}, λ::Float64)
	x = sort(samples)
	n = length(x)
	if n == 0
		error("KS test requires at least 1 sample")
	end

	Dplus = 0.0
	Dminus = 0.0
	for i in 1:n
		F = 1 - exp(-λ * x[i])
		Dplus = max(Dplus, i / n - F)
		Dminus = max(Dminus, F - (i - 1) / n)
	end
	max(Dplus, Dminus)
end

function ks_pvalue_asymptotic(D::Float64, n::Int)
	if D <= 0
		return 1.0
	end
	λn = (sqrt(n) + 0.12 + 0.11 / sqrt(n)) * D
	s = 0.0
	for k in 1:200
		term = (-1)^(k - 1) * exp(-2 * k^2 * λn^2)
		s += term
		if abs(term) < 1e-12
			break
		end
	end
	clamp(2 * s, 0.0, 1.0)
end

function monte_carlo_binomial_lower_pvalue(k_obs::Int, n::Int, p::Float64; n_sim::Int = 20000)
	if n <= 0
		return 1.0
	end
	counts = Vector{Int}(undef, n_sim)
	for s in 1:n_sim
		c = 0
		for _ in 1:n
			c += rand() < p
		end
		counts[s] = c
	end
	mean(counts .<= k_obs)
end

function main()
	city_pattern = length(ARGS) >= 1 ? ARGS[1] : "תל אביב"
	input_jld2 = length(ARGS) >= 2 ? ARGS[2] : "alerts.jld2"
	output_png = length(ARGS) >= 3 ? ARGS[3] : "tel_aviv_interval_distribution.png"
	start_date = length(ARGS) >= 4 ? Date(ARGS[4], dateformat"y-m-d") : nothing
	end_date = length(ARGS) >= 5 ? Date(ARGS[5], dateformat"y-m-d") : nothing
	refractory_start_min = length(ARGS) >= 6 ? parse(Float64, ARGS[6]) : 30.0
	refractory_end_min = length(ARGS) >= 7 ? parse(Float64, ARGS[7]) : 60.0
	scan_step_min = length(ARGS) >= 8 ? parse(Float64, ARGS[8]) : 10.0
	scan_window_width_min = length(ARGS) >= 9 ? parse(Float64, ARGS[9]) : 30.0

	contains_hebrew(s::AbstractString) = any(c -> '\u0590' <= c <= '\u05FF', s)
	reverse_text(s::AbstractString) = join(reverse(collect(s)))

	function reverse_hebrew_segments(s::String)
		replace(s, r"[\u0590-\u05FF\s]+" => x -> reverse_text(x))
	end
	title_city = contains_hebrew(city_pattern) ? reverse_hebrew_segments(city_pattern) : city_pattern

	if start_date !== nothing && end_date !== nothing && start_date > end_date
		error("start_date must be <= end_date")
	end
	if refractory_start_min < 0 || refractory_end_min <= refractory_start_min
		error("refractory window must satisfy 0 <= start < end")
	end
	if scan_step_min <= 0 || scan_window_width_min <= 0
		error("scan_step_min and scan_window_width_min must be > 0")
	end

	alerts = jldopen(input_jld2, "r") do file
		file["alerts"]
	end

	all_names = names(alerts)

	function get_col(df::DataFrame, name::String)
		if name in names(df)
			return df[!, name]
		elseif Symbol(name) in names(df)
			return df[!, Symbol(name)]
		else
			return nothing
		end
	end

	city_col = get_col(alerts, "city")
	date_col = get_col(alerts, "date")
	alert_id_col = get_col(alerts, "alert_id")

	if city_col === nothing || date_col === nothing || alert_id_col === nothing
		error("alerts table must contain columns: alert_id, city, date. Found: $(all_names)")
	end

	city_strings = String.(city_col)
	mask = occursin.(city_pattern, city_strings)

	alert_dates = Date.(date_col)
	if start_date !== nothing
		mask .&= alert_dates .>= start_date
	end
	if end_date !== nothing
		mask .&= alert_dates .<= end_date
	end

	city_alerts = DataFrame(
		alert_id = alert_id_col[mask],
		city = city_col[mask],
		date = date_col[mask],
	)

	if nrow(city_alerts) == 0
		error("No alerts found for city pattern/date range selection")
	end

	event_times = combine(groupby(city_alerts, :alert_id), :date => minimum => :date)
	sort!(event_times, :date)

	if nrow(event_times) < 2
		error("Need at least 2 alert events to compute intervals for city pattern: \"$city_pattern\"")
	end

	intervals = diff(event_times.date)
	interval_minutes = Float64.(Dates.value.(intervals) ./ 60000)
	positive_intervals = interval_minutes[interval_minutes .> 0]

	if isempty(positive_intervals)
		error("No positive inter-arrival times found after filtering")
	end

	mean_interval = mean(positive_intervals)
	median_interval = median(positive_intervals)
	λ_hat = 1 / mean_interval

	ks_D = ks_statistic_exponential(positive_intervals, λ_hat)
	ks_p = ks_pvalue_asymptotic(ks_D, length(positive_intervals))

	at_risk_mask = positive_intervals .>= refractory_start_min
	n_at_risk = count(at_risk_mask)
	if n_at_risk == 0
		error("No intervals survive until refractory_start_min=$(refractory_start_min) min")
	end

	events_in_window = count((positive_intervals .>= refractory_start_min) .& (positive_intervals .< refractory_end_min))
	empirical_conditional_prob = events_in_window / n_at_risk

	window_width = refractory_end_min - refractory_start_min
	poisson_conditional_prob = 1 - exp(-λ_hat * window_width)
	refractory_ratio = empirical_conditional_prob / poisson_conditional_prob
	refractory_lower_p = monte_carlo_binomial_lower_pvalue(events_in_window, n_at_risk, poisson_conditional_prob)

	max_scan_start = maximum(positive_intervals) - scan_window_width_min
	scan_starts = max_scan_start >= 0 ? collect(0.0:scan_step_min:max_scan_start) : Float64[]

	scan_df = DataFrame(
		window_start_min = Float64[],
		window_end_min = Float64[],
		n_at_risk = Int[],
		events_in_window = Int[],
		empirical_conditional_prob = Float64[],
		poisson_conditional_prob = Float64[],
		refractory_ratio = Float64[],
		lower_tail_pvalue = Float64[],
	)

	for start_min in scan_starts
		end_min = start_min + scan_window_width_min
		n_risk = count(positive_intervals .>= start_min)
		if n_risk == 0
			continue
		end

		k = count((positive_intervals .>= start_min) .& (positive_intervals .< end_min))
		emp_prob = k / n_risk
		pois_prob = 1 - exp(-λ_hat * (end_min - start_min))
		ratio = emp_prob / pois_prob
		p_lower = monte_carlo_binomial_lower_pvalue(k, n_risk, pois_prob)

		push!(scan_df, (
			start_min,
			end_min,
			n_risk,
			k,
			emp_prob,
			pois_prob,
			ratio,
			p_lower,
		))
	end

	x_min = minimum(positive_intervals)
	x_max = maximum(positive_intervals)
	x_grid = collect(range(x_min, x_max; length = 500))
	theory_cdf = 1 .- exp.(-λ_hat .* x_grid)

	x_sorted = sort(positive_intervals)
	n = length(x_sorted)
	empirical_cdf = collect(1:n) ./ n

	p_hist = histogram(
		log10.(positive_intervals);
		bins = 10,
		xlabel = "log10(Interval between alerts in minutes)",
		ylabel = "Count",
		title = "Alert interval distribution: $title_city",
		legend = false,
		alpha = 0.8,
		left_margin = 10mm,
		bottom_margin = 10mm,
	)

	p_cdf = plot(
		x_grid,
		theory_cdf;
		label = "Exponential CDF (Poisson)",
		xlabel = "Interval (minutes)",
		ylabel = "CDF",
		title = "KS D=$(round(ks_D, digits=4)), p=$(round(ks_p, digits=4))",
		linewidth = 2,
		left_margin = 8mm,
		bottom_margin = 10mm,
	)
	scatter!(p_cdf, x_sorted, empirical_cdf; label = "Empirical CDF", markersize = 2, alpha = 0.5)

	p = plot(p_hist, p_cdf; layout = (1, 2), size = (1300, 560))

	savefig(p, output_png)

	root, ext = splitext(output_png)
	scan_table_path = "$(root)_refractory_scan.csv"
	scan_plot_path = "$(root)_refractory_scan$(ext)"
	CSV.write(scan_table_path, scan_df)

	if nrow(scan_df) > 0
		p_ratio = plot(
			scan_df.window_start_min,
			scan_df.refractory_ratio;
			xlabel = "Window start (minutes)",
			ylabel = "Empirical / Poisson",
			title = "Refractory ratio scan: $title_city",
			marker = :circle,
			markersize = 3,
			linewidth = 1.5,
			label = "Ratio",
			left_margin = 10mm,
			bottom_margin = 10mm,
		)
		hline!(p_ratio, [1.0]; linestyle = :dash, label = "Poisson baseline")

		p_pval = plot(
			scan_df.window_start_min,
			scan_df.lower_tail_pvalue;
			xlabel = "Window start (minutes)",
			ylabel = "Lower-tail p-value",
			title = "Lower-than-Poisson p-values",
			marker = :circle,
			markersize = 3,
			linewidth = 1.5,
			label = "p-value",
			left_margin = 10mm,
			bottom_margin = 10mm,
		)
		hline!(p_pval, [0.05]; linestyle = :dash, label = "0.05")

		p_scan = plot(p_ratio, p_pval; layout = (2, 1), size = (900, 800))
		savefig(p_scan, scan_plot_path)
	end

	println("City pattern: $city_pattern")
	println("Date range: $(start_date === nothing ? "(none)" : string(start_date)) to $(end_date === nothing ? "(none)" : string(end_date))")
	println("Rows matched: $(nrow(city_alerts))")
	println("Unique events: $(nrow(event_times))")
	println("Intervals total: $(length(interval_minutes))")
	println("Intervals > 0: $(length(positive_intervals))")
	println("Mean interval (min): $(round(mean_interval, digits=2))")
	println("Median interval (min): $(round(median_interval, digits=2))")
	println("Poisson-process fit: λ = $(round(λ_hat, digits=6)) per minute")
	println("KS statistic D: $(round(ks_D, digits=6))")
	println("KS p-value (asymptotic): $(round(ks_p, digits=6))")
	println("Refractory window (min): [$(refractory_start_min), $(refractory_end_min))")
	println("At-risk intervals (X >= start): $n_at_risk")
	println("Events in refractory window: $events_in_window")
	println("Empirical conditional P(start <= X < end | X >= start): $(round(empirical_conditional_prob, digits=6))")
	println("Poisson conditional probability for same window: $(round(poisson_conditional_prob, digits=6))")
	println("Refractory ratio (empirical / poisson): $(round(refractory_ratio, digits=6))")
	println("One-sided p-value for lower-than-poisson (Monte Carlo): $(round(refractory_lower_p, digits=6))")
	println("Scan step (min): $(scan_step_min)")
	println("Scan window width (min): $(scan_window_width_min)")
	println("Scan windows evaluated: $(nrow(scan_df))")
	println("Saved scan table to: $scan_table_path")
	if nrow(scan_df) > 0
		println("Saved scan plot to: $scan_plot_path")
	else
		println("Scan plot skipped (no valid windows)")
	end
	println("Saved plot to: $output_png")
end

main()
