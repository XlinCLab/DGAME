using BSplineKit
using CategoricalArrays
using CSV
using DataFrames
using JLD2
using JSON3
using StatsModels
using Unfold


function _set_categorical_levels!(df::DataFrame)
    if :condition in names(df)
        df.condition = categorical(df.condition)
        levels!(df.condition, ["no_conflict", "conflict"]; allowmissing=true)
    end
    if :fix_at in names(df)
        df.fix_at = categorical(df.fix_at)
        levels!(df.fix_at, ["elsewhere", "other", "target"]; allowmissing=true)
    end
    return df
end


function _filter_missing_events(df::DataFrame)
    required = Dict(
        "prev" => Symbol[],
        "next" => Symbol[],
        "fixation" => [:condition, :fix_at, :trial_time, :saccAmpl],
        "D" => [:condition, :trial],
        "N" => [:condition, :trial_time, :trial],
    )
    keep = trues(nrow(df))
    for (idx, row) in enumerate(eachrow(df))
        t = row.type
        cols = get(required, t, Symbol[])
        for col in cols
            if ismissing(row[col])
                keep[idx] = false
                break
            end
        end
    end
    return df[keep, :]
end


function _merge_overlaps(winrej::Vector{Tuple{Int,Int}})
    if isempty(winrej)
        return winrej
    end
    sorted = sort(winrej, by=x -> x[1])
    merged = [sorted[1]]
    for (start, stop) in sorted[2:end]
        last_start, last_stop = merged[end]
        if start <= last_stop
            merged[end] = (min(last_start, start), max(last_stop, stop))
        else
            push!(merged, (start, stop))
        end
    end
    return merged
end


function _continuous_artifact_detect(data::AbstractMatrix, srate::Real;
        amplitudeThreshold::Real=150, # microvolts
        windowsize::Real=2000,
        stepsize::Real=100,
        channels::AbstractVector=collect(1:size(data, 1)),
        boundary_latencies::AbstractVector=Int[])
    winpnts = floor(Int, windowsize * srate / 1000)
    stepnts = floor(Int, stepsize * srate / 1000)
    pnts = size(data, 2)
    latebound = [0, pnts]
    if !isempty(boundary_latencies)
        lb = [0; round.(Int, boundary_latencies .- 0.5); pnts]
        lb = clamp.(lb, 0, pnts)
        latebound = sort(unique(lb))
    end

    winrej = Vector{Tuple{Int,Int}}()
    for q in 1:(length(latebound) - 1)
        bp1 = latebound[q] + 1
        bp2 = latebound[q + 1]
        j = bp1
        while j <= bp2 - (winpnts - 1)
            t1 = j + 1
            t2 = j + winpnts - 1
            bad = false
            for ch in channels
                segment = @view data[ch, t1:t2]
                vdiff = maximum(segment) - minimum(segment)
                if abs(vdiff) > amplitudeThreshold
                    bad = true
                    break
                end
            end
            if bad
                push!(winrej, (t1, t2))
            end
            j += stepnts
        end
    end
    return _merge_overlaps(winrej)
end


function _apply_artifact_exclusion!(
        data::AbstractMatrix,
        events::DataFrame,
        winrej::Vector{Tuple{Int,Int}},
        srate::Real,
        timelimits::Tuple{Real,Real}
    )
    if isempty(winrej)
        return events
    end
    pnts = size(data, 2)
    bad_mask = falses(pnts)
    for (start, stop) in winrej
        start = max(start, 1)
        stop = min(stop, pnts)
        bad_mask[start:stop] .= true
    end
    data[:, bad_mask] .= 0.0

    keep = trues(nrow(events))
    pre = Int(round(timelimits[1] * srate))
    post = Int(round(timelimits[2] * srate))
    for (idx, row) in enumerate(eachrow(events))
        if ismissing(row.latency)
            keep[idx] = false
            continue
        end
        lat = Int(round(row.latency))
        start = max(lat + pre, 1)
        stop = min(lat + post, pnts)
        if start <= stop && any(bad_mask[start:stop])
            keep[idx] = false
        end
    end
    return events[keep, :]
end


function _build_design(srate::Real)

    basis_for(name) = firbasis(τ = (-0.5, 1.5), sfreq = srate, name = String(name))

    f_prev = @formula(0 ~ 1)
    f_next = @formula(0 ~ 1)
    f_fix = @formula(0 ~ 1 + condition * fix_at * trial_time + spl(saccAmpl, 5))
    f_d = @formula(0 ~ 1 + condition + trial)
    f_n = @formula(0 ~ 1 + condition * trial_time + trial)

    return [
        "prev"     => (f_prev, basis_for("prev")),
        "next"     => (f_next, basis_for("next")),
        "fixation" => (f_fix,  basis_for("fixation")),
        "D"        => (f_d,    basis_for("D")),
        "N"        => (f_n,    basis_for("N")),
    ]
end


function _export_beta_csv(model, out_csv::AbstractString)
    try
        tbl = Unfold.coeftable(model)
        CSV.write(out_csv, tbl;
            transform = (col, val) -> something(val, missing)
        )
    catch err
        @warn "Unable to export beta CSV: $err"
    end
end


function _export_ufresult_struct(
        out_path::AbstractString,
        model,
        fir_basis,
        chan_names::Vector{String},
        event_order::Vector{String},
    )
    # Extract coefficients from the fitted Unfold model
    coefs = coef(model)
    beta = Float64.(coefs)  # [channels × times × betas]

    # Extract the time axis from the FIRBasis object directly
    times = collect(fir_basis.times)

    # Minimal chanlocs struct array compatible with downstream reconstruction logic
    chanlocs_struct = [Dict("labels" => ch) for ch in chan_names]

    # Build explicit coefficient metadata so downstream reconstruction can refer to
    # coefficients by event-qualified name instead ordered position
    # Names specified as "event::coefficient" because coefficient labels such as
    # "(Intercept)" repeat across event blocks
    coef_table = Unfold.coeftable(model)
    coefnames_by_event = Dict{String, Vector{String}}()
    coef_order = String[]
    for ev in event_order
        ev_mask = coef_table.eventname .== ev
        ev_coefnames = String.(unique(coef_table.coefname[ev_mask]))
        coefnames_by_event[ev] = ev_coefnames
        append!(coef_order, ["$(ev)::$(name)" for name in ev_coefnames])
    end

    # Package into the same structure expected by downstream steps, plus explicit metadata
    ufresult = Dict(
        "beta" => beta,
        "times" => times,
        "chanlocs" => chanlocs_struct,
    )
    ufmeta = Dict(
        "event_order" => event_order,
        "coef_order" => coef_order,
        "coefnames_by_event" => coefnames_by_event,
    )

    # Persist a lightweight struct (betas/times/chanlocs) for downstream steps
    JLD2.@save out_path ufresult ufmeta
end


function run_unfold_step_g_from_arrays(
    data::AbstractMatrix,
    srate::Real,
    events_csv::AbstractString,
    chan_names::Vector{String},
    out_dir::AbstractString,
    subject_id::AbstractString,
)
    df = CSV.read(events_csv, DataFrame)
    _set_categorical_levels!(df)
    df = _filter_missing_events(df)
    df.type = Symbol.(coalesce.(df.type, ""))

    boundary_latencies = collect(skipmissing(df.latency[df.type .== :boundary]))
    winrej = _continuous_artifact_detect(
        data,
        srate;
        amplitudeThreshold=150, # microvolts
        windowsize=2000,
        stepsize=100,
        channels=collect(1:size(data, 1)),
        boundary_latencies=boundary_latencies
    )
    df = _apply_artifact_exclusion!(data, df, winrej, srate, (-0.5, 1.5))

    design = _build_design(srate)
    contrasts = Dict(
        # Reference level matches MATLAB Unfold's `cfgDesign.categorical` specification,
        # where no_conflict is the first (reference) level for condition.
        :condition => StatsModels.DummyCoding(base = "no_conflict"),
        # fix_at reference level must be set explicitly to match MATLAB's global
        # `cfgDesign.codingschema = 'reference'` which defaults to the first level
        # alphabetically — "elsewhere" — for all categorical predictors.
        :fix_at    => StatsModels.DummyCoding(base = "elsewhere"),
    )
    df.type = String.(df.type)
    model = fit(
        UnfoldModel,
        design,
        df,
        data;
        eventcolumn = :type,
        contrasts = contrasts,
    )

    mkpath(out_dir)
    out_model = joinpath(out_dir, string(subject_id, "_ufresult.jld2"))
    @info "Saving Unfold model to $out_model"
    JLD2.@save out_model model

    # Export beta CSV
    out_csv = joinpath(out_dir, string(subject_id, "_beta_dc.csv"))
    _export_beta_csv(model, out_csv)

    # Export UFRESULT-like struct — extract the FIRBasis from the design directly
    out_struct = joinpath(out_dir, string(subject_id, "_ufresult_struct.jld2"))
    fir_basis = design[findfirst(p -> p.first == "fixation", design)].second[2]
    event_order = [String(pair.first) for pair in design]
    _export_ufresult_struct(out_struct, model, fir_basis, chan_names, event_order)

    return out_model
end
