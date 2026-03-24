using BSplineKit
using CategoricalArrays
using CSV
using DataFrames
using FileIO
using JLD2
using MAT
using StatsModels
using Unfold


function _getfieldvalue(obj, name::AbstractString, default=nothing)
    if obj === nothing
        return default
    end
    if obj isa Dict
        return get(obj, name, get(obj, Symbol(name), default))
    end
    sym = Symbol(name)
    if hasproperty(obj, sym)
        return getproperty(obj, sym)
    end
    return default
end


function _to_string(value)
    if value === nothing
        return nothing
    elseif value isa AbstractString
        return value
    elseif value isa AbstractVector && length(value) > 0
        return _to_string(value[1])
    else
        return string(value)
    end
end


function _to_float(value)
    if value === nothing
        return missing
    end
    try
        return Float64(value)
    catch
        return missing
    end
end


function _load_eeg(set_path::AbstractString)
    mat = MAT.matread(set_path)
    eeg = haskey(mat, "EEG") ? mat["EEG"] : nothing
    if eeg !== nothing
        return eeg
    end
    eeg = Dict{String, Any}()
    for (key, value) in mat
        if key in ("__header__", "__version__", "__globals__")
            continue
        end
        eeg[key] = value
    end
    if isempty(eeg)
        error("No EEG struct found in $set_path")
    end
    return eeg
end


function _load_eeg_data(eeg, set_path::AbstractString)
    data = _getfieldvalue(eeg, "data")
    nbchan = Int(round(_getfieldvalue(eeg, "nbchan", 0)))
    pnts = Int(round(_getfieldvalue(eeg, "pnts", 0)))
    trials = Int(round(_getfieldvalue(eeg, "trials", 1)))
    if data isa AbstractString
        base = _getfieldvalue(eeg, "filepath", dirname(set_path))
        datfile = _getfieldvalue(eeg, "datfile", data)
        fdt_path = joinpath(base, datfile)
        nvals = nbchan * pnts * trials
        raw = Vector{Float32}(undef, nvals)
        open(fdt_path, "r") do io
            read!(io, raw)
        end
        data_arr = reshape(raw, nbchan, pnts * trials)
    else
        data_arr = Array(data)
        if ndims(data_arr) == 3
            data_arr = reshape(data_arr, size(data_arr, 1), size(data_arr, 2) * size(data_arr, 3))
        end
    end
    return Float64.(data_arr)
end


function _as_event_list(events)
    if events === nothing
        return Any[]
    end
    if events isa AbstractVector
        return collect(events)
    end
    return [events]
end


function _events_to_dataframe(events)
    rows = Vector{NamedTuple}()
    for ev in events
        push!(
            rows,
            (
                type = _to_string(_getfieldvalue(ev, "type")),
                latency = _to_float(_getfieldvalue(ev, "latency")),
                condition = _to_string(_getfieldvalue(ev, "condition")),
                trial = _to_float(_getfieldvalue(ev, "trial")),
                trial_time = _to_float(_getfieldvalue(ev, "trial_time")),
                saccAmpl = _to_float(_getfieldvalue(ev, "saccAmpl")),
                fix_at = _to_string(_getfieldvalue(ev, "fix_at")),
                set = _to_string(_getfieldvalue(ev, "set")),
                pattern = _to_string(_getfieldvalue(ev, "pattern")),
            )
        )
    end
    df = DataFrame(rows)
    return df
end


function _set_categorical_levels!(df::DataFrame)
    if :condition in names(df)
        df.condition = categorical(df.condition)
        levels!(df.condition, ["no_conflict", "conflict"]; allowmissing=true)
    end
    if :fix_at in names(df)
        df.fix_at = categorical(df.fix_at)
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
        amplitudeThreshold::Real=150,
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
    basis = firbasis(τ = (-0.5, 1.5), sfreq = srate)
    if !isdefined(Unfold, :spl)
        error("Unfold.spl not available. Ensure BSplineKit is installed and Unfold is loaded with spline support.")
    end
    f_prev = @formula(0 ~ 1)
    f_next = @formula(0 ~ 1)
    f_fix = @formula(0 ~ 1 + condition * fix_at * trial_time + spl(saccAmpl, 5))
    f_d = @formula(0 ~ 1 + condition + trial)
    f_n = @formula(0 ~ 1 + condition * trial_time + trial)
    return (
        [
            :prev => (f_prev, basis),
            :next => (f_next, basis),
            :fixation => (f_fix, basis),
            :D => (f_d, basis),
            :N => (f_n, basis),
        ],
        basis
    )
end


function _export_beta_csv(model, out_csv::AbstractString)
    try
        tbl = Unfold.coeftable(model)
        CSV.write(out_csv, tbl)
    catch err
        @warn "Unable to export beta CSV: $err"
    end
end


function _export_matlab_ufresult(model, basis, chanlocs, out_mat::AbstractString)
    coefs = coef(model)
    beta = Float64.(coefs)
    ufresult = Dict(
        "beta" => beta,
        "times" => collect(basis.times),
        "chanlocs" => chanlocs,
    )
    MAT.matwrite(out_mat, Dict("ufresult" => ufresult))
end


function run_unfold_step_g(set_path::AbstractString, out_dir::AbstractString, subject_id::AbstractString)
    eeg = _load_eeg(set_path)
    data = _load_eeg_data(eeg, set_path)
    srate = Float64(_getfieldvalue(eeg, "srate", 0))
    events_raw = _getfieldvalue(eeg, "event")
    events = _as_event_list(events_raw)
    df = _events_to_dataframe(events)
    _set_categorical_levels!(df)
    df = _filter_missing_events(df)
    df.type = Symbol.(coalesce.(df.type, ""))

    boundary_latencies = collect(skipmissing(df.latency[df.type .== :boundary]))
    winrej = _continuous_artifact_detect(
        data,
        srate;
        amplitudeThreshold=150,
        windowsize=2000,
        stepsize=100,
        channels=collect(1:size(data, 1)),
        boundary_latencies=boundary_latencies
    )
    df = _apply_artifact_exclusion!(data, df, winrej, srate, (-0.5, 1.5))

    design, basis = _build_design(srate)
    contrasts = Dict(:condition => StatsModels.DummyCoding(base = "no_conflict"))
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
    FileIO.save(out_model, model; compress = true)

    out_csv = joinpath(out_dir, string(subject_id, "_beta_dc.csv"))
    _export_beta_csv(model, out_csv)

    out_mat = joinpath(out_dir, string(subject_id, "_ufresult.mat"))
    chanlocs = _getfieldvalue(eeg, "chanlocs")
    _export_matlab_ufresult(model, basis, chanlocs, out_mat)

    return out_model
end
