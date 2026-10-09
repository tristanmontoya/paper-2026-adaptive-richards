# Replot one or more saved conservation-study directories; never runs a simulation.

module RichardsConservationPlots

using CairoMakie, HydroTrixi, LaTeXStrings
using Dates, Printf
using DelimitedFiles

include("study_config.jl")

const RTOL_FONT = assetpath("fonts", "DejaVuSansMono.ttf")
const RTOL_LABELS = [rich(rich("rtol", font = RTOL_FONT), " = 10",
                          superscript(String(last(split(tolerance.tag, 'e')))))
                     for tolerance in RELTOLS]
const SOLUTION_YLIMS = Dict(:haverkamp => (-0.65, -0.15),
                            :new_mexico => (-10.5, -0.5))
const GROUPS = ((name = :mixed, cases = ((name = :mixed, transfer = :water_content),
                                         (name = :mixed_pressure_head_transfer,
                                          transfer = :pressure_head))),
                (name = :pressure_head,
                 cases = ((name = :pressure_head_pressure_head_transfer,
                           transfer = :pressure_head),
                          (name = :pressure_head_water_content_transfer,
                           transfer = :water_content))))

function collect_tables(directories)
    files = Dict{String, String}()
    for directory in directories
        data_directory = joinpath(abspath(directory), "data")
        for path in sort(readdir(data_directory; join = true))
            if !endswith(path, ".dat")
                continue
            end
            name = basename(path)
            if haskey(files, name)
                error("Duplicate conservation table: $name")
            end
            files[name] = path
        end
    end

    tables = []
    expected = Set{String}()
    for benchmark in BENCHMARKS, group in GROUPS
        transfer_tables = []
        for case in group.cases
            study_name = data_stem(benchmark.name, case.name)
            paths = [get(files, data_filename(study_name, tolerance), nothing)
                     for tolerance in RELTOLS]
            present = .!isnothing.(paths)
            if any(present)
                paths = String.(paths[present])
                for path in paths
                    push!(expected, basename(path))
                end
                # Preserve tolerance identities when intermediate files are absent
                push!(transfer_tables, (; stem = study_name, transfer = case.transfer,
                                        paths, tolerance_indices = findall(present)))
            end
        end
        if !isempty(transfer_tables)
            push!(tables, (; name = data_stem(benchmark.name, group.name),
                           final_time = benchmark.final_time,
                           solution_ylims = SOLUTION_YLIMS[benchmark.name], transfer_tables))
        end
    end
    if isempty(tables)
        error("No conservation studies found")
    end
    unexpected = setdiff(Set(keys(files)), expected)
    if !isempty(unexpected)
        error("Unexpected conservation tables: $(join(sort!(collect(unexpected)), ", "))")
    end

    return tables
end

function plot_study(table, output_directory, bias_axis)
    time_ticks = if table.final_time == 360.0
        collect(0.0:60.0:360.0)
    else
        collect(0.0:20_000.0:80_000.0)
    end
    tick_labels = [@sprintf("%.0f", time) for time in time_ticks]
    legend_position = (:right, :top)
    paths = String[]
    labels = Any[]
    colors = Any[]
    linestyles = Symbol[]
    palette = Makie.wong_colors()
    water_content_tolerances = [i for transfer_table in table.transfer_tables
                               if transfer_table.transfer === :water_content
                               for i in transfer_table.tolerance_indices]
    for transfer_table in table.transfer_tables
        linestyle = if transfer_table.transfer === :water_content
            :solid
        else
            :dash
        end
        for (i, path) in zip(transfer_table.tolerance_indices, transfer_table.paths)
            # Label each tolerance once, preferring its water-content curve
            label = if transfer_table.transfer === :water_content ||
                       !(i in water_content_tolerances)
                RTOL_LABELS[i]
            else
                nothing
            end
            push!(paths, path)
            push!(labels, label)
            push!(colors, palette[i])
            push!(linestyles, linestyle)
        end
    end

    bias_path = joinpath(output_directory, "$(table.name)_mass_bias.pdf")
    bias_fig = plot_mass_bias_magnitude(paths; output_path = bias_path,
                                        time_column = "time_s",
                                        mass_balance_column = "mass_bias_m", labels,
                                        colors, linestyles, legend_position,
                                        xticks = (time_ticks, tick_labels),
                                        yticks = bias_axis.ticks,
                                        xlims = (0.0, table.final_time),
                                        ylims = bias_axis.limits)
    # Apply the manuscript terminology to the mass-bias plotter's fixed labels
    ax = bias_fig.content[1]
    ax.xlabel = L"Time $t$ (s)"
    ax.ylabel = L"Absolute bias error $|\epsilon_{\mathrm{b}}(t)|$ (m)"
    save(bias_path, bias_fig)
    output_paths = [bias_path]

    # The manuscript compares time steps for water-content solution transfer only
    water_content_paths = [path for transfer_table in table.transfer_tables
                           if transfer_table.transfer === :water_content
                           for path in transfer_table.paths]
    if !isempty(water_content_paths)
        step_path = joinpath(output_directory, "$(table.name)_time_steps.pdf")
        plot_time_steps(water_content_paths; output_path = step_path,
                        step_column = "accepted_step", time_column = "time_s",
                        dt_column = "dt_s",
                        xlabel = L"Time $t$ (s)",
                        ylabel = L"Time step size $\Delta t$ (s)",
                        labels = RTOL_LABELS[water_content_tolerances],
                        colors = palette[water_content_tolerances],
                        legend_position = (:left, :top),
                        xticks = (time_ticks, tick_labels),
                        xlims = (0.0, table.final_time))
        push!(output_paths, step_path)
    end
    return output_paths
end

function plot_snapshots(table, output_directory)
    output_paths = String[]
    skipped = String[]
    for transfer_table in table.transfer_tables
        for (i, path) in zip(transfer_table.tolerance_indices, transfer_table.paths)
            for stage in (:half, :full)
                filename = snapshot_filename(transfer_table.stem, RELTOLS[i], stage)
                pdf_filename = replace(filename, ".dat" => ".pdf")
                snapshot_path = joinpath(dirname(dirname(path)), "snapshots", filename)
                if !isfile(snapshot_path)
                    push!(skipped, pdf_filename)
                    continue
                end
                values = readdlm(snapshot_path, Float64; comments = true)
                data = (; time = values[1, 1], x = values[:, 2], values = values[:, 3],
                         mesh_vertices_x = filter(isfinite, values[:, 4]))
                output_path = joinpath(output_directory, pdf_filename)
                plot_solution_1d(data; output_path,
                                  xlabel = L"Depth $z$ (m)",
                                  ylabel = L"Pressure head $\psi$ (m)",
                                  ylims = table.solution_ylims,
                                  show_element_boundaries = true)
                push!(output_paths, output_path)
            end
        end
    end
    if !isempty(skipped)
        println("Skipped $(length(skipped)) snapshot PDFs without saved tables:")
        foreach(path -> println("  $path"), skipped)
    end
    return output_paths
end

function plot_conservation(directories...)
    if isempty(directories)
        error("Supply at least one saved run directory")
    end
    tables = collect_tables(directories)
    paths = [path for table in tables for transfer_table in table.transfer_tables
             for path in transfer_table.paths]
    bias_axis = mass_bias_magnitude_axis(paths; time_column = "time_s",
                                         mass_balance_column = "mass_bias_m")

    id = Dates.format(now(UTC), "yyyymmddTHHMMSSsssZ")
    output_directory = joinpath(abspath(first(directories)), "figures", id)
    mkpath(dirname(output_directory))
    mkdir(output_directory)

    paths = String[]
    for table in tables
        append!(paths, plot_study(table, output_directory, bias_axis))
        append!(paths, plot_snapshots(table, output_directory))
    end
    println("Figures: $output_directory")
    return paths
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    if isempty(ARGS)
        error("Usage: julia --project=. $(@__FILE__) RUN_DIRECTORY [RUN_DIRECTORY ...]")
    end
    RichardsConservationPlots.plot_conservation(ARGS...)
end
