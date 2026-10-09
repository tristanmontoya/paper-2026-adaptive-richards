# Run all eight conservation studies or one selected study; see the root README.md.

module RichardsConservation

using HydroTrixi, Trixi, SciMLBase
using Dates, LinearAlgebra

include("study_config.jl")

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
# Solver options and output names for each formulation and transfer choice
const CASES = (
    (name = :mixed, formulation = :mixed, transfer = :water_content,
     overrides = (;)),
    (name = :mixed_pressure_head_transfer, formulation = :mixed,
     transfer = :pressure_head, overrides = (; form = MixedForm(; transfer_state = true))),
    (name = :pressure_head_pressure_head_transfer, formulation = :pressure_head,
     transfer = :pressure_head,
     overrides = (; form = PressureHeadForm(),
                   error_control_block = state_variable_block,
                   error_control_mapping = water_content)),
    (name = :pressure_head_water_content_transfer, formulation = :pressure_head,
     transfer = :water_content,
     overrides = (; form = PressureHeadForm(; transfer_variables = water_content),
                   error_control_block = state_variable_block,
                   error_control_mapping = water_content)),
)
const STUDIES = [(benchmark = benchmark.name, formulation = config.formulation,
                  transfer = config.transfer)
                 for benchmark in BENCHMARKS for config in CASES]

# Keep the elixirs' global assignments separate from the conservation driver.
module InfiltrationElixir
result = nothing
end

function study_configuration(study)
    return first(config for config in CASES
                 if config.formulation == study.formulation &&
                    config.transfer == study.transfer)
end

function study_name(study)
    return data_stem(study.benchmark, study_configuration(study).name)
end

function save_snapshots(result, study, tolerance, snapshot_directory, final_time)
    component = study.formulation == :pressure_head ? 1 : 2
    for (stage, time) in ((:half, final_time / 2), (:full, final_time))
        index = findfirst(==(time), result.sol.t)
        data = solution_data_1d(result; index, component)
        path = joinpath(snapshot_directory,
                        snapshot_filename(study_name(study),
                                                                     tolerance, stage))
        open(path, "w") do io
            println(io, "#time_s depth_m pressure_head_m mesh_edge_m")
            for i in eachindex(data.x)
                mesh_edge = i <= length(data.mesh_vertices_x) ? data.mesh_vertices_x[i] : NaN
                println(io, join((data.time, data.x[i], data.values[i], mesh_edge), ' '))
            end
        end
    end
    return nothing
end

function solve_case(study, tolerance, data_directory, snapshot_directory)
    configuration = study_configuration(study)
    name = study_name(study)
    partial_path = joinpath(data_directory,
                            ".$(name)_rtol$(tolerance.tag).partial")
    final_path = joinpath(data_directory,
                          data_filename(name, tolerance))
    elixir = joinpath(ROOT, "code", "elixirs",
                      "elixir_richards_celia_$(study.benchmark).jl")
    final_time = first(benchmark.final_time for benchmark in BENCHMARKS
                       if benchmark.name == study.benchmark)

    amr_options = study.benchmark == :haverkamp ? (; amr = true) : (;)
    redirect_stdout(devnull) do
        Trixi.trixi_include(InfiltrationElixir, elixir;
                            amr_options..., configuration.overrides...,
                            reltol = tolerance.value,
                            analysis_interval = 1, save_analysis = true,
                            output_directory = data_directory,
                            analysis_filename = basename(partial_path),
                            dense = false, saveat = [final_time / 2, final_time],
                            save_mesh_history = true)
    end
    result = InfiltrationElixir.result
    sol = result.sol
    if !SciMLBase.successful_retcode(sol)
        error("Conservation solve failed; partial data retained at $partial_path")
    end

    history = accepted_step_history(partial_path; include_initial = true)
    _, biases = mass_bias_history(partial_path)
    open(final_path, "w") do io
        println(io, "#accepted_step time_s dt_s mass_bias_m")
        for row in zip(history.steps, history.times, history.dts, biases)
            println(io, join(row, ' '))
        end
    end
    rm(partial_path)
    save_snapshots(result, study, tolerance, snapshot_directory, final_time)
    return final_path
end

function run_conservation(studies = STUDIES; output_directory = nothing)
    if !all(study -> study in STUDIES, studies)
        error("Unknown conservation study")
    end
    BLAS.set_num_threads(1)

    id = Dates.format(now(UTC), "yyyymmddTHHMMSSsssZ")
    output = isnothing(output_directory) ?
             joinpath(ROOT, "results", "richards_conservation", id) :
             abspath(output_directory)
    mkpath(dirname(output))
    mkdir(output)
    mkdir(joinpath(output, "data"))
    mkdir(joinpath(output, "snapshots"))

    for study in studies, tolerance in RELTOLS
        println("Running $(study_name(study)) with reltol=$(tolerance.tag)")
        solve_case(study, tolerance, joinpath(output, "data"),
                   joinpath(output, "snapshots"))
    end
    println("Results: $output")
    return output
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    if !(length(ARGS) in (0, 3))
        error("Usage: julia --project=. $(@__FILE__) " *
              "[haverkamp|new_mexico mixed|pressure_head water_content|pressure_head]")
    end
    studies = if isempty(ARGS)
        RichardsConservation.STUDIES
    else
        [(benchmark = Symbol(ARGS[1]), formulation = Symbol(ARGS[2]),
          transfer = Symbol(ARGS[3]))]
    end
    RichardsConservation.run_conservation(studies)
end
