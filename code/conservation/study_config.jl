# Shared benchmark parameters and filenames for simulation and figure generation

const BENCHMARKS = ((name = :haverkamp, final_time = 360.0),
                    (name = :new_mexico, final_time = 86_400.0))
const RELTOLS = ((value = 1.0e-5, tag = "1e-5"),
                 (value = 1.0e-7, tag = "1e-7"),
                 (value = 1.0e-9, tag = "1e-9"))

data_stem(benchmark, case) = "richards_celia_$(benchmark)_$(case)"
data_filename(stem, tolerance) = "$(stem)_rtol$(tolerance.tag).dat"
snapshot_filename(stem, tolerance, stage) =
    "$(stem)_rtol$(tolerance.tag)_$(stage)_solution_mesh.dat"
