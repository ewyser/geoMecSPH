export setup_paths

"""
    setup_paths(base_path::String) -> String

Create output directory structure for simulation results.

# Arguments
- `base_path`: Base directory path (default: "./out/")

# Returns
- `path`: Output directory path
"""
function setup_paths(base_path::String="./out/")
    if !isdir(base_path)
        mkpath(base_path)
        @info "Created output directory: $(base_path)"
    end
    return base_path
end
