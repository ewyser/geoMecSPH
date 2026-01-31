"""
    superInc(lists::Vector{String}; root::String="", tree::Bool=false) -> Vector{String}

Recursively includes all `.jl` files in the given directories and their subdirectories.

# Arguments
- `lists::Vector{String}`: List of directory names to include.
- `root::String=""`: Root directory for the modules.
- `tree::Bool=false`: If true, displays a tree structure of included files.

# Returns
- `Vector{String}`: Messages summarizing the inclusion status for each directory.
"""
function superInc(lists::Vector{String}; root::String="", tree::Bool=false)
    sucess = ["superInc() jls parser:"]
    for (k, dir) ∈ enumerate(lists)
        dict = superDir(joinpath(root, dir))
        # Collect and include all .jl files in this subtree
        function collect_and_include_jls(d)
            files = String[]
            for (k, v) ∈ d
                if isa(v, Dict)
                    append!(files, collect_and_include_jls(v))
                elseif endswith(k, ".jl")
                    include(v)
                    push!(files, k)
                end
            end
            return files
        end
        jls_files = collect_and_include_jls(dict)
        push!(sucess, "✓ $(dir)")
    end
    return sucess
end

"""
    superDir(DIR::String) -> Dict

Recursively builds a nested Dict representing the directory structure and `.jl` files.

# Arguments
- `DIR::String`: Directory to search for Julia files and subdirectories.

# Returns
- `Dict`: Nested dictionary where keys are directory or file names.
"""
function superDir(DIR::String)
    d = Dict{String,Any}()
    for entry ∈ readdir(DIR; join=true)
        name = splitpath(entry)[end]
        if isdir(entry)
            d[name] = superDir(entry)
        elseif endswith(name, ".jl")
            d[name] = entry
        end
    end
    return d
end
