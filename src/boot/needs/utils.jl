"""
    welcome_log() -> Nothing

Display welcome message for geoMecSPH package.
"""
function welcome_log()
    printstyled("┌ Welcome to geoMecSPH 👻\n", color=:green, bold=true)
    printstyled("│", color=:green, bold=true); println(" Geomechanical SPH/MPM solver")
    printstyled("│", color=:green, bold=true); println(" New comer ? Try:")
    printstyled("│", color=:green, bold=true); println("   using geoMecSPH")
    printstyled("└", color=:green, bold=true); println("   sim()")
    return nothing
end

"""
    exit_log(msg::String) -> Nothing

Display exit message.
"""
function exit_log(msg::String)
    @info msg
    return nothing
end
