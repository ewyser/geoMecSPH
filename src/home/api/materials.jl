export elasticity_matrix

"""
    elasticity_matrix(E, ν) -> (K, G, Del)

Compute elastic material parameters and stiffness matrix.

# Arguments
- `E`: Young's modulus [Pa]
- `ν`: Poisson's ratio [-]

# Returns
- `K`: Bulk modulus [Pa]
- `G`: Shear modulus [Pa]
- `Del`: 4x4 elastic stiffness matrix [Pa]
"""
function elasticity_matrix(E, ν)
    # Compute elastic moduli
    G = E / (2.0 * (1.0 + ν))    # Shear modulus
    K = E / (3.0 * (1.0 - 2.0 * ν))  # Bulk modulus
    
    # Construct stiffness matrix (plane strain)
    Del = [K + 4/3*G  K - 2/3*G  K - 2/3*G  0.0;
           K - 2/3*G  K + 4/3*G  K - 2/3*G  0.0;
           K - 2/3*G  K - 2/3*G  K + 4/3*G  0.0;
           0.0        0.0        0.0        2*G]
    
    return (K, G, Del)
end
