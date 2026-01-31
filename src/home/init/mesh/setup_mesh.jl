export setup_mesh

"""
    setup_mesh(nel, lx, lz, typeD) -> (Mesh, Boundary)

Initialize computational mesh and boundary conditions.

# Arguments
- `nel`: Number of elements [nelx, nelz]
- `lx`: Domain length in x-direction
- `lz`: Domain length in z-direction  
- `typeD`: Data type (Float64 or Float32)

# Returns
- `meD`: Mesh structure
- `bc`: Boundary condition structure
"""
function setup_mesh(nel, lx, lz, typeD)
    # number of nodes
    nno = [nel + 1, nel + 1, (nel + 1) * (nel + 1)]
    nn = 16  # nodes per element
    
    # domain dimensions
    L = [lx, lz]
    h = L ./ nel
    
    # node coordinates
    xn = collect(typeD, LinRange(0, lx, nno[1]))
    zn = collect(typeD, LinRange(0, lz, nno[2]))
    
    xn = (xn' .* ones(typeD, nno[2], 1))
    zn = (ones(typeD, nno[1], 1)' .* zn)
    xn = vec(xn)
    zn = vec(zn)
    
    # element-to-node topology
    e2n = e2N(nno, nel, nn)
    
    # boundary conditions
    xB = [minimum(xn) + 2 * h[1], maximum(xn) - 2 * h[1], 0.0, Inf]
    bcx = vcat(findall(x -> x <= xB[1], xn), findall(x -> x >= xB[2], xn))
    bcz = findall(x -> x <= xB[3], zn)
    
    # create structures
    meD = Mesh(nel, nno, nn, L, h, xn, zn, e2n, xB)
    bc = Boundary(bcx, bcz)
    
    return (meD, bc)
end

"""
    e2N(nno, nel, nn) -> Array{Int64}

Generate element-to-node connectivity matrix.

# Arguments
- `nno`: Number of nodes
- `nel`: Number of elements
- `nn`: Nodes per element

# Returns
- `e2n`: Element-to-node connectivity array
"""
function e2N(nno, nel, nn)
    gnum = reverse(reshape(1:(nno[3]), nno[2], nno[1]), dims=1)
    e2n = zeros(nel[3], nn)
    iel = 1
    
    for i in 1:nel[1]  # nelx
        for j in 1:nel[2]  # nelz
            if (i > 1 && i < nel[1] && j > 1 && j < nel[2])
                e2n[iel, 1] = gnum[j-1, i-1]
                e2n[iel, 2] = gnum[j-0, i-1]
                e2n[iel, 3] = gnum[j+1, i-1]
                e2n[iel, 4] = gnum[j+2, i-1]

                e2n[iel, 5] = gnum[j-1, i]
                e2n[iel, 6] = gnum[j-0, i]
                e2n[iel, 7] = gnum[j+1, i]
                e2n[iel, 8] = gnum[j+2, i]

                e2n[iel, 9] = gnum[j-1, i+1]
                e2n[iel, 10] = gnum[j-0, i+1]
                e2n[iel, 11] = gnum[j+1, i+1]
                e2n[iel, 12] = gnum[j+2, i+1]

                e2n[iel, 13] = gnum[j-1, i+2]
                e2n[iel, 14] = gnum[j-0, i+2]
                e2n[iel, 15] = gnum[j+1, i+2]
                e2n[iel, 16] = gnum[j+2, i+2]
            end
            iel = iel + 1
        end
    end
    
    return convert(Array{Int64}, e2n)
end
