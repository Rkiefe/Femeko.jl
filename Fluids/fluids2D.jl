#=
    2D viscous fluid simulation
=#

include("../src/gmsh_wrapper.jl")
include("../src/FEM.jl")

# using GLMakie

function main(meshSize=0.0, localSize=0.0, showGmsh=false)

    # Simulation settings
    viscosity::Float64 = 1.0 

    # Create model
    gmsh.initialize()
    
    # Add an obstacle
    cells = []
    # id = addRectangle([0,0,0], [5, 3], cells)
    id = addDisk([-5,0,0], 1.0, cells)

    # Add a container
    box = addRectangle([0,0,0], [20, 5])
    # box = addDisk([0,0,0], 4)

    # Combine the geometries
    gmsh.model.occ.fragment(vcat(cells,[(2,box)]), [])
    gmsh.model.occ.synchronize()

    # Generate mesh
    mesh::MESH = Mesh2D(cells, meshSize, localSize, 2)

    println("\nNumber of elements ",size(mesh.t,2))
    println("Number of Inside elements ",length(mesh.InsideElements))
    println("Number of nodes ",size(mesh.p,2))
    println("Number of Inside nodes ",length(mesh.InsideNodes))
    println("Number of surface elements ",size(mesh.surfaceT,2))
    println("Mesh Order: ", mesh.order)

    # Run Gmsh GUI
    if showGmsh
       gmsh.fltk.run()
    end
    gmsh.fltk.finalize()

    # Get the number of mesh nodes (not counting the midpoints)
    nVertices = length(unique(vec(mesh.t[1:3,:])))

    # Global Stiffness matrix
    A = spzeros(mesh.nv, mesh.nv)

    # Local stiffness matrix
    Ak::Matrix{Float64} = zeros(36, mesh.nt) # 6 x 6
    temp::Matrix{Float64} = zeros(6,6)

    for k in 1:mesh.nt
        nds = @view mesh.t[:,k]

        temp .= 0
        for i in 1:length(nds)
            Si = quadraticBasis2D(mesh.p, nds, nds[i])
            
            for j in i:length(nds)
                Sj = quadraticBasis2D(mesh.p, nds, nds[j])

                # 6 node quadrature
                aux::Float64 = 0.0
                for n in 1:6
                    dxi::Float64 = Si[2] + 2*Si[4]*mesh.p[1,nds[n]] + Si[5]*mesh.p[2,nds[n]]
                    dyi::Float64 = Si[3] + Si[5]*mesh.p[1,nds[n]] + 2*Si[6]*mesh.p[2,nds[n]] 
                    dxj::Float64 = Sj[2] + 2*Sj[4]*mesh.p[1,nds[n]] + Sj[5]*mesh.p[2,nds[n]]
                    dyj::Float64 = Sj[3] + Sj[5]*mesh.p[1,nds[n]] + 2*Sj[6]*mesh.p[2,nds[n]]
                    aux += dxi*dxj + dyi*dyj
                end 
                aux /= 6

                temp[i,j] = aux*mesh.VE[k]
                temp[j,i] = temp[i,j] # It is symmetric
            end
        end

        Ak[:,k] = vec(temp)

    end # Local stiffness matrix

    # Update sparse global matrix
    n = 0
    for i in 1:6
        for j in 1:6
            n += 1
            A += sparse(mesh.t[i,:],mesh.t[j,:],Ak[n,:],mesh.nv,mesh.nv)
        end
    end


    return

    # Testing quadratic basis functions
    nds = @view mesh.t[:, 1]
    S = quadraticBasis2D(mesh.p, nds, nds[3])

    for nd in nds
        x = mesh.p[1, nd]
        y = mesh.p[2, nd]

        f = S[1] + S[2]*x + S[3]*y + S[4]*x^2 + S[5]*x*y + S[6]*y^2
        println(f)
    end

    return

    # Element centroids
    centroids::Matrix{Float64} = zeros(2, mesh.nt)
    for k in 1:mesh.nt
        nds = mesh.t[1:3,k]
        centroids[1,k] = sum(mesh.p[1,nds])/3
        centroids[2,k] = sum(mesh.p[2,nds])/3
    end

    
    
    # Mass matrix
    Mlocal::Matrix{Float64} = 1/12 *[2 1 1;
                                     1 2 1;
                                     1 1 2]

    M = spzeros(mesh.nv,mesh.nv)
    Mk::Matrix{Float64} = zeros(9, mesh.nt)
    for k in 1:mesh.nt
        nds = @view mesh.t[:,k]
        Mk[:,k] = mesh.VE[k]*Mlocal[:];
    end

    # Update sparse global matrix
    n = 0
    for i in 1:3
        for j in 1:3
            n += 1
            M += sparse(mesh.t[i,:],mesh.t[j,:],Mk[n,:],mesh.nv,mesh.nv)
        end
    end

    # fig = Figure()
    # ax = Axis(fig[1, 1], aspect = DataAspect(), title="Heat simulation")
    # scatterPlot = scatter!(ax, 
    #     mesh.p[1,:],
    #     mesh.p[2,:],
    #     color = T, 
    #     colormap=:thermal, 
    #     colorrange = (minimum(T), maximum(T)),
    #     markersize=5) 

    # Colorbar(fig[1, 2], scatterPlot, label="Temperature")
    
    # display(fig)
end

main(0.0, 0.0, false)