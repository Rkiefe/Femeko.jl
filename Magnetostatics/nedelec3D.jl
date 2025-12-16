#=
        Calculate the magnetostatic interaction between non-linear magnetic materials and a source field

    Using: 
        - Nedelec shape functions (instead of linear Lagrange)
        - The magnetostatic vector potential (instead of the scalar potential)
        - Magnetic Reluctance (instead of permeability)

    Why:
        - Should be much more stable for non-linear media
        - Can include Surface and Volume electric currents

    This example also includes the implementation of the Newton-Raphson method 
    with Nedelec shape elements

    In my experience, using a single F-P iteration and then immediately using
    the N-R method is best. Specially since there is now a line search method to
    reduce the residue of the next iteration
=#

include("../src/Femeko.jl")
include("../src/magneticProperties.jl")

using IterativeSolvers
using GLMakie

# Updates the Newton-Raphson iteration with a line search
# where the full step size is reduced to minimize the new residue
function lineSearch(u, du, A, RHS, minStep::Float64=1e-4)
    # Update the solution with a weight: u_new = u + alf*du
    alf = 1.0 # Initial weight (full step)

    # Residual before update of u
    residual_old = norm(RHS - A*u) 

    # New solution trial
    u_trial = u + du
    residual_new = norm(RHS - A*u_trial)
    
    # Reduce the step size until the new solution is more accurate than the old solution
    while residual_new > residual_old && alf > minStep
        alf *= 0.9 # Reduce the size of the weight by 10 %
        
        # New trial solution
        u_trial = u + alf*du
        residual_new = norm(RHS - A*u_trial)
    end
    
    # Update the solution
    if alf < minStep
        # Don't update u otherwise it will increase the residue
        println("Warning: N-R could not decrease the residue further")
    else
        u .= u_trial
    end

end # Adapt the N-R step size based on the residual

function main(meshSize=0.0, localSize=0.0, showGmsh=false)

    # vacuum magnetic permeability
    mu0::Float64 = pi*4e-7

    # Temperature
    T::Float64 = 300.0

    # Applied field | Tesla
    Bext::Vector{Float64} = [0.5, 
                             0.0,
                             0.0]

    # Convergence criteria
    picardDeviation::Float64 = 1.1*norm(Bext) # Only let 1 F-P iteration
    maxDeviation::Float64 = 1e-4
    maxAtt::Int32 = 100

    # Load Gd data
    # data = DATA()
    # loadMaterial( data,
    #              "Materials", # Folder with materials
    #              "Gd_MFT",    # Data folder of target material
    #              "Gd",        # Material name
    #              7.9,         # density g/cm3
    #              T)

    # Load Fe data
    data = DATA()
    data.density = 7.874 # g/cm3
    data.HofM = vec(readdlm("Materials/Pure_Iron_FEMM/H_Fe_extrap.dat"))  # A/m
    data.B = vec(readdlm("Materials/Pure_Iron_FEMM/B_Fe_extrap.dat"))     # T
    materialPermeability(data) # Get the permeability and its derivative

    # Prepare the spline for interpolation over the dataset
    spl = Spline1D(data.B, data.nu; bc = "nearest")

    # Derivative of the magnetic reluctance
    spl_dnu = Spline1D( data.B,   # x
                        gradient(data.B, data.nu) # y
                        ; bc = "nearest"
                       )

    # Create model
    gmsh.initialize()
    cells = []

    # Import step model
    # box = importCAD("../STEP_Models/Fennec_Fox.step", cells, true) # true -> make a bounding shell

    # Add magnetic geometry
    id = addCuboid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], cells)
    # id = addSphere([0.0, 0.0, 0.0], 0.5, cells)

    # Add a container
    box = addSphere([0.0, 0.0, 0.0], 5.0)

    # Combine the geometries
    shell_id, _ = unifyModel(cells, box)

    # Generate mesh
    mesh = Mesh(cells, meshSize, localSize, false, 2)

    # Edges (mid-points)
    edges = unique(mesh.t[5:10, :])
    ne = length(edges)

    # Map the global edge label to an ordered, local edge label
    # (the global edge label is not ordered with the node label)
    global2local_edge::Vector{Int32} = zeros(maximum(edges))
    for e in 1:ne
        edge = edges[e] # Global edge ID (its unique)
        global2local_edge[edge] = e
    end

    println("\nNumber of elements: ", mesh.nt)
    println("Number of Inside elements ", mesh.nInside)
    println("Number of edges: ", ne)
    println("Number of vertices: ", mesh.nv)
    println("Number of Inside nodes ", mesh.nInsideNodes)
    println("Number of surface elements ", mesh.ne)
    println("")

    # Run Gmsh GUI
    if showGmsh
        gmsh.option.setNumber("Mesh.VolumeFaces", 1)
        gmsh.option.setNumber("Mesh.Clip", 1)
        gmsh.option.setNumber("General.ClipWholeElements", 1)
        gmsh.fltk.run()
    end
    gmsh.fltk.finalize()

    # Local stiffness matrix
    Ak::Matrix{Float64} = nedelecStiffness(mesh)

    # Load vector
    RHS::Vector{Float64} = zeros(ne)
    for k in 1:mesh.nt
        nds = @view mesh.t[1:4, k] # Nodes of the linear volume element
        
        # Hat shape element for each of the 4 nodes
        hat::Matrix{Float64} = zeros(4, 4) # a,b,c,d for each node
        for i in 1:4
            hat[1, i], 
            hat[2, i], 
            hat[3, i], 
            hat[4, i] = abcd(mesh.p, nds, nds[i])
        end

        for ie in 1:6 # For each edge of the tetrahedron
            
            # Global edge label
            edge = mesh.t[4+ie, k]
            i = global2local_edge[edge] # Local edge label
            
            # Global node labels of the edge
            ndi, ndj = NodesFromLocalEdge(mesh, k, ie)

            # Length of edge
            edgeLength = norm(mesh.p[1:3, nds[ndj]] - mesh.p[1:3, nds[ndi]])

            # 1st order Lagrange basis function (hat function)
            _, bi, ci, di = hat[:, ndi]
            _, bj, cj, dj = hat[:, ndj]

            curlN = 2.0*edgeLength*cross([bi, ci, di], [bj, cj, dj])

            RHS[i] += mesh.VE[k]*dot(Bext, curlN)
        end
    
    end # Loop over the volume element labels

    # Relative magnetic reluctance
    nu::Vector{Float64} = ones(mesh.nt)

    # Prepare output
    u::Vector{Float64} = zeros(ne) # Magnetic vector potential

    Bfield::Matrix{Float64} = zeros(3, mesh.nt)
    B::Vector{Float64} = zeros(mesh.nt)
    Bold::Vector{Float64} = zeros(mesh.nt)

    att::Int32 = 0
    div::Float64 = Inf
    while att < maxAtt && div > picardDeviation
        att += 1
        Bold .= B

        # Global stiffness matrix
        A = spzeros(ne, ne)

        # Update sparse global stiffness matrix
        n = 0
        for i in 1:6
            edge1 = global2local_edge[mesh.t[4+i, :]]
            
            for j in 1:6
                n += 1

                edge2 = global2local_edge[mesh.t[4+j, :]]
                A += sparse(edge1, edge2, Ak[n,:].*nu, ne, ne)
            end
        end

        # Solve the magnetostatic vector potential
        u = cg(A, RHS) # ; verbose=true

        # Get the magnetic flux B
        Bfield .= 0.0
        for k in 1:mesh.nt
            nds = @view mesh.t[1:4, k] # Nodes of the linear volume element
            
            # Hat shape element for each of the 4 nodes
            hat::Matrix{Float64} = zeros(4, 4) # a,b,c,d for each node
            for i in 1:4
                hat[1, i], 
                hat[2, i], 
                hat[3, i], 
                hat[4, i] = abcd(mesh.p, nds, nds[i])
            end

            for ie in 1:6 # For each edge of the tetrahedron
                
                # Global edge label
                edge = mesh.t[4+ie, k]
                i = global2local_edge[edge] # Local edge label
                
                # Global node labels of the edge
                ndi, ndj = NodesFromLocalEdge(mesh, k, ie)

                # Length of edge
                edgeLength = norm(mesh.p[1:3, nds[ndj]] - mesh.p[1:3, nds[ndi]])

                # 1st order Lagrange basis function (hat function)
                _, bi, ci, di = hat[:, ndi]
                _, bj, cj, dj = hat[:, ndj]

                curlN = 2.0*edgeLength*cross([bi, ci, di], [bj, cj, dj])

                Bfield[:, k] += u[i].*curlN
            end
        end

        # Norm of flux field B
        B .= 0.0
        for k in 1:mesh.nt
            B[k] = norm(Bfield[:, k])
        end

        # Update reluctance
        nu[mesh.InsideElements] = spl(B[mesh.InsideElements])

        if any(x->!isfinite(x), nu)
            error("Nans in the interpolation")
        end

        div = maximum(abs.(B .- Bold))
        println(att, " , |B(n)-B(n-1)| = ", div)

    end

    println("Newton-Raphson iteration method")

    # Newton-Raphson
    dnu::Vector{Float64} = zeros(mesh.nt)
    dnu[mesh.InsideElements] = spl_dnu(B[mesh.InsideElements])
    while div > maxDeviation && att < maxAtt

        att += 1
        Bold .= B

        # Global sparse stiffness matrix
        A = spzeros(ne, ne)
        n = 0
        for i in 1:6
            edge1 = global2local_edge[mesh.t[4+i, :]]
            
            for j in 1:6
                n += 1

                edge2 = global2local_edge[mesh.t[4+j, :]]
                A += sparse(edge1, edge2, Ak[n, :].*nu, ne, ne)
            end
        end

        # Tangential stiffness matrix
        At = nedelecTangentialStiffness(mesh, global2local_edge, 
                                        dnu, ne, Bfield, B)

        # Correction to the magnetic scalar potential
        du = cg(A+At, RHS-A*u; verbose=false)

        # Update the vector potential
        # u .+= du # Direct update
        lineSearch(u, du, A, RHS) # Line search (more stable)
        residue = norm(RHS - A*u)

        # Get the magnetic flux B
        Bfield .= 0.0
        for k in 1:mesh.nt
            nds = @view mesh.t[1:4, k] # Nodes of the linear volume element
            
            # Hat shape element for each of the 4 nodes
            hat::Matrix{Float64} = zeros(4, 4) # a,b,c,d for each node
            for i in 1:4
                hat[1, i], 
                hat[2, i], 
                hat[3, i], 
                hat[4, i] = abcd(mesh.p, nds, nds[i])
            end

            for ie in 1:6 # For each edge of the tetrahedron
                
                # Global edge label
                edge = mesh.t[4+ie, k]
                i = global2local_edge[edge] # Local edge label
                
                # Global node labels of the edge
                ndi, ndj = NodesFromLocalEdge(mesh, k, ie)

                # Length of edge
                edgeLength = norm(mesh.p[1:3, nds[ndj]] - mesh.p[1:3, nds[ndi]])

                # 1st order Lagrange basis function (hat function)
                _, bi, ci, di = hat[:, ndi]
                _, bj, cj, dj = hat[:, ndj]

                curlN = 2.0*edgeLength*cross([bi, ci, di], [bj, cj, dj])

                Bfield[:, k] += u[i].*curlN
            end
        end

        # Norm of flux field B
        B .= 0.0
        for k in 1:mesh.nt
            B[k] = norm(Bfield[:, k])
        end

        # Update magnetic reluctance
        nu[mesh.InsideElements] .= spl(B[mesh.InsideElements])
        
        # d/dB nu
        dnu[mesh.InsideElements] = spl_dnu(B[mesh.InsideElements])

        if any(x->!isfinite(x), nu) || any(x->!isfinite(x), dnu)
            error("Nans in the interpolation")
        end

        # Check deviation from previous result
        div = maximum(abs.(B .- Bold))
        println(att, " , |B(n)-B(n-1)| = ", div, " , |y-Ax| = ", residue)

        # Check if Au = RHS by the residue
        # if residue < 1e-7 && att > 2
        #     println("condition Au = RHS reached within tol = ", residue)
        #     break
        # end

    end # Newton iteration

    Hfield = zeros(3, mesh.nt)
    H = zeros(mesh.nt)
    for k in 1:mesh.nt
        H[k] = nu[k]/mu0 * B[k]
        Hfield[:, k] = nu[k]/mu0 .* Bfield[:, k]
    end

    # Magnetization (A/m)
    chi::Vector{Float64} = 1.0./nu .- 1.0 # Magnetic susceptibility

    M::Vector{Float64} = chi.*H # Magnetization intensity

    # Magnetization vector field
    Mfield::Matrix{Float64} = zeros(3, mesh.nt)
    Mfield[1, :] = chi.*Hfield[1, :]
    Mfield[2, :] = chi.*Hfield[2, :]
    Mfield[3, :] = chi.*Hfield[3, :]

    # M (emu/g)
    M_emug = M./(data.density*1e3)
    Mfield_emug = Mfield./(data.density*1e3)

    # Plot results
    println("\nGenerating plots...")
    elements =  
                # 1:mesh.nt               # All elements of the mesh
                mesh.InsideElements     # Only the magnetic region

    # Element centroids
    centroids::Matrix{Float64} = zeros(3, mesh.nt)
    for k in 1:mesh.nt
        nds = mesh.t[:, k]
        centroids[:, k] = mean(mesh.p[:, nds], 2)
    end

    fig = Figure()
    ax = Axis3(fig[1,1], aspect = :data)
    
    graph = arrows3d!(  ax
                        , centroids[1, elements]
                        , centroids[2, elements]
                        , centroids[3, elements]
                        , mu0*Mfield[1, elements] 
                        , mu0*Mfield[2, elements]
                        , mu0*Mfield[3, elements]
                        , color = M_emug[elements]
                        , lengthscale = 0.1
                        , colormap = :CMRmap,  # :CMRmap :viridis :redsblues :turbo :rainbow
                    )

    # Add color bar
    Colorbar(fig[1, 2], graph, label = "M (emu/g)")

    wait(display(fig))
end

main(5.0, 0.1, false)
