#=
        Calculate the magnetostatic interaction between non-linear magnetic materials and a source field

    Using: 
        + Nedelec shape functions (instead of linear Lagrange)
        + The magnetostatic vector potential (instead of the scalar potential)
        + Magnetic Reluctance (instead of permeability)
        + Newton-Raphson iteration method

    Why:
        + Should be much more stable for non-linear media
            + Its actually much better for simulating iron structures
            - Its much worse at simulating Gd at low temperatures (280 K already breaks)
        + Can include Surface and Volume electric currents

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

function main(meshSize=0.0, localSize=0.0, showGmsh=false, verbose=true)

    # vacuum magnetic permeability
    mu0::Float64 = pi*4e-7

    # Temperature
    T::Float64 = 290.0

    # Applied field | Tesla
    Bext::Vector{Float64} = [0.5, 
                             0.0,
                             0.0]

    # Convergence criteria
    tol::Float64 = 1e-4 # Stop when B changes less than this (tesla)
    maxAtt::Int32 = 100 # Maximum number of iterations

    # Load Gd data automatically from a folder
    # data = DATA()
    # loadMaterial( data,
    #              "Materials", # Folder with materials
    #              "Gd_MFT",    # Data folder of target material
    #              "Gd",        # Material name
    #              7.9,         # density g/cm3
    #              T)

    # Load Fe data manually
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

    println("\nNumber of elements: ", mesh.nt)
    println("Number of Inside elements ", mesh.nInside)
    println("Number of edges: ", mesh.ne)
    println("Number of vertices: ", mesh.nv)
    println("Number of Inside nodes ", mesh.nInsideNodes)
    println("Number of surface elements ", mesh.ns)
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
    Ak::Matrix{Float64} = nedelecLocalStiffness(mesh)

    # Load vector
    RHS::Vector{Float64} = zeros(mesh.ne)
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

        for i in 1:6 # For each edge of the tetrahedron
            
            # Global node labels of current edge
            ndi, ndj = NodesFromLocalEdge(mesh, k, i)

            # Length of edge
            edgeLength = norm(mesh.p[1:3, nds[ndj]] - mesh.p[1:3, nds[ndi]])

            # 1st order Lagrange basis function (hat function)
            _, bi, ci, di = hat[:, ndi]
            _, bj, cj, dj = hat[:, ndj]

            # Curl element from Lagrange element
            curlN = 2.0*edgeLength*cross([bi, ci, di], [bj, cj, dj])

            # Global edge label
            edge = mesh.t[4+i, k]
            e = mesh.edge2localMap[edge] # Local edge label

            # Update the load vector on the local edge
            RHS[e] += mesh.VE[k]*dot(Bext, curlN)
        end

    end # Loop over the volume element labels

    # Relative magnetic reluctance
    nu::Vector{Float64} = ones(mesh.nt)
    dnu::Vector{Float64} = zeros(mesh.nt)

    # Prepare output
    u::Vector{Float64} = zeros(mesh.ne) # Magnetic vector potential on each edge

    Bfield = zeros(3, mesh.nt)            # 3D vector field B
    B::Vector{Float64} = zeros(mesh.nt)        # Norm of B
    Bold::Vector{Float64} = zeros(mesh.nt)     # Norm of B from previous iteration

    # Newton-Raphson iteration method
    att::Int32 = 0
    div::Float64 = Inf
    while div > tol && att < maxAtt

        att += 1
        Bold .= B

        # Global sparse stiffness matrix
        A = spzeros(mesh.ne, mesh.ne)
        n = 0
        for i in 1:6
            edge1 = mesh.edge2localMap[mesh.t[4+i, :]]
            
            for j in 1:6
                n += 1

                edge2 = mesh.edge2localMap[mesh.t[4+j, :]]
                A += sparse(edge1, edge2, Ak[n, :].*nu, mesh.ne, mesh.ne)
            end
        end # Global stiffness matrix

        # Tangential stiffness matrix
        if att > 1
            At = nedelecTangentialStiffness(mesh, dnu, Bfield, B)
        else
            At = spzeros(mesh.ne, mesh.ne)
        end

        # Correction to the magnetic scalar potential
        du = cg(A+At, RHS-A*u; verbose=false)

        # Update the vector potential
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

            for i in 1:6 # For each edge of the tetrahedron
                
                # Global node labels of the edge
                ndi, ndj = NodesFromLocalEdge(mesh, k, i)

                # Length of edge
                edgeLength = norm(mesh.p[1:3, nds[ndj]] - mesh.p[1:3, nds[ndi]])

                # 1st order Lagrange basis function (hat function)
                _, bi, ci, di = hat[:, ndi]
                _, bj, cj, dj = hat[:, ndj]

                # Curl element from Lagrange element
                curlN = 2.0*edgeLength*cross([bi, ci, di], [bj, cj, dj])

                # Global edge label
                edge = mesh.t[4+i, k]
                e = mesh.edge2localMap[edge] # Local edge label

                # Update vector field
                Bfield[:, k] += u[e].*curlN
            end
        end # Vector field B

        # Norm of B
        B .= 0.0
        for k in 1:mesh.nt
            B[k] = norm(Bfield[:, k])
        end

        # Update magnetic reluctance
        nu[mesh.InsideElements] .= spl(B[mesh.InsideElements])
        dnu[mesh.InsideElements] = spl_dnu(B[mesh.InsideElements])

        # Check for nans
        if any(x->!isfinite(x), nu) || any(x->!isfinite(x), dnu)
            error("Nans in the interpolation")
        end

        # Check deviation from previous result
        div = maximum(abs.(B .- Bold))
        verbose ? println(att, " , |B(n)-B(n-1)| = ", div, " , |y-Ax| = ", residue) : nothing

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
