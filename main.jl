using LinearAlgebra, SparseArrays, IterativeSolvers

using Gmsh
include("gmsh_wrapper.jl")

function areaTriangle(xt,yt,zt)
    Atr = 0.5*sqrt(det([xt';yt';[1 1 1]])^2 + det([yt';zt';[1 1 1]])^2 + det([zt';xt';[1 1 1]])^2);
    return Atr
end # Area of the 3D triangle

function abcd(p,nodes,nd)
    n1,n2,n3 = nodes[nodes .!= nd]
    x = @view p[1, [nd, n1, n2, n3]]
    y = @view p[2, [nd, n1, n2, n3]]
    z = @view p[3, [nd, n1, n2, n3]]

    M = [1.0 x[1] y[1] z[1];
         1.0 x[2] y[2] z[2];
         1.0 x[3] y[3] z[3];
         1.0 x[4] y[4] z[4]]
    r = M\[1;0;0;0]

    return r[1],r[2],r[3],r[4]
end # Basis function coef.

function stiffnessMatrix(mesh,chi)
    A = zeros(mesh.nv,mesh.nv)
    for k in 1:mesh.nt
        nds = @view mesh.t[:,k]

        b = zeros(4,1); c = zeros(4,1); d = zeros(4,1);
        for i in 1:length(nds)
            _,b[i],c[i],d[i] = abcd(mesh.p,nds,nds[i])
        end
        # A[nds,nds] += chi[k].*mesh.VE[k].*(b*b' + c*c' + d*d')
    
        # Calculate the local matrix contribution
        local_mat = chi[k] * mesh.VE[k] * (b*b' + c*c' + d*d')
        
        # Add to the global sparse matrix
        for i in 1:length(nds)
            for j in 1:length(nds)
                A[nds[i], nds[j]] += local_mat[i,j]
            end
        end
    end

    return A
end

function lagrange(mesh)
    C = zeros(mesh.nv,1);
    for k in 1:mesh.nt
        nds = mesh.t[:,k];   # Nodes of that element

        # For each node of that element
        for i in 1:length(nds)
            C[nds[i]] = C[nds[i]] + mesh.VE[k]/4;
        end
    end
    return C
end # Integral of basis function

function BoundaryIntegral(mesh,F,normal,shell_id)
    RHS = zeros(mesh.nv,1);
    for s in 1:mesh.ne

        # Only integrate over the outer shell
        if !(mesh.surfaceT[4,s] in shell_id)
            continue
        end

        # Nodes of the surface triangle
        nds = @view mesh.surfaceT[1:3,s];
        
        # Area of surface triangle
        areaT = areaTriangle(mesh.p[1,nds],mesh.p[2,nds],mesh.p[3,nds]);
        
        RHS[nds] .+= dot(normal[:,s],F)*areaT/3;
    end

    return RHS
end # Integral of uniform field over surface


function main(meshSize=0,localSize=0,saveMesh=false)
    #=
        Makes a model with cubes and spheres and refines the mesh on the spheres
    
        Input:
            meshSize  - Mesh size (0 = let gmsh choose)
            localSize - Size of mesh in every volume beyond the container (0 for no local refinement)
            saveMesh  - Save mesh to a FEMCE compatible format 

    =#
    

    # Applied field
    mu0 = pi*4e-7      # vacuum magnetic permeability
    Hext = [1,0,0]     # T
    
    # Magnetic susceptibility
    suscetibility = 3

    # Create a geometry
    gmsh.initialize()

    # >> Model
    # Create an empty container
    box = addCuboid([0,0,0],[4,4,4])

    # List of cells inside the container
    cells = []

    # Add 1st sphere
    addSphere([0,0,0],0.5,cells)

    # Fragment to make a unified geometry
    _, fragments = gmsh.model.occ.fragment([(3, box)], cells)
    gmsh.model.occ.synchronize()

    # Update container volume ID
    box = fragments[1][1][2]

    # Generate Mesh
    mesh = Mesh(cells,meshSize,localSize,saveMesh)
    
    # Get bounding shell surface id
    shell_id = gmsh.model.getAdjacencies(3, box)[2]

    # Must remove the surface Id of the interior surfaces
    shell_id = shell_id[1:6] # All other, are interior surfaces

    # gmsh.fltk.run()
    gmsh.finalize()

    println("Number of elements ",mesh.nt)
    println("Number of nodes ",mesh.nv)
    println("Number of surface elements ",mesh.ne)
    println("Number of Inside nodes ",length(mesh.InsideNodes))

    # FEM
    chi = zeros(mesh.nt,1);
    chi[mesh.InsideElements] .= suscetibility

    # List of all surface triangle normals
    normal = zeros(3,mesh.ne);
    for i in 1:mesh.ne
        normal[:,i] = normal_surface(mesh.p,@view mesh.surfaceT[1:3,i]);
    end

    # Boundary conditions
    RHS = BoundaryIntegral(mesh,Hext,normal,shell_id)

    # Stiffness matrix
    A = @inbounds stiffnessMatrix(mesh,chi)

    # Lagrange multiplier technique
    Lag = lagrange(mesh)

    # Extend the matrix for the Lagrange multiplier technique
    mat = [A Lag;Lag' 0]

    # Magnetic scalar potential
    
# 1
    G = svd(mat)
    u = G\[-RHS;0]

# 2
    # u = pinv(mat)*[-RHS;0]


    println(norm(mat*u+[RHS;0]))
    println(u[end])

    u = u[1:mesh.nv]

end # end of main

main(10,0,true)
