#=
    Wrapper to Gmsh model occ kernel, simplifies model creation
=#

using Gmsh

# Add 2d rectangle from its center
function addRectangle(position, dimensions, cells=[])

    x = position[1] - dimensions[1]/2
    y = position[2] - dimensions[2]/2

    id = gmsh.model.occ.addRectangle(x, y, 0.0, 
                                    dimensions[1], dimensions[2])

    push!(cells,(2, id))

    # Sync kernel before exiting
    gmsh.model.occ.synchronize()
    
    return id
end

# Add 2d disk (circular)
function addDisk(position, radius, cells=[])
    id = gmsh.model.occ.addDisk(position[1],
                                position[2],
                                0.0,
                                radius,
                                radius)
    push!(cells,(2, id))
    
    # Sync kernel before exiting
    gmsh.model.occ.synchronize()
    
    return id
end

# Add a cuboid based on its center
function addCuboid(position, dimensions, cells=[])
    #=
        Makes a cuboid based on its centroid position
        Updates the cells list in case this cuboid is not meant to be
        the container of the simulation
    =#

    r = position - dimensions/2
    box = gmsh.model.occ.addBox(r[1], r[2], r[3], dimensions[1], dimensions[2], dimensions[3])
    
    append!(cells, [(3, box)])

    # Sync kernel before exiting
    gmsh.model.occ.synchronize()

    return box
end # Make a cuboid based on its center

# Add a sphere
function addSphere(position, radius, cells=[])
    #=
        Inputs:
            Position vector
            radius value
            cells <- a list of volumes (cells) that are inside a container
    =# 
    
    # Add a sphere to the current model
    sphere = gmsh.model.occ.addSphere(position[1],position[2],position[3],radius)

    # If sphere is not the container
    append!(cells, [(3, sphere)])

    # Sync kernel before exiting
    gmsh.model.occ.synchronize()
    
    return sphere
end # Make a sphere

function addCylinder(position, axis, radius, cells=[])
    # 'position' defines the center of the bottom base of the cylinder
    # 'axis' defines the direction of the cylinder (and length)

    id = gmsh.model.occ.addCylinder( position[1], 
                                     position[2], 
                                     position[3], 
                                     axis[1], 
                                     axis[2], 
                                     axis[3],
                                     radius)

    append!(cells, [(3, id)])

    # Sync kernel before exiting
    gmsh.model.occ.synchronize()

    return id
end # Make a cylinder

# Create container based on current model surface (3D only)
function makeContainer(scale::Float64=5.0)

    # Get all surface entities
    surfaces = gmsh.model.getEntities(2)

    # Initialize min/max coordinates
    x_min = Inf
    y_min = Inf
    z_min = Inf
    x_max = -Inf
    y_max = -Inf
    z_max = -Inf

    # Find global bounding box of the STL
    for s in surfaces
        bb = gmsh.model.getBoundingBox(s[1], s[2])
        x_min = min(x_min, bb[1])
        y_min = min(y_min, bb[2])
        z_min = min(z_min, bb[3])
        x_max = max(x_max, bb[4])
        y_max = max(y_max, bb[5])
        z_max = max(z_max, bb[6])
    end

    Lx = x_max-x_min
    Ly = y_max-y_min
    Lz = z_max-z_min
    L = maximum([Lx,Ly,Lz])

    # Container position and dimensions
    center = [(x_min + x_max)/2, (y_min + y_max)/2, (z_min + z_max)/2]
    dimensions = scale*L

    box = addSphere(center,dimensions)

    return box
end # Create container based on current model surface

# Import cad geometry file
function importCAD(file::String, cells, setContainer::Bool=false, scale::Float64=5.0)
    #= 
        Import cad geometry file and create a container
        if there is none
    =#

    # Import CAD file (BREP, STEP or IGES)
    gmsh.model.occ.importShapes(file)
    volume = gmsh.model.occ.healShapes()
    _,volume = volume[1]

    gmsh.model.occ.synchronize()
    cells = append!(cells,[(3,volume)])

    # Make a container for the stl file
    if setContainer
        box = makeContainer(scale)
        return box
    end

end # Import cad geometry file

# Boolean subtract the geometries
function occCut(cell1, cell2)

    outDimTags = gmsh.model.occ.cut(cell1, cell2)[1]
    gmsh.model.occ.synchronize()

    return outDimTags
end

# Unify the volumes and get the surface IDs of the bounding shell
# function unifyModel(cells, box=-1)
#     # Works for both 2D and 3D geometries
    
#     dim = cells[1][1] # Get dimension of the model

#     inputBox::Bool = box > 0 # Check if a bounding cell was provided

#     # Fragment to make a unified geometry
#     if box > 0
#         fragments, _ = gmsh.model.occ.fragment(vcat(cells,[(dim, box)]), [])
#     else
#         fragments, _ = gmsh.model.occ.fragment(cells, [])
#     end

#     gmsh.model.occ.synchronize()

#     # Update cell ids
#     newCells = box > 0 ? fragments[1:end-1] : fragments[1:end]

#     if length(newCells) > length(cells)
#         println("Warning: 'unifyModel' outputs more cells it received \n")

#         # Update 'cells' with the new cells up to the same number of cells
#         for i in 1:length(cells)
#             cells[i] = (cells[i][1], newCells[i][2])
#         end

#         # Push the new cells
#         for i in length(cells)+1:length(newCells)
#             push!(cells, newCells[i])
#         end

#     else
#         cells .= newCells
#     end
    
#     # Set the box to the last volume
#     box = fragments[end][2]

#     # Get bounding shell surface id
#     shell_id = gmsh.model.getBoundary([(dim, box)], false, false, false) # (dim, tag)
#     shell_id = [s[2] for s in shell_id] # tag

#     # Volume surface ids
#     if inputBox
#         internal_surfaces = gmsh.model.getBoundary(cells, false, false, false) # (dim, tag)
    
#     else # The last cell is actually a box
#         internal_surfaces = gmsh.model.getBoundary(cells[1:end-1], false, false, false) # (dim, tag)
    
#     end
#     internal_surfaces = [s[2] for s in internal_surfaces] # tag


#     shell_id = setdiff(shell_id, internal_surfaces) # Only the outer surfaces

#     # Must separate return logic to not mess the logic of older code
#     # otherwise shell_id will be a tuple in some cases
#     if inputBox
#         return shell_id
#     else
#         return shell_id, box
#     end

# end # Unify the volumes

function unifyModel(cells, box = -1) # cells (dim, tag); box (tag of bounding shell)

    # The new cell IDs follow the order 
    # in which the cells were originally added

    dim = cells[1][1] # Get dimension of the model

    # Check if a bounding shell cell ID was provided
    inputBox::Bool = false
    if box > 0
        inputBox = true
    end

    # Unify the model by fragment()
    if !inputBox
        _, tagsMap = gmsh.model.occ.fragment(cells, [])
    else
        _, tagsMap = gmsh.model.occ.fragment(vcat(cells, [(dim, box)]), [])
    end

    gmsh.model.occ.synchronize()

    # Update the cells ID, preserving the order of when they were added to the model 
    for (i, tags) in enumerate(tagsMap)
        
        if tags[1][2] != box
            cells[i] = tags[1]
        end

        if length(tags) > 1 && box < 0 # Found the cell ID of the bounding shell
            box = tags[1][2]
        end
    end

    # Get bounding shell surface id
    shell_id = gmsh.model.getBoundary([(dim, box)], false, false, false) # (dim, tag)
    shell_id = [s[2] for s in shell_id] # Remove the dim, leave only the tag

    # Internal cells
    if !inputBox # No bounding shell cell ID provided
        
        innerCells = []

        # Find which cells are not the bounding shell
        for (i, cell) in enumerate(cells)
            if cell[2] != box
                push!(innerCells, cell)
            end
        end
    
    else # Every cell in 'cells' is an internal cell
        innerCells = cells

    end # Internal cells

    # Surface tags of cells inside bounding shell
    internal_surfaces = gmsh.model.getBoundary(innerCells, false, false, false) # (dim, tag)
    internal_surfaces = [s[2] for s in internal_surfaces] # tag

    # Only the outer surfaces
    shell_id = setdiff(shell_id, internal_surfaces) 

    return shell_id, box  
end


# Save variable to .txt
function save2file(fileName,input)
    # Saves matrix to a .txt file
    open(fileName, "w") do io
        for row in eachrow(input)
            println(io, join(row, " , "))  # Space-separated
        end
    end
end # Save matrix to .txt file

# Calculates the size of the variable and prints if verbose
function memory(var, verbose::Bool=true)
    mem::Float64 = Base.summarysize(var)/(1024*1024)
    
    if verbose
        println("Var size (MB): ", mem)
    end

    return mem
end