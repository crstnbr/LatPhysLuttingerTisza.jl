################################################################################
#
#   FUNCTIONS FOR CONSTRUCTION OF K SPACE POINTS
#
#   FILE CONTAINS
#   - general building of points
#   - high symmetry mesh
#
#   ALL FUNCTIONS ARE NOT EXPORTED
#
################################################################################






################################################################################
#
#   FUNCTIONS FOR GENERAL CONSTRUCTION OF K SPACE POINTS
#
################################################################################

# contruct line between vertices with a certain resolution
function getPointsOnLine(
            point1 :: Vector{<:Real},
            point2 :: Vector{<:Real},
            resolution :: Integer
        ) :: Vector{Vector{Float64}}

    # the array of points which is later returned
    line = Vector{Vector{Float64}}(undef, resolution)
    # the multipliers between the two points
    alphas = range(0.0, stop=1.0, length=resolution)

    # set all line points
    for i in 1:resolution
        line[i] = (point1.*(1-alphas[i])) .+ (point2.*alphas[i])
    end

    # return the array of points
    return line
end

# get random points inside a box
function getPointsRandomInBox(
            center      :: Vector{<:Real},
            dimensions  :: Vector{<:Real},
            N           :: Integer
        )

    # the list of points to return later
    points = Vector{Vector{Float64}}(undef, N)
    # the dimension of points
    d = length(dimensions)

    # set all points
    for i in 1:N
        # set the point to a random point
        points[i] = rand(Float64,d)
        # shift the point into the box
        points[i] .*= 2
        points[i] .-= 1
        points[i] .*= dimensions
        points[i] .+= center
    end

    # return the array of points
    return points
end

# compute the center as half the diagonal from a set of vertices forming a BZ surface
function getPointCenter(
            points :: Vector{<:Vector{<:Real}}
        ) :: Vector{Float64}

    # return sum of all points and divide by number of points
    return sum(points) ./ length(points)
end
function getPointCenter(
            points :: Vector{<:Real} ...
        ) :: Vector{Float64}

    # return sum of all points and divide by number of points
    return sum(points) ./ length(points)
end





################################################################################
#
#   FUNCTIONS FOR HIGH SYMMETRY MESH
#
################################################################################

# contructor for symmetric mesh of given brillouin zone
function getPointsHighSymmetryMesh(
            bz              :: BZ,
            line_resolution :: Integer
            ;
            include_bz_corners :: Bool           = true,
            include_face_centers :: Bool         = true,
            include_bz_edges :: Bool             = true,  # edges of faces (i.e. lines connecting the corners)
            include_lines_in_faces :: Bool       = true,
            include_lines_to_face_center :: Bool = true,
            include_lines_to_gamma :: Bool       = true,
            include_random_face_points :: Bool   = true
        ) :: Vector{Vector{Float64}} where {D,B, P<:AbstractReciprocalPoint{D}, R<:AbstractReciprocalUnitcell{P,B}, BZ<:AbstractBrillouinZone{R}}

    # get the gamma point
    gamma = point(gammaPoint(reciprocalUnitcell(bz)))
    # the list of points that form the mesh
    mesh = Vector{Float64}[]

    # iterate over all faces of the BZ
    for f in faces(bz)

        # obtain the relevant corners
        corners = unique([corners(bz)[i] for i in f])
        # points that connect to gamma
        points_to_gamma = Vector{Float64}[]

        # INCLUDE CORNER POINTS IN MESH
        if include_bz_corners
            for c in corners
                # put the corner into the mesh
                push!(mesh, c)
                # put the corner into the list of points that connect to the gamma point
                push!(points_to_gamma, c)
            end
        end

        # INCLUDE FACE CENTERS IN MESH
        if include_face_centers
            # push the center into the mesh
            push!(mesh, getPointCenter(corners))
            # push the center into the list of connecting points
            push!(points_to_gamma, getPointCenter(corners))
        end

        # INCLUDE EDGES IN MESH
        if include_bz_edges
            # the edge that wraps in the point list
            line = getPointsOnLine(corners[end], corners[1], line_resolution)
            for p in line
                # put the point into the mesh
                push!(mesh, p)
            end
            # put the start, end and center of line into the gamma point connecting point list
            push!(points_to_gamma, corners[end])
            push!(points_to_gamma, corners[1])
            push!(points_to_gamma, getPointCenter(corners[1], corners[end]))
            # all other edges
            for i in 1:length(corners)-1
                # construct line
                line = getPointsOnLine(corners[i], corners[i+1], line_resolution)
                for p in line
                    # put the point into the mesh
                    push!(mesh, p)
                end
                # put the start, end and center of line into the gamma point connecting point list
                push!(points_to_gamma, corners[i])
                push!(points_to_gamma, corners[i+1])
                push!(points_to_gamma, getPointCenter(corners[i], corners[i+1]))
            end
        end

        # INCLUDE LINES IN THE FACE
        if include_lines_in_faces
            # iterate over all combination of points
            for i in 1:length(corners)
            for j in i+1:length(corners)
                # build the line from i to j
                line = getPointsOnLine(corners[i], corners[j], line_resolution)
                for p in line
                    push!(mesh, p)
                end
                # put the start, end and center of line into the gamma point connecting point list
                push!(points_to_gamma, corners[i])
                push!(points_to_gamma, corners[j])
                push!(points_to_gamma, getPointCenter(corners[i], corners[j]))
            end
            end
        end

        # INCLUDE LINES TO THE FACE CENTER
        if include_lines_to_face_center
            # build the face center
            center = getPointCenter(corners)
            # iterate over all points
            for i in 1:length(corners)
                # build the line from i to j
                line = getPointsOnLine(corners[i], center, line_resolution)
                for p in line
                    push!(mesh, p)
                end
            end
        end

        # INCLUDE LINES CONNECTING TO GAMMA POINT
        if include_lines_to_gamma
            # get a unique list of connecting points
            unique!(points_to_gamma)
            # build lines from these points to the gamma point (if they are not the gamma point)
            for ptg in points_to_gamma
                if dot(ptg.-gamma, ptg.-gamma) < 1e-8
                    continue
                end
                line = getPointsOnLine(ptg, gamma, line_resolution)
                for p in line
                    push!(mesh, p)
                end
            end
        end

    end

    # remove redundant points
    unique!(mesh)

    # return the mesh
    return mesh
end
