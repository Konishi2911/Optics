using Plots

abstract type AbstractLens end

struct PlanoConvexLens <: AbstractLens
    diameter::Float64
    carvature_radius::Float64
    thickness::Float64
    refractive_index::Float64

    elements::Vector{AbstractGeometricElement}

    is_mirrored::Bool
end

struct PlanoConcaveLens <: AbstractLens
    diameter::Float64
    carvature_radius::Float64
    center_thickness::Float64
    refractive_index::Float64

    elements::Vector{AbstractGeometricElement}

    is_mirrored::Bool
end


function PlanoConvexLens(diameter::Float64, carvature_radius::Float64, thickness::Float64, refractive_index::Float64; is_mirrored::Bool = false)
    center = thickness - carvature_radius
    edge_thickness = thickness - carvature_radius + sqrt(carvature_radius^2 - (0.5 * diameter)^2)

    if is_mirrored
        theta_s = pi - asin(0.5 * diameter / carvature_radius)
        theta_e = pi + asin(0.5 * diameter / carvature_radius)
        elements = [
            Segment([0., -0.5 * diameter], [0, 0.5 * diameter]),
            Segment([0, 0.5 * diameter], [-edge_thickness, 0.5 * diameter]),
            Arc([-center, 0], carvature_radius, theta_s, theta_e),
            Segment([-edge_thickness, -0.5 * diameter], [0, -0.5 * diameter]),
        ]
    else
        theta_s = -asin(0.5 * diameter / carvature_radius)
        theta_e = asin(0.5 * diameter / carvature_radius)
        elements = [
            Segment([0., 0.5 * diameter], [0, -0.5 * diameter]),
            Segment([0, -0.5 * diameter], [edge_thickness, -0.5 * diameter]),
            Arc([center, 0], carvature_radius, theta_s, theta_e),
            Segment([edge_thickness, 0.5 * diameter], [0, 0.5 * diameter]),
        ]
    end

    return PlanoConvexLens(
        diameter,
        carvature_radius,
        thickness,
        refractive_index,
        elements,
        is_mirrored
    )
end

function PlanoConcaveLens(diameter::Float64, carvature_radius::Float64, center_thickness::Float64, refractive_index::Float64; is_mirrored::Bool = false)
    center = center_thickness + carvature_radius
    theta = asin(0.5 * diameter / carvature_radius)
    edge_thickness = carvature_radius * (1 - cos(theta)) + center_thickness

    if is_mirrored 
        theta_s = -theta
        theta_e = theta
        elements = [
            Segment([0., 0.5 * diameter], [0, -0.5 * diameter]),
            Segment([0, -0.5 * diameter], [-edge_thickness, -0.5 * diameter]),
            Arc([-center, 0], carvature_radius, theta_s, theta_e),
            Segment([-edge_thickness, 0.5 * diameter], [0, 0.5 * diameter]),
        ]
    else
        theta_s = pi - theta
        theta_e = pi + theta
        elements = [
            Segment([0., -0.5 * diameter], [0, 0.5 * diameter]),
            Segment([0, 0.5 * diameter], [edge_thickness, 0.5 * diameter]),
            Arc([center, 0], carvature_radius, theta_s, theta_e),
            Segment([edge_thickness, -0.5 * diameter], [0, -0.5 * diameter]),
        ]
    end

    return PlanoConcaveLens(
        diameter,
        carvature_radius,
        center_thickness,
        refractive_index,
        elements,
        is_mirrored
    )
end

function geom2d(lens::PlanoConvexLens, offset::Vector{Float64} = [0.0, 0.0])
    geom_fn = (t) -> begin
        if 0 <= t < 0.25
            s = t / 0.25
            return geom2d(lens.elements[1])(s) .+ offset
        elseif 0.25 <= t < 0.5
            s = (t - 0.25) / 0.25
            return geom2d(lens.elements[2])(s) .+ offset
        elseif 0.5 <= t < 0.75
            s = (t - 0.5) / 0.25
            return geom2d(lens.elements[3])(s) .+ offset
        elseif 0.75 <= t <= 1
            s = (t - 0.75) / 0.25
            return geom2d(lens.elements[4])(s) .+ offset
        else
            error("t must be in the range [0, 1]")
        end
    end

    return geom_fn
end

function geom2d(lens::PlanoConcaveLens, offset::Vector{Float64} = [0.0, 0.0])
    geom_fn = (t) -> begin
        if 0 <= t < 0.25
            s = t / 0.25
            return geom2d(lens.elements[1])(s) .+ offset
        elseif 0.25 <= t < 0.5
            s = (t - 0.25) / 0.25
            return geom2d(lens.elements[2])(s) .+ offset
        elseif 0.5 <= t < 0.75
            s = (t - 0.5) / 0.25
            return geom2d(lens.elements[3])(s) .+ offset
        elseif 0.75 <= t <= 1
            s = (t - 0.75) / 0.25
            return geom2d(lens.elements[4])(s) .+ offset
        else
            error("t must be in the range [0, 1]")
        end
    end

    return geom_fn
end