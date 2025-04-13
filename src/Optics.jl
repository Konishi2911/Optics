module Optics

include("Geom.jl")
include("Elements.jl")
include("Rays.jl")

struct OpticalSystem
    lenses::Vector{Tuple{AbstractLens, Vector{Float64}}}
    walls::Vector{Wall}
end


function plot!(ax, os::OpticalSystem)
    # Plot the optical system
    # ax: Plots.Plot
    # os: OpticalSystem
    t = 0:0.01:1
    for (lens, pos) in os.lenses
        f = Optics.geom2d(lens, pos)
        f_x = [p[1] for p in f.(t)]
        f_y = [p[2] for p in f.(t)]
        lines!(ax, f_x, f_y, label="lens", color=:red)
    end
    for wall in os.walls
        w = Optics.geom2d(wall)
        w_x = [p[1] for p in w.(t)]
        w_y = [p[2] for p in w.(t)]
        lines!(ax, w_x, w_y, label="wall", color=:gray)
    end
end

function analyze_rays(source::PointSource, os::OpticalSystem, n_rays::Int = 100)
    # Calculate the refraction of a ray through a lens
    # r: PointSource
    # l: AbstractLens
    # lens_offset: offset of the lens center from the origin
    # Returns the refracted ray

    rays = []
    for i in n_rays
        # Generate a random ray from the source
        angle = rand() * 2 * π
        direction = [cos(angle), sin(angle)]
        ray = [SingleRay(source.position, direction, nothing)]
        push!(rays, ray)
    end

    # Analyze the rays through the optical system
    analyze_rays!(rays, os)

    return rays
end




function analyze_rays!(rays::Vector{SingleRay}, os::OpticalSystem)
    # Calculate the refraction of a ray through a lens
    # r: Vector{SingleRay}
    # l: AbstractLens
    # lens_offset: offset of the lens center from the origin
    # Returns the refracted ray

    for i in os.lenses |> eachindex
        (lens, lens_offset) = os.lenses[i]

        # Analyze the ray refraction
        Optics.calc_refract!(rays, lens, lens_offset)
    end

    # Check if the ray crosses the wall
    for k in os.walls |> eachindex
        wall = os.walls[k]

        Optics.calc_refract!(rays, wall)
    end
end


function calc_refract!(rays::Vector{SingleRay}, l::AbstractLens, lens_offset::Vector{Float64} = [0.0, 0.0])
    # Calculate the refraction of a ray through a lens
    # r: Vector{SingleRay}
    # l: AbstractLens
    # lens_offset: offset of the lens center from the origin
    # Returns the refracted ray

    while true
        offset_ray = SingleRay(rays[end].origin - lens_offset, rays[end].direction, rays[end].medium)

        # Check if the ray intersects with the lens
        intersected_points = []
        for (i, elem) in l.elements |> enumerate
            p = intersect_point(offset_ray, elem)
            if p !== nothing
                push!(intersected_points, (i, p))
            end
        end

        if isempty(intersected_points)
            # No intersection found, break the loop
            break
        end

        # Find the closest intersection point
        min_id = argmin([norm(p[2] - offset_ray.origin) for p in intersected_points])
        elem = l.elements[intersected_points[min_id][1]]
        p_i = intersected_points[min_id][2]  
        
        # local coordinate of intersection point on the element
        t = find_local(elem, p_i)

        # normal vector at the intersection point
        n = normal(elem, t)
        if dot(offset_ray.direction, n) < 0
            n = -n
        end

        if offset_ray.medium !== nothing 
            r_index_i = offset_ray.medium.refractive_index
            r_index_o = 1.0
        else
            r_index_i = 1.0
            r_index_o = l.refractive_index
        end

        # Calculate the angle of incidence
        dot_nd = dot(offset_ray.direction, n)
        cross_nd = offset_ray.direction[1] * n[2] - offset_ray.direction[2] * n[1]
        theta_i = atan(cross_nd, dot_nd)
        theta_r = asin(sin(theta_i) * r_index_i / r_index_o)
        d_theta = theta_r - theta_i

        # Calculate the new direction of the ray
        d = [
            offset_ray.direction[1] * cos(d_theta) + offset_ray.direction[2] * sin(d_theta),
            -offset_ray.direction[1] * sin(d_theta) + offset_ray.direction[2] * cos(d_theta)
        ]

        # Calculate the new origin of the ray
        p_o = [
            p_i[1] + lens_offset[1],
            p_i[2] + lens_offset[2]
        ] + d * 1e-6  # small offset to avoid numerical issues

        new_medium = offset_ray.medium === nothing ? l : nothing
        push!(rays, SingleRay(p_o, d, new_medium))
    end
end

function calc_refract!(rays::Vector{SingleRay}, wall::Wall)
    last_ray = rays[end]
    # Check if the ray intersects with the wall
    p = intersect_point(last_ray, Segment(wall.sp, wall.ep))
    if p !== nothing
        n = [0.0, 0.0]
        push!(rays, SingleRay(p, n, nothing))
        return rays
    else
        return rays
    end
end


end # module Optics