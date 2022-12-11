import FromFile: @from
using Debugger
using Grassmann, ProgressMeter, StaticArrays, LinearAlgebra
using Distributions
using IterTools
Point, Color = (SVector{3, Float64}, SVector{3, Float64})
red, white, blue, black =
    Color(1, 0, 0), Color(1, 1, 1), Color(0.5, 0.7, 1.0), Color(0, 0, 0)
PPM_FILE = "/Users/alokbeniwal/.julia/dev/RayTracingOneWeekend/target/raytraced.ppm"

normsq(x) = norm(x)^2

# hittable
abstract type Solid end

abstract type Material end

struct Sphere{R <: Real} <: Solid
    center::Point
    radius::R
    material::Material
end
struct Camera
    origin::Point
    horizontal::Point
    vertical::Point
    lower_left_corner::Point
end

function ray(u, v, cam::Camera)::Ray
    Ray(cam.origin, cam.lower_left_corner + u * cam.horizontal + v * cam.vertical)
end

struct HitRecord
    pt::Point
    normal::Point
    t::Float64
    material::Material
    front_face::Bool
end

struct Lambertian <: Material
    albedo::Color
end

struct Metal <: Material
    albedo::Color
end

# Ray = Tuple{Point, Point} # start, end
struct Ray
    origin::Point
    direction::Point # really a vector, affine space
end

function reflect(v, n)
    # v + 2(v ⋅ n)n
    # n + 2(v ⋅ n)v
    v - 2(v ⋅ n)n
end

function scatter(ray_in::Ray, hit_record::HitRecord, material::Lambertian)
    r = rand_unit_vector()

    scatter_direction = if hit_record.normal ≈ -r
        hit_record.normal
    else
        hit_record.normal + r
    end

    # TODO: return scattered and attenuation, the extra args are meant to be mutated. not on my watch.
    scattered = Ray(hit_record.pt, scatter_direction)
    return scattered
end

function scatter(ray_in::Ray, hit_record::HitRecord, material::Metal)
    # TODO: ray_in may be normalized(?)
    reflected = reflect(normalize(ray_in.direction), hit_record.normal)

    # TODO: return scattered and attenuation, the extra args are meant to be mutated. not on my watch.
    scattered = Ray(hit_record.pt, reflected)
    return scattered
end

at(ray::Ray, t::Real) = Ray(ray.origin, t * ray.direction)

# TODO: each ray can be done totally parallel to the others
# TODO: why doesn't t_min of 0.01 work?
function colorize(
    world::AbstractArray{<:Solid},
    ray::Ray;
    t_min = 0.5,#TODO: fix
    t_max = Inf,
    depth = 50,
)::Color
    if depth ≤ 0
        return black
    end

    hit_record = hit(world, ray; t_min = t_min, t_max = t_max)
    if hit_record !== nothing
        scattered = scatter(ray, hit_record, hit_record.material)
        # TODO: revert
        if (scattered.direction ⋅ hit_record.normal) > 0.0
            return hit_record.material.albedo .* colorize(
                world,
                scattered;
                t_min = t_min,
                t_max = t_max,
                depth = depth - 1,
            )
        else
            return black
        end
    end
    # TODO: understand this block
    t = 0.5 * (normalize(ray.direction).y + 1.0)
    (1.0 - t) * white + t * blue
end

function test_img(h, w)::Matrix{Color}
    mat = zeros(Color, (h, w))

    for row in h:-1:1, col in 1:w # Column major, so outer loop is *first*.
        r, g, b = (col - 1) / (w - 1), (row - 1) / (h - 1), 0.25
        mat[row, col] = Color(r, g, b)
    end
    mat
end

# return record or nothing
function hit(sphere::Sphere, ray::Ray; t_min = -Inf, t_max = Inf)::Union{HitRecord, Nothing}
    oc = ray.origin - sphere.center
    a = normsq(ray.direction)
    half_b = oc ⋅ ray.direction
    c = normsq(oc) - sphere.radius^2
    discriminant = half_b^2 - a * c
    if discriminant < 0
        return nothing
    end

    root = -(half_b + sqrt(discriminant)) / a
    if root < t_min || root > t_max
        # try other root
        root = -(half_b - sqrt(discriminant)) / a
        if root < t_min || root > t_max
            return nothing
        end
    end

    pt = at(ray, root).direction
    normal = (pt - sphere.center) / sphere.radius
    front_face = normal ⋅ ray.direction < 0.0

    HitRecord(pt, front_face ? normal : -normal, root, sphere.material, front_face)
end

function hit(
    objects::AbstractArray{<:Solid},
    ray::Ray;
    t_min = -Inf, # TODO: fix
    t_max = Inf,
)::Union{HitRecord, Nothing}
    closest, closest_t = nothing, t_max
    for object in objects
        hit_record = hit(object, ray; t_min = t_min, t_max = closest_t)
        if hit_record !== nothing && hit_record.t < closest_t
            closest_t, closest = hit_record.t, hit_record
        end
    end

    closest
end

function rand_unit_vector()::Point
    while true
        p = 2 * rand(Point) .- 1 # center to -1, 1
        if normsq(p) >= 1
            continue
        else
            return normalize(p)
        end
    end
end

function main()
    # Image
    aspect_ratio, img_width, viewport_height, focal_len = 16 / 9, 400, 2.0, 1.0
    # TODO: increase samples
    samples_per_pixel, max_depth = 100, 50
    img_height = Int(img_width ÷ aspect_ratio)
    viewport_width = viewport_height * aspect_ratio
    ground, center, left, right = Lambertian(Color(0.8, 0.8, 0.0)),
    Lambertian(Color(0.7, 0.3, 0.3)),
    # Lambertian(Color(0.8, 0.8, 0.8)),
    # Lambertian(Color(0.8, 0.6, 0.2))
    Metal(Color(0.8, 0.8, 0.8)),
    Metal(Color(0.8, 0.6, 0.2))
    # World
    world = [
        Sphere(Point(0.0, -100.5, -1.0), 100.0, ground),
        Sphere(Point(0.0, 0.0, -1.0), 0.5, center),
        Sphere(Point(-1.0, 0.0, -1.0), 0.5, left),
        Sphere(Point(1.0, 0.0, -1.0), 0.5, right),
    ]

    # Camera
    origin, horizontal, vertical =
        zero(Point), Point([viewport_width, 0.0, 0.0]), Point([0.0, viewport_height, 0.0])
    camera = Camera(
        origin,
        horizontal,
        vertical,
        origin - horizontal / 2 - vertical / 2 - Point(0.0, 0.0, focal_len),
    )

    function write_ppm(matrix::Matrix{Color}, filename = PPM_FILE)
        rows, cols = size(matrix)

        function write_color(io, c, samples_per_pixel)
            # TODO try writing to pixels without averaging only at end
            # the ^.5 = sqrt is gamma correction
            r, g, b = floor.(Int, 256 * clamp.((c / samples_per_pixel) .^ 0.5, 0.0, 0.9))

            write(io, "$r $g $b \n")
        end

        open(filename, "w") do io
            write(
                io,
                """
          P3
          $(cols) $(rows)
          255
          """,
            )
            for row in img_height:-1:1, col in 1:img_width
                # antialias
                pixel_color = red #TODO vary this, check for NaNs
                for (jitter_u, jitter_v) in ((rand(), rand()) for _ in 1:samples_per_pixel)
                    # TODO: check (u,v), since not 0 indexed
                    u, v = [
                        ((col - 1 + jitter_u) / (img_width - 1)),
                        ((row - 1 + jitter_v) / (img_height - 1)),
                    ]

                    r = ray(u, v, camera)
                    pixel_color += colorize(world, r, depth = max_depth)
                end
                write_color(io, pixel_color, samples_per_pixel)
            end
        end
    end

    write_ppm(zeros(Color, (img_height, img_width)))
end
main()

# TODO: intersect ray and viewport (line meet bivector)
# TODO: can broadcast to color in parallel
