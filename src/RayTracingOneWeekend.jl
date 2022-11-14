import FromFile: @from
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

struct Sphere{R <: Real} <: Solid
    center::Point
    radius::R
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
    front_face::Bool
end

# Ray = Tuple{Point, Point} # start, end
struct Ray
    origin::Point
    direction::Point # really a vector, affine space
end

at(ray::Ray, t::Real) = Ray(ray.origin, t * ray.direction)

# TODO: each ray can be done totally parallel to the others
function colorize(
    world::AbstractArray{<:Solid},
    ray::Ray;
    t_min = 0.01,
    t_max = Inf,
    depth = 50,
)::Color
    if depth ≤ 0
        return black
    end

    hit_record = hit(world, ray; t_min = t_min, t_max = t_max)
    if hit_record !== nothing
        tgt = hit_record.pt + hit_record.normal + rand_unit_vec()
        # reflect/absorb 1/2
        return colorize(world, Ray(hit_record.pt, tgt - hit_record.pt); depth = depth - 1) /
               2
    end
    t = (normalize(ray.direction).y + 1.0) / 2
    (1 - t) * white + t * blue
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
    front_face = normal ⋅ ray.direction < 0

    HitRecord(pt, front_face ? normal : -normal, root, front_face)
end

function hit(
    objects::AbstractArray{<:Solid},
    ray::Ray;
    t_min = -Inf,
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

function rand_unit_vec()::Point
    while true
        p = 2 * rand(Point) .- 1 # center to -1, 1
        if normsq(p) < 1
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

    # World
    world = [Sphere(Point(0, 0, -1), 0.5), Sphere(Point(0, -100.5, -1), 100)]

    # Camera
    origin, horizontal, vertical =
        zero(Point), Point([viewport_width, 0, 0]), Point([0, viewport_height, 0])
    camera = Camera(
        origin,
        horizontal,
        vertical,
        origin - horizontal / 2 - vertical / 2 - Point(0, 0, focal_len),
    )

    function write_ppm(matrix::Matrix{Color}, filename = PPM_FILE)
        rows, cols = size(matrix)

        function write_color(io, c, samples_per_pixel)
            # the sqrt is gamma correction
            r, g, b = Int.(floor.(256 * clamp.(sqrt.(c / samples_per_pixel), 0.0, 0.999)))

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
                pixel_color = black
                for (jitter_u, jitter_v) in ((rand(), rand()) for _ in 1:samples_per_pixel)
                    # TODO: check (u,v), since not 0 indexed
                    u, v = (col - 1 + jitter_u) / (img_width - 1),
                    (row - 1 + jitter_v) / (img_height - 1)

                    r = ray(u, v, camera)
                    pixel_color += colorize(world, r, depth = max_depth)
                end
                write_color(io, pixel_color, samples_per_pixel)
            end
        end
    end

    write_ppm(test_img(img_height, img_width))
end

main()

# TODO: intersect ray and viewport (line meet bivector)
# TODO: can broadcast to color in parallel
