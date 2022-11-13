import FromFile: @from
using Grassmann, ProgressMeter, StaticArrays, LinearAlgebra
using IterTools
Point, Color = (SVector{3, Float64}, SVector{3, Float64})
red = Color(1, 0, 0)
PPM_FILE = "/Users/alokbeniwal/.julia/dev/RayTracingOneWeekend/target/raytraced.ppm"

normsq(x) = norm(x)^2

# hittable
abstract type Solid end

struct Sphere{R <: Real} <: Solid
    center::Point
    radius::R
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

at(ray, t) = Ray(ray.origin, t * ray.direction)

# TODO: each ray can be done totally parallel to the others
function color(ray::Ray)
    # hardcode test for sphere
    hit_record = hit(Sphere(Point(0, 0, -1), 0.5), ray)
    if hit_record !== nothing && hit_record.t > 0
        N = normalize(at(ray, hit_record.t).direction - Point(0, 0, -1))
        return Color(N .+ 1) ./ 2
    end
    t = 0.5 * (ray.direction[2] + 1.0)
    t = (normalize(ray.direction)[2] + 1.0) / 2 # unit y direction
    # TODO: ones is white, other is blue
    white, blue = ones(Color), Color([0.5, 0.7, 1.0])
    (1 - t) * white + t * blue
end

function test_img(h, w)
    mat = zeros(Color, (h, w))

    for row in h:-1:1, col in 1:w # Column major, so outer loop is *first*.
        r, g, b = (col - 1) / (w - 1), (row - 1) / (h - 1), 0.25
        mat[row, col] = Color([r, g, b])
    end
    mat
end

# TODO: can broadcast to color in parallel
function color(rays) end

# return record or nothing
function hit(sphere::Sphere, ray::Ray, t_min = -Inf, t_max = Inf)::Union{HitRecord, Nothing}
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

    pt = at(ray, root).origin
    normal = (pt - sphere.center) / sphere.radius
    front_face = normal ⋅ ray.direction < 0

    HitRecord(pt, front_face ? normal : -normal, root, front_face)
end

function hit(
    xs::AbstractArray{<:Solid},
    ray::Ray,
    t_min = -Inf,
    t_max = Inf,
)::Union{HitRecord, Nothing}
    closest = nothing
    closest_t = t_max
    for x in xs
        hit_record = hit(x, ray, t_min, closest_t)
        if hit_record !== nothing && hit_record.t < closest_t
            closest = hit_record
            closest_t = hit_record.t
        end
    end
    closest
end

# TODO: intersect ray and viewport (line meet bivector)
function main()
    aspect_ratio = 16 / 9
    img_width = 400
    img_height = Int(img_width ÷ aspect_ratio)
    viewport_height = 2.0
    viewport_width = viewport_height * aspect_ratio
    focal_len = 1.0

    origin = zero(Point)

    horizontal = Point([viewport_width, 0, 0])
    vertical = Point([0, viewport_height, 0])
    inwards = Point(0, 0, focal_len)

    lower_left_corner = origin - horizontal / 2 - vertical / 2 - inwards

    function write_ppm(matrix::Matrix{Color}, filename = PPM_FILE)
        rows, cols = size(matrix)
        function write_color(io, c)
            r, g, b = Int.(floor.(255.999 * c))
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
                # TODO: check (u,v), since not 0 indexed
                u, v = (col - 1) / (img_width - 1), (row - 1) / (img_height - 1)
                ray = Ray(origin, lower_left_corner + u * horizontal + v * vertical)
                c = color(ray)
                write_color(io, c)
            end
        end
    end

    write_ppm(test_img(img_height, img_width))
end

main()
