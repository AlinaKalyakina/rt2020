#include "Primitives.h"

#include <utility>

vec3f cross(vec3f v1, vec3f v2) {
    return vec3f(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

Light::Light(const vec3f &p, const float i) : position(p), intensity(i) {}

Material::Material(const vec4 &a, const vec3f &color, float spec, float r) : refractive_index(r), diff_spec_refl_refr(a),
                                                                            color(color),
                                                                            specular_exponent(spec) {}

Material::Material() : refractive_index(1), diff_spec_refl_refr(0.2, 0.8, 0, 0), color(0.5, 0, 0), specular_exponent(10) {}


Hit::Hit() = default;

Hit::Hit(const vec3f &hit_point, float d, const vec3f &norm, const Material &mat) : hit(true), point(hit_point), dist(d),
                                                                                  n(norm),
                                                                                  material(mat) {}

Hit::operator bool() { return hit; }

Sphere::Sphere(const vec3f &c, const float radius, const Material &m) : material(m), center(c), r(radius) {}

Hit Sphere::ray_intersect(const Ray &ray) const {
    vec3f L = center - ray.orig;
    float tca = dot(L,ray.dir);
    float d2 = dot(L,L) - tca * tca;
    if (d2 > r * r)
        return {};
    float thc = std::sqrt(r * r - d2);
    float t0 = tca - thc;
    float t1 = tca + thc;
    if (t0 < 0) t0 = t1;
    if (t0 < 0) {
        return {};
    } else {
        vec3f hit = ray.orig + ray.dir * t0;
        return {hit, t0, (hit - center).normalize(), material};
    }
}

float Sphere::dist(const vec3f &point) const {
    return (center-point).norm()-r;
}

Material Sphere::get_material(const vec3f &point) const {
    return material;
}


HorPlane::HorPlane(float d, Material m) : material(m), y(d) {}

Hit HorPlane::ray_intersect(const Ray &ray) const {
    auto d = y - ray.orig.y;
    auto t = d / ray.dir.y;
    if (t > 0 && t < MAX_DIST) { //watch in one side and n
        auto hit_point = ray.orig + ray.dir * t;
        return {hit_point, t, vec3f(0, -d / std::fabs(d), 0), get_material(hit_point)};
    } else {
        return {};
    }
}

HorPlane::HorPlane(float d, Material mat, std::vector<vec3f> tex, int tex_width, int tex_height) :
        y(d), material(mat), texture(std::move(tex)), tex_h(tex_height), tex_w(tex_width) {}

float HorPlane::dist(const vec3f &point) const {
    return point.y - y;
}

Material HorPlane::get_material(const vec3f &point) const {
    if (tex_w == -1) {
        return material;
    } else {
        auto dist = point.norm();
        int x_coord = (long long) (point.x * dist / (1e-3 * 200)) % tex_w;
        int z_coord = (long long) (point.z * dist / (1e-3 * 200)) % tex_h;
        if (x_coord < 0) {
            x_coord += tex_w;
        }
        if (z_coord < 0) {
            z_coord += tex_h;
        }
        auto ret_mat = Material(material);
        ret_mat.color = texture[x_coord + tex_w * z_coord];
        return ret_mat;

    }
}
Ray::Ray(const vec3f &o, const vec3f &d) : dir(d), orig(o) {}


Background::Background(float t, std::vector<vec3f> tex, int tex_width, int tex_height):
        r(t), texture(std::move(tex)), tex_h(tex_height), tex_w(tex_width), env(Sphere(vec3f(0., 0., 0.), t, Material())){

}

vec3f Background::get_color(const Ray &ray) const {
    //env.ray_intersect(ray);
    vec3f p = env.ray_intersect(ray).point;
    int i = (atan2(p.z, p.x)/(2*M_PI)+1)*tex_w;
    i =  (i +2*tex_w)%tex_w;
    int j = acos(p.y / r) / M_PI * tex_h;
    j = (j+tex_h)%tex_h;
    return texture[i+j*tex_w];
}

Triangle::Triangle(const vec3f &a, const vec3f &b, const vec3f &c, const Material& m): v0(a), v1(b), v2(c),
    material(m), n(cross(v1-v0, v2-v0).normalize()) {

}

Hit Triangle::ray_intersect(const Ray &ray) const {
    vec3f edge1 = v1-v0;
    vec3f edge2 = v2-v0;
    vec3f pvec = cross(ray.dir, edge2);
    float det = dot(edge1,pvec);
    if (det<EPS/100) return {};

    vec3f tvec = ray.orig - v0;
    float u = dot(tvec,pvec);
    if (u < 0 || u > det) return {};

    vec3f qvec = cross(tvec, edge1);
    float v = dot(ray.dir,qvec);
    if (v < 0 || u + v > det) return {};

    float tnear = dot(edge2,qvec) * (1./det);
    if (tnear <- EPS/100) return {};

    Hit hit;
    hit.hit = true;
    hit.point = ray.orig + ray.dir*tnear;
    hit.material = material;
    hit.dist = tnear;
    hit.n = dot(n,ray.dir) < 0?n:-n;
    return hit;
}

float Triangle::dist(const vec3f &point) const {
    return 0;
//    TODO
}

Material Triangle::get_material(const vec3f &point) const {
    //TODO
    return Material();
}

float sdBox( vec3f p, vec3f b ) {
    vec3f d = abs(p) - b;
    return fmin(fmax(d.x,fmax(d.y,d.z)),0.0) + vec3f(fmax(d.x,0.0), fmax(d.y,0.0), fmax(d.z,0.0)).norm();
}
//
//float Fractal::dist(const vec3f &p) const {
//    float d = sdBox(p, vec3f(1, 1, 1));
//    float s = 1.0, da, db, dc, c;
//    vec3f a, r;
//    for (int i = 0; i < 4; i++) {
//        a = mod( p*s, 2.0 )-1.0;
//        s *= 3;
//        r = abs(vec3f(1, 1,1)- abs(a)*3);
//        da = fmax(r.x,r.y);
//        db = fmax(r.y,r.z);
//        dc = fmax(r.z,r.x);
//        c = (fmin(dc, fmin(da, db)) - 1)/s;
//        d = fmax(d, c);
//    }
//    return d;
//}
Hit Cone::ray_intersect(const Ray &ray) const {
    return Hit();
}


float Cone::dist(const vec3f &p) const {
//    float sdCone( in vec3 p, in vec2 c, float h )
//    {
//        // c is the sin/cos of the angle, h is height
//        // Alternatively pass q instead of (c,h),
//        // which is the point at the base in 2D
//        vec2 q = h*vec2(c.x/c.y,-1.0);
//
//        vec2 w = vec2( length(p.xz), p.y );
//        vec2 a = w - q*clamp( dot(w,q)/dot(q,q), 0.0, 1.0 );
//        vec2 b = w - q*vec2( clamp( w.x/q.x, 0.0, 1.0 ), 1.0 );
//        float k = sign( q.y );
//        float d = min(dot( a, a ),dot(b, b));
//        float s = max( k*(w.x*q.y-w.y*q.x),k*(w.y-q.y)  );
//        return sqrt(d)*sign(s);
//    }


    vec3f point = p - pos;
    vec2f q = vec2f(c.x/c.y,-1.0)*h;
    vec2f w = vec2f(vec2f(point.x, point.z).norm(), point.y );
    vec2f a = w - q*clamp( dot(w,q)/dot(q,q), 0.0f, 1.0f );
    vec2f b = w - q*vec2f( clamp( w.x/q.x, 0.0f, 1.0f), 1.0 );
    float k = sign(q.y);
    float d = fmin(dot( a, a ),dot(b, b));
    float s = fmax( k*(w.x*q.y-w.y*q.x),k*(w.y-q.y)  );
    return sqrt(d)*sign(s);
   // return 0;
}

Material Cone::get_material(const vec3f &point) const {
    return material;
}

Cone::Cone(const vec3f &_pos, float _h, const vec2f &_c, const Material &mat) : pos(_pos), h(_h), c(_c), material(mat) {}
