#include "Primitives.h"

#include <utility>
#include <iostream>


vec3f cross(vec3f v1, vec3f v2) {
    return vec3f(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

Light::Light(const vec3f &p, const float i) : position(p), intensity(i) {}

Material::Material(const vec4 &a, const vec3f &color, float spec, float r) : refractive_index(r), diff_spec_refl_refr(a),
                                                                            color(color),
                                                                            specular_exponent(spec) {}

Material::Material() : refractive_index(1), diff_spec_refl_refr(0.2, 0.8, 0, 0), color(0.1, 0, 0), specular_exponent(10) {}


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
        int x_coord = (long long) (point.x * dist / (EPS * 200)) % tex_w;
        int z_coord = (long long) (point.z * dist / (EPS * 200)) % tex_h;
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
    int i = (atan2f(p.z, p.x)/(2*M_PI)+1)*tex_w;
    i =  (i +100*tex_w)%tex_w;
    int j = acosf(p.y / r) / M_PI * tex_h;
    //std::cout <<p.y<< " " <<  j << std::endl;
    j =(j+100*tex_h)%tex_h;
    if (i < 0) {
        i = 0;
    }
    if (j < 0) {
        j = 0;
    }

    return texture[i+j*tex_w];
}


float sdBox( vec3f p, vec3f b ) {
    vec3f d = abs(p) - b;
    return fmin(fmax(d.x,fmax(d.y,d.z)),0.0) + vec3f(fmax(d.x,0.0), fmax(d.y,0.0), fmax(d.z,0.0)).norm();
}
Hit Cone::ray_intersect(const Ray &ray) const {
    return {};
}


float Cone::dist(const vec3f &p) const {
    vec3f point = p - pos;
    vec2f q = vec2f(c.x/c.y,-1.0)*h;
    vec2f w = vec2f(vec2f(point.x, point.z).norm(), point.y );
    vec2f a = w - q*clamp( dot(w,q)/dot(q,q), 0.0f, 1.0f );
    vec2f b = w - q*vec2f( clamp( w.x/q.x, 0.0f, 1.0f), 1.0 );
    float k = sign(q.y);
    float d = fminf(dot( a, a ),dot(b, b));
    float s = fmaxf( k*(w.x*q.y-w.y*q.x),k*(w.y-q.y)  );
    return sqrtf(d)*sign(s);
}

Material Cone::get_material(const vec3f &point) const {
    return material;
}

Cone::Cone(const vec3f &_pos, float _h, const vec2f &_c, const Material &mat) : pos(_pos), h(_h), c(_c), material(mat) {}

Hit Box::ray_intersect(const Ray &ray) const {
    return {};
}

Material Box::get_material(const vec3f &point) const {
    return material;
}

float Box::dist(const vec3f &point) const {
    auto p = point - pos;
    return sdBox(p, b);
    //vec3f d = abs(p) - b;
    //return fmin(fmax(d.x,fmax(d.y,d.z)),0.0) + vec3f(fmax(d.x,0.0), fmax(d.y,0.0), fmax(d.z,0.0)).norm();
}

Box::Box(const vec3f &_pos, const vec3f &_dims, const Material &_material) : pos(_pos), b(_dims), material(_material) {}

Hit Fractal::ray_intersect(const Ray &ray) const {
    return {};
}
//
float Fractal::dist(const vec3f &point) const {
//    return 0;
    vec3f p = point - pos;
    float d = sdBox(p, vec3f(3, 3, 3));
    float s = 1.0, da, db, dc, c;
    vec3f a, r;
    for (int i = 0; i < 4; i++) {
        a = mod( p*s, 2.0 )-1.0f;
        s *= 3;
        r = abs(1 - abs(a)*3);
        da = fmaxf(r.x,r.y);
        db = fmaxf(r.y,r.z);
        dc = fmaxf(r.z,r.x);
        c = (fminf(dc, fminf(da, db)) - 1)/s;
        d = fmaxf(d, c);
    }
    return d;
}

Material Fractal::get_material(const vec3f &point) const {
    return material;
}

Fractal::Fractal(const vec3f &_pos, const Material &mat) : pos(_pos), material(mat){}
