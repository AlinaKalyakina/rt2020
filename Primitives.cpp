#include "Primitives.h"

#include <utility>


Light::Light(const vec3 &p, const float i) : position(p), intensity(i) {}

Material::Material(const vec4 &a, const vec3 &color, float spec, float r) : refractive_index(r), diff_spec_refl_refr(a),
                                                                            color(color),
                                                                            specular_exponent(spec) {}

Material::Material() : refractive_index(1), diff_spec_refl_refr(1, 0, 0, 0), color(), specular_exponent() {}


Hit::Hit() = default;

Hit::Hit(const vec3 &hit_point, float d, const vec3 &norm, const Material &mat) : hit(true), point(hit_point), dist(d),
                                                                                  n(norm),
                                                                                  material(mat) {}

Hit::operator bool() { return hit; }

Sphere::Sphere(const vec3 &c, const float r, const Material &m) : material(m), center(c), r(r) {}

Hit Sphere::ray_intersect(const Ray &ray) const {

    vec3 L = center - ray.orig;
    float tca = L * ray.dir;
    float d2 = L * L - tca * tca;
    if (d2 > r * r)
        return {};
    float thc = std::sqrt(r * r - d2);
    float t0 = tca - thc;
    float t1 = tca + thc;
    if (t0 < 0) t0 = t1;
    if (t0 < 0) {
        return {};
    } else {
        vec3 hit = ray.orig + ray.dir * t0;
        return {hit, t0, (hit - center).normalize(), material};
    }
}


HorPlane::HorPlane(float d, Material m) : material(m), y(d) {}

Hit HorPlane::ray_intersect(const Ray &ray) const {
    auto d = y - ray.orig.y;
    auto t = d / ray.dir.y;
    if (t > 0 && t < MAX_DIST) { //watch in one side and n
        auto hit_point = ray.orig + ray.dir * t;
        if (tex_w == -1) {
            return {hit_point, t, vec3(0, -d / std::fabs(d), 0), material};
        } else {
            auto dist = hit_point.norm();
            int x_coord = (long long)(hit_point.x *dist/ (EPS*200)) % tex_w;
            int z_coord = (long long)(hit_point.z *dist/ (EPS*200)) % tex_h;
            if (x_coord < 0) {
                x_coord += tex_w;
            }
            if (z_coord < 0) {
                z_coord += tex_h;
            }
            auto ret_mat = Material(material);
            ret_mat.color = texture[x_coord + tex_w * z_coord];
            return {hit_point, t, vec3(0, -d / std::fabs(d), 0), ret_mat};

        }
    } else {
        return {};
    }
}

HorPlane::HorPlane(float d, Material mat, std::vector<vec3> tex, int tex_width, int tex_height) :
        y(d), material(mat), texture(std::move(tex)), tex_h(tex_height), tex_w(tex_width) {}

Ray::Ray(const vec3 &o, const vec3 &d) : dir(d), orig(o) {}


Background::Background(float t, const std::vector<vec3> &tex,int tex_width, int tex_height):
        r(t), texture(tex), tex_h(tex_height), tex_w(tex_width), env(Sphere(vec3(0.,0.,0.), t, Material())){

}

vec3 Background::get_color(const Ray &ray) const {
    //env.ray_intersect(ray);
    vec3 p = env.ray_intersect(ray).point;
    int i = (atan2(p.z, p.x)/(2*M_PI)+1)*tex_w;
    i %= tex_w;
    int j = acos(p.y / r) / M_PI * tex_h;
    return texture[i+j*tex_w];
}
