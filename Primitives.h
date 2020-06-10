//
// Created by alina on 06.06.2020.
//

#ifndef RT_PRIMITIVES_H
#define RT_PRIMITIVES_H

#include <vector>
#include "my_math.h"

#define EPS 1e-3
#define MAX_DIST 200

struct Ray {
    vec3 dir;
    vec3 orig;
    Ray(const vec3 &o, const vec3 &d);
};

struct Light {
    Light(const vec3 &p, float i);
    vec3 position;
    float intensity;
};

struct Material {
    vec4 diff_spec_refl_refr;
    vec3 color;
    float specular_exponent;
    float refractive_index;
    Material(const vec4 &a, const vec3 &color, float spec, float r);
    Material();
};

struct Hit {
    bool hit = false;
    vec3 point;
    float dist{};
    vec3 n;
    Material material;

    Hit();
    Hit(const vec3 &hit_point, float d, const vec3& norm, const Material &mat);
    operator bool ();
};

struct Primitive {
    //Material material;
    //explicit Primitive(const Material &m):material(m) {}
    virtual Hit ray_intersect(const Ray& ray) const = 0;
    virtual float dist(const vec3& point) const = 0;
    virtual Material get_material(const vec3& point) const = 0;
    virtual ~Primitive() = default;
};

struct Sphere: Primitive {
    vec3 center;
    float r;
    Material material;
    Sphere(const vec3 &c, float r, const Material &m);
    Hit ray_intersect(const Ray& ray) const override;
    float dist(const vec3& point) const override;
    Material get_material(const vec3& point) const override ;
};

struct HorPlane : Primitive {
    float y = -4;
    Material material;
    std::vector<vec3> texture;
    int tex_w = -1, tex_h=-1;
    HorPlane(float d, Material mat);
    HorPlane(float d, Material mat, std::vector<vec3> tex, int tex_width, int tex_height);
    Hit ray_intersect(const Ray &ray) const override;
    float dist(const vec3& point) const override;
    Material get_material(const vec3& point) const override;

};


struct Background {

    float r;
    std::vector<vec3> texture;
    int tex_w{}, tex_h{};
    const Sphere env;

    Background(float d, std::vector<vec3> tex, int tex_width, int tex_height);
    vec3 get_color(const Ray &ray) const;

};
#endif //RT_PRIMITIVES_H
