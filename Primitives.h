//
// Created by alina on 06.06.2020.
//

#ifndef RT_PRIMITIVES_H
#define RT_PRIMITIVES_H
#include "my_math.h"

#define EPS 1e-03
#define MAX_DIST 1000

struct Light {
    Light(const vec3 &p, float i);
    vec3 position;
    float intensity;
};

struct Material {
    float refractive_index;
    vec4 albedo;
    vec3 diffuse_color;
    float specular_exponent;
    Material(float r, const vec4 &a, const vec3 &color, float spec);
    Material();
};

struct Hit {
    bool hit = false;
    vec3 point;
    float dist{};
    vec3 n;
    Material material;

    Hit();
    Hit(vec3 hit_point, float d, vec3 norm, Material mat);
    operator bool ();
};

struct Primitive {
    Material material;
    explicit Primitive(const Material &m):material(m) {}
    virtual Hit ray_intersect(const vec3 &orig, const vec3 &dir) const = 0;
    virtual ~Primitive() = default;;
};

struct Sphere: Primitive {
    vec3 center;
    float radius;

    Sphere(const vec3 &c, float r, const Material &m);

    Hit ray_intersect(const vec3 &orig, const vec3 &dir) const override;

};

struct HorPlane: Primitive {
    float y = -4;
    explicit HorPlane(float d, Material mat);
    Hit ray_intersect(const vec3 &orig, const vec3 &dir) const override;
};

#endif //RT_PRIMITIVES_H
