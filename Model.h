//
// Created by alina on 11.06.2020.
//

#ifndef RT_MODEL_H
#define RT_MODEL_H

#include <vector>
#include <string>
//#include "my_math.h"
#include "Primitives.h"


struct Model : Primitive {
private:
    std::vector<vec3f> verts;
    std::vector<vec3i> faces;
    Material material = Material(vec4(0.9,  0, 0, 0), vec3f(.84, .21, .09), 1, 1.);
public:
    Model(const char *filename);

    int nverts() const;                          // number of vertices
    int nfaces() const;                          // number of triangles

    bool ray_triangle_intersect(const int &fi, const vec3f &orig, const vec3f &dir, float &tnear) const;
    Hit ray_intersect(const Ray& ray) const override;
    float dist(const vec3f& point) const override;
    Material get_material(const vec3f& point) const override;

    const vec3f &point(int i) const;                   // coordinates of the vertex i
    vec3f &point(int i);                   // coordinates of the vertex i
    int vert(int fi, int li) const;              // index of the vertex for the triangle fi and local index li
    void get_bbox(vec3f &min, vec3f &max); // bounding box for all the vertices, including isolated ones
};


//std::ostream& operator<<(std::ostream& out, Model &m);

#endif //RT_MODEL_H
