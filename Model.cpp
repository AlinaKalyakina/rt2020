//
// Created by alina on 11.06.2020.
//

#include "Model.h"
#include <iostream>
#include <cassert>
#include <fstream>
#include <sstream>

vec3f cross(vec3f v1, vec3f v2) {
    return vec3f(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

// fills verts and faces arrays, supposes .obj file to have "f " entries without slashes
Model::Model(const char *filename) : verts(), faces() {
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) {
        std::cerr << "Failed to open " << filename << std::endl;
        return;
    }
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            vec3f v;
            for (int i=0;i<3;i++) iss >> v[i];
            v.z -= 5;
            verts.push_back(v);
        } else if (!line.compare(0, 2, "f ")) {
            vec3i f;
            int idx, cnt=0;
            iss >> trash;
            while (iss >> idx) {
                idx--; // in wavefront obj all indices start at 1, not zero
                f[cnt++] = idx;
            }
            if (3==cnt) faces.push_back(f);
        }
    }
    std::cerr << "# v# " << verts.size() << " f# "  << faces.size() << std::endl;

    vec3f min, max;
    get_bbox(min, max);
}

// Moller and Trumbore
bool Model::ray_triangle_intersect(const int &fi, const vec3f &orig, const vec3f &dir, float &tnear) const {
    vec3f edge1 = point(vert(fi,1)) - point(vert(fi,0));
    vec3f edge2 = point(vert(fi,2)) - point(vert(fi,0));
    vec3f pvec = cross(dir, edge2);
    float det = edge1*pvec;
    if (det<1e-5) return false;

    vec3f tvec = orig - point(vert(fi,0));
    float u = tvec*pvec;
    if (u < 0 || u > det) return false;

    vec3f qvec = cross(tvec, edge1);
    float v = dir*qvec;
    if (v < 0 || u + v > det) return false;

    tnear = edge2*qvec * (1./det);
    return tnear>1e-5;
}


int Model::nverts() const {
    return (int)verts.size();
}

int Model::nfaces() const {
    return (int)faces.size();
}

void Model::get_bbox(vec3f &min, vec3f &max) {
    min = max = verts[0];
    for (int i=1; i<(int)verts.size(); ++i) {
        for (int j=0; j<3; j++) {
            min[j] = std::min(min[j], verts[i][j]);
            max[j] = std::max(max[j], verts[i][j]);
        }
    }
//    std::cerr << "bbox: [" << min << " : " << max << "]" << std::endl;
}

const vec3f &Model::point(int i) const {
    assert(i>=0 && i<nverts());
    return verts[i];
}


vec3f &Model::point(int i) {
    assert(i>=0 && i<nverts());
    return verts[i];
}


int Model::vert(int fi, int li) const {
    assert(fi>=0 && fi<nfaces() && li>=0 && li<3);
    return faces[fi][li];
}


Hit Model::ray_intersect(const Ray &ray) const {
    float duck_dist = std::numeric_limits<float>::max();
    float cur_dist = std::numeric_limits<float>::max();
    int best_t = 0;
    Hit hit;
    for (int t=0; t<faces.size(); t++) {
        float dist;
        if (ray_triangle_intersect(t, ray.orig, ray.dir, cur_dist) && cur_dist<duck_dist) {
            duck_dist = cur_dist;
            hit.point = ray.orig + ray.dir*duck_dist;
            hit.dist = duck_dist;
            hit.hit = true;
            vec3f v0 = point(vert(t, 0));
            vec3f v1 = point(vert(t, 1));
            vec3f v2 = point(vert(t, 2));
            hit.n = cross(v1-v0, v2-v0).normalize();
            hit.material =  material;
            }
    }

    return hit;
}


float Model::dist(const vec3f &point) const {
    return 0;
}

Material Model::get_material(const vec3f &point) const {
    return material;
}

//std::ostream& operator<<(std::ostream& out, Model &m) {
//    for (int i=0; i<m.nverts(); i++) {
//        out << "v " << m.point(i) << std::endl;
//    }
//    for (int i=0; i<m.nfaces(); i++) {
//        out << "f ";
//        for (int k=0; k<3; k++) {
//            out << (m.vert(i,k)+1) << " ";
//        }
//        out << std::endl;
//    }
//    return out;
//}
