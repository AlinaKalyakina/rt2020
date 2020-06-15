//
// Created by alina on 11.06.2020.
//

#include "Model.h"
#include <iostream>
#include <cassert>
#include <fstream>
#include <sstream>
#include <utility>




vec3f cross(const vec3f& v1, vec3f v2) {
    return vec3f(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

// fills verts and faces arrays, supposes .obj file to have "f " entries without slashes
Model::Model(const char *filename) : verts(), faces() {

    meshes = load_obj(filename);

    vec3f min, max;
    //get_bbox(min, max);
}

inline void str_split(const std::string &in,
                      std::vector<std::string> &out,
                      const std::string &token)
{
    out.clear();

    std::string temp;

    for (int i = 0; i < int(in.size()); i++)
    {
        std::string test = in.substr(i, token.size());

        if (test == token)
        {
            if (!temp.empty())
            {
                out.push_back(temp);
                temp.clear();
                i += (int)token.size() - 1;
            }
            else
            {
                out.push_back("");
            }
        }
        else if (i + token.size() >= in.size())
        {
            temp += in.substr(i, token.size());
            out.push_back(temp);
            break;
        }
        else
        {
            temp += in[i];
        }
    }
}


Mesh load_mesh(std::ifstream &in, std::map<std::string, Material> &materials, const std::vector<vec3f>& verts,
        const std::vector<vec3f>& normals, const std::vector<vec2f>& tex_coord, std::string& line){

    //std::string line;
    std::string tmp;
    std::vector<std::string> tmp_str_vec;
    std::vector<vec3i> faces;
    //vec3f v;
    std::string name;
    std::string trash_s;
    std::vector<Vertex> vertices;
    Material material;
    std::string token;
    while (!in.eof()) {
        //auto token = firstToken(line);
        if (line[0] == '#') {
            std::getline(in, line);
            continue;
        }
        std::istringstream iss(line.c_str());
        char trash;
        iss >> token;
        //str_split(line, tmp, " ");
        if (token == "usemtl") {
            iss >> tmp;
            material = materials[tmp];
        } else if (token == "f") {
            vec3i f;
            int cnt = 0;
            int idx_pos = 0, idx_n = -1, idx_tex = -1;
            //std::vector<Vertex> tmp_v;
            //std::
            while (iss >> tmp) { //foreach v1/vt1/vn1
                Vertex new_v;
                str_split(tmp, tmp_str_vec, "/");
                switch (tmp_str_vec.size()) {
                    case 1:
                        idx_pos = atoi(tmp_str_vec[0].c_str()) - 1;
                        break;
                    case 2:
                        idx_pos = atoi(tmp_str_vec[0].c_str()) - 1;
                        idx_tex = atoi(tmp_str_vec[1].c_str()) - 1;
                        break;
                    case 3:
                        idx_pos = atoi(tmp_str_vec[0].c_str()) - 1;
                        idx_tex = atoi(tmp_str_vec[1].c_str()) - 1;
                        idx_n = atoi(tmp_str_vec[2].c_str()) - 1;
                        break;
                }
                new_v.pos = verts[idx_pos];
                if (idx_n != -1) {
                    new_v.n = normals[idx_n];
                }
                if (idx_tex != -1) {
                    new_v.tex = tex_coord[idx_tex];
                }
                vertices.push_back(new_v);
                cnt++;
                if (cnt >= 3) { // have more than 3 vertives
                    faces.emplace_back(vec3i(vertices.size() - cnt, vertices.size() - 2, vertices.size() - 1));
                    if (idx_n == -1) { //no normals
                        vertices[vertices.size()-1].n = cross(vertices[vertices.size() - 2].pos-vertices[vertices.size() - cnt].pos,
                                        vertices[vertices.size() - 1].pos-vertices[vertices.size() - cnt].pos).normalize();

                    } //no normals
                }
            }
        } else if (token=="s") {

        } else { //faces ended
            return Mesh("unnamed", vertices, faces, material);
        }
        std::getline(in, line);

    }
    return Mesh("unnamed", vertices, faces, material);
}


std::vector<Mesh> load_obj(const char *filename) {
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) {
        std::cerr << "Failed to open " << filename << std::endl;
        exit(-1);
    }
    std::string line;
    std::map<std::string, Material> materials;
    //std::vector<std::string> tmp;
    std::vector<Mesh> meshes;
    std::vector<vec3f> verts;
    std::vector<vec3f> normals;
    std::vector<vec2f> tex_coord;
    vec3f v;
    std::string token;
    std::string tmp;
    std::getline(in, line);
    Mesh mesh;
    while (!in.eof()) {
        //std::getline(in, line);
        //auto token = firstToken(line);
        if (line[0] == '#') {
            std::getline(in, line);
            continue;
        }
        std::istringstream iss(line.c_str());
        //char trash;
        iss >> token;
        //str_split(line, tmp, " ");
        if (token == "v") {
            //iss >> trash;
            vec3f v;
            for (int i = 0; i < 3; i++) {
                iss >> v[i];
                v[i] /= 10;
            }
            std::swap(v.y, v.z);// = v.z, v.y;
            v.z -= 25;
            v.y -= 10;
            verts.emplace_back(v);
        } else if (token == "vn") {
            //iss >> trash_s;
            for (int i = 0; i < 3; i++) iss >> v[i];
            normals.emplace_back(v);
        } else if (token == "vt") {
            //iss >> trash >> trash;
            for (int i = 0; i < 2; i++) iss >> v[i];
            tex_coord.emplace_back(v.x, v.y);
        } else if (token == "g" || token == "f") {
            //iss >>trash;
            if (token == "g")  std::getline(in, line);
            mesh = load_mesh(in, materials, verts, normals, tex_coord, line);
            if (token == "g")
                iss >> mesh.name;
            else
                mesh.name = std::string("a");
            meshes.push_back(mesh);
            continue;
        } else if (token == "mtllib") {
            iss >> tmp;
            std::vector<std::string> temp;
            str_split(std::string(filename), temp, "/");

            std::string path_to_model;

            if (temp.size() != 1) {
                for (int i = 0; i < temp.size() - 1; i++) {
                    path_to_model += temp[i] + "/";
                }
            }

            //path_to_model += tmp;
            materials = load_materials((path_to_model + tmp).c_str(), path_to_model);

        } else if (token != ""){
            std::cerr << "unknown token " << token << std::endl;
        }
        token = "";
        std::getline(in, line);
    }
    //materials.clear();
    return std::move(meshes);
}


std::map<std::string, Material> load_materials(const char* path, const std::string &dir) {
    std::ifstream file(path);
    std::map<std::string, Material> materials;


    // If the file is not found return false
    if (!file.is_open()) {
        std::cerr << "Can't open .mtl";
        materials["default"] = Material();
        return std::move(materials);
    }

    //Material tempMaterial;

    bool listening = false;

    // Go through each line looking for material variables
    std::string curline;
    std::string token;
    std::string cur_name;
    Material tmp_material;
    while (std::getline(file, curline))
    {
        std::istringstream iss(curline.c_str());
        iss >> token;
        // new material and material name
        if (token == "newmtl")
        {
            if (listening) {
                // Generate the material

                // Push Back loaded Material
                materials[cur_name] = tmp_material;
            }
            iss >> cur_name;
            listening = true;
            tmp_material = Material();
        }
//        else if (token == "Ka")
//        {
//
//            std::vector<std::string> temp;
//            algorithm::split(algorithm::tail(curline), temp, " ");
//
//            if (temp.size() != 3)
//                continue;
//
//            tempMaterial.Ka.X = std::stof(temp[0]);
//            tempMaterial.Ka.Y = std::stof(temp[1]);
//            tempMaterial.Ka.Z = std::stof(temp[2]);
//        }
        // Diffuse Color
        else if (token == "Kd")
        {
            iss >> tmp_material.color.x;
            iss >> tmp_material.color.y;
            iss >> tmp_material.color.z;
            tmp_material.diff_spec_refl_refr[0] = 1;
        }
        // Specular Color
        else if (token == "Ks")
        {
            //    vec4 diff_spec_refl_refr;
            iss >> tmp_material.diff_spec_refl_refr[1];

        }
        // Specular Exponent
        else if (token == "Ns")
        {
            iss >> tmp_material.specular_exponent;
        }
        // Optical Density
        else if (token == "Ni")
        {
            iss >> tmp_material.refractive_index;
            //tempMaterial.Ni = std::stof(algorithm::tail(curline));
        }
        // Dissolve
//        if (algorithm::firstToken(curline) == "d")
//        {
//            tempMaterial.d = std::stof(algorithm::tail(curline));
//        }
//        // Illumination
//        if (algorithm::firstToken(curline) == "illum")
//        {
//            tempMaterial.illum = std::stoi(algorithm::tail(curline));
//        }
//        // Ambient Texture Map
//        if (algorithm::firstToken(curline) == "map_Ka")
//        {
//            tempMaterial.map_Ka = algorithm::tail(curline);
//        }
        // Diffuse Texture Map
//        if (algorithm::firstToken(curline) == "map_Kd")
//        {
//            tempMaterial.map_Kd = algorithm::tail(curline);
//        }
//        // Specular Texture Map
//        if (algorithm::firstToken(curline) == "map_Ks")
//        {
//            tempMaterial.map_Ks = algorithm::tail(curline);
//        }
//        // Specular Hightlight Map
//        if (algorithm::firstToken(curline) == "map_Ns")
//        {
//            tempMaterial.map_Ns = algorithm::tail(curline);
//        }
//        // Alpha Texture Map
//        if (algorithm::firstToken(curline) == "map_d")
//        {
//            tempMaterial.map_d = algorithm::tail(curline);
//        }
//        // Bump Map
//        if (algorithm::firstToken(curline) == "map_Bump" || algorithm::firstToken(curline) == "map_bump" || algorithm::firstToken(curline) == "bump")
//        {
//            tempMaterial.map_bump = algorithm::tail(curline);
//        }
    }
    if (listening) {
        materials[cur_name] = tmp_material;
    }
    return std::move(materials);

}


// Moller and Trumbore
bool Model::ray_triangle_intersect(const int &fi, const vec3f &orig, const vec3f &dir, float &tnear) const {
    return false;
//    vec3f edge1 = point(vert(fi,1)) - point(vert(fi,0));
//    vec3f edge2 = point(vert(fi,2)) - point(vert(fi,0));
//    vec3f pvec = cross(dir, edge2);
//    float det = dot(edge1,pvec);
//    if (det<1e-5) return false;
//
//    vec3f tvec = orig - point(vert(fi,0));
//    float u = dot(tvec,pvec);
//    if (u < 0 || u > det) return false;
//
//    vec3f qvec = cross(tvec, edge1);
//    float v = dot(dir,qvec);
//    if (v < 0 || u + v > det) return false;
//    return tnear>1e-5;
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
    float best_dist = std::numeric_limits<float>::max();
    //float cur_dist = std::numeric_limits<float>::max();
    //int best_t = 0;
    Hit cur_hit, best_hit;
    for (auto &mesh:meshes) {
        cur_hit = mesh.ray_intersect(ray);
        if (cur_hit.hit && cur_hit.dist < best_dist) {
            best_hit = cur_hit;
            best_dist = best_hit.dist;
        }
    }
    return best_hit;
//    for (int t=0; t<faces.size(); t++) {
//        float dist;
//        if (ray_triangle_intersect(t, ray.orig, ray.dir, cur_dist) && cur_dist < best_dist) {
//            best_dist = cur_dist;
//            hit.point = ray.orig + ray.dir * best_dist;
//            hit.dist = best_dist;
//            hit.hit = true;
//            vec3f v0 = point(vert(t, 0));
//            vec3f v1 = point(vert(t, 1));
//            vec3f v2 = point(vert(t, 2));
//            hit.n = cross(v1-v0, v2-v0).normalize();
//            hit.material =  material;
//            }
//    }

    //return hit;
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
Mesh::Mesh(std::string _name, std::vector<Vertex> _vertices, std::vector<vec3i> _faces, Material _material) :
name(std::move(_name)), vertices(std::move(_vertices)), faces(std::move(_faces)), material(_material){

}

Hit Mesh::ray_intersect(const Ray &ray) const {
    return Hit();
//    float best_dist = std::numeric_limits<float>::max();
//    vec3i best_face(-1, -1, -1);
//    float tnear;
//    for (const auto& f:faces){
//        auto v0 = vertices[f.x];
//        auto v1 = vertices[f.y];
//        auto v2 = vertices[f.z];
//        vec3f edge1 = v1.pos - v0.pos;
//        vec3f edge2 = v2.pos - v0.pos;
//        vec3f pvec = cross(ray.dir, edge2);
//        float det = dot(edge1,pvec);
//        if (det<1e-5) continue;
//
//        vec3f tvec = ray.orig - v0.pos;
//        float u = dot(tvec,pvec);
//        if (u < 0 || u > det) continue;
//
//        vec3f qvec = cross(tvec, edge1);
//        float v = dot(ray.dir,qvec);
//        if (v < 0 || u + v > det) continue;
//        tnear = dot(edge2,qvec) * (1./det);
//
//        if (tnear <= EPS*1e-2) {
//            return {};
//        }
//
//        if (best_dist > tnear) {
//            best_face = f;
//            best_dist = tnear;
//        }
//
//        //return tnear>1e-5;
//
//    }
//    if (best_face.x < 0){
//        return {};
//    }
//    Hit hit;
//    hit.hit = true;
//    hit.dist = best_dist;
//    hit.point = ray.orig + ray.dir*best_dist;
//    hit.material = material;
//    auto v0 = vertices[best_face.x];
//    auto v1 = vertices[best_face.y];
//    auto v2 = vertices[best_face.z];
//    auto a = (v0.pos-hit.point).norm();
//    auto b = (v1.pos-hit.point).norm();
//    auto c = (v2.pos-hit.point).norm();
//    hit.n = cross(v1.pos-v0.pos, v2.pos-v0.pos).normalize();
////    hit.n = (v0.n*(b*c) + v1.n*(a*c) + v2.n*(a*b))*(1/(b*c + a*c +a*b));
//    if (dot(hit.n, ray.dir) > 0) {
//        hit.n = -hit.n;
//    }
//    return std::move(hit);
}
