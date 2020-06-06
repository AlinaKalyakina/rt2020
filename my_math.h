//
// Created by alina on 06.06.2020.
//

#ifndef RT_MY_MATH_H
#define RT_MY_MATH_H
#include <cmath>
struct vec3
{
    float x, y, z;

    vec3() : x(0), y(0), z(0) {}
    vec3(float a, float b, float c) : x(a), y(b), z(c) {}
    vec3(const float* ptr) : x(ptr[0]), y(ptr[1]), z(ptr[0]) {}

    float norm() { return std::sqrt(x*x+y*y+z*z); }
    vec3& normalize() {
        *this = (*this) * (1.0 / norm());
        return *this;
    }
    vec3 operator*(const float t) const {
        auto res = vec3(*this);
        res.x *= t;
        res.y *= t;
        res.z *= t;
        return res;
    }
    vec3 operator-(const vec3& t) const {
        auto res = vec3(*this);
        res.x -= t.x;
        res.y -= t.y;
        res.z -= t.z;
        return res;
    }
    float operator*(const vec3& t) const {
        float res = x* t.x + y * t.y + z * t.z;
        return res;
    }
    vec3 operator-() const {
        auto res = vec3((*this));
        return res*(-1.0);
    }
    vec3 operator+(const vec3& t) const {
        auto res = vec3(*this);
        res.x += t.x;
        res.y += t.y;
        res.z += t.z;
        return res;
    }

    void crop(float t) {
        x = std::min(x, t);
        y = std::min(y, t);
        z = std::min(z, t);

    }
};

struct vec4
{
    vec4() : x(0), y(0), z(0), w(0) {}
    vec4(float a, float b, float c, float d) : x(a), y(b), z(c), w(d) {}
    vec4(vec3 a, float b) {
        x = a.x;
        y = a.y;
        z = a.z;
        w = b;
    }
//    vec3 xyz(){
//
//    }
    float x, y, z, w;
};

struct mat4 {
    mat4() { identity(); }

    mat4(const float arr[16]) {
        row[0] = vec4(arr[0], arr[1], arr[2], arr[3]);
        row[1] = vec4(arr[4], arr[5], arr[6], arr[7]);
        row[2] = vec4(arr[8], arr[9], arr[10], arr[11]);
        row[3] = vec4(arr[12], arr[13], arr[14], arr[15]);
    }

    void identity() {
        row[0] = vec4(1, 0, 0, 0);
        row[1] = vec4(0, 1, 0, 0);
        row[2] = vec4(0, 0, 1, 0);
        row[3] = vec4(0, 0, 0, 1);
    }

    static inline vec4 mul(mat4 m, vec4 v)
    {
        vec4 res;
        res.x = m.row[0].x*v.x + m.row[0].y*v.y + m.row[0].z*v.z + m.row[0].w*v.w;
        res.y = m.row[1].x*v.x + m.row[1].y*v.y + m.row[1].z*v.z + m.row[1].w*v.w;
        res.z = m.row[2].x*v.x + m.row[2].y*v.y + m.row[2].z*v.z + m.row[2].w*v.w;
        res.w = m.row[3].x*v.x + m.row[3].y*v.y + m.row[3].z*v.z + m.row[3].w*v.w;
        return res;
    }

    vec4 row[4];
};

#endif //RT_MY_MATH_H
