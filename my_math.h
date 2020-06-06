//
// Created by alina on 06.06.2020.
//

#ifndef RT_MY_MATH_H
#define RT_MY_MATH_H
#include <cmath>
#include <algorithm>
struct vec3
{
    float x, y, z;

    vec3() : x(0), y(0), z(0) {}
    vec3(float a, float b, float c) : x(a), y(b), z(c) {}
    explicit vec3(const float* ptr) : x(ptr[0]), y(ptr[1]), z(ptr[0]) {}

    float norm() { return std::sqrt(x*x+y*y+z*z); }
    vec3& normalize() {
        *this = (*this) * (1.0f / norm());
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

    void crop(float min, float max) {
        x = std::max(min, std::min(x, max));
        y = std::max(min, std::min(y, max));
        z = std::max(min, std::min(z, max));

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


#endif //RT_MY_MATH_H
