#ifndef RT_MY_MATH_H
#define RT_MY_MATH_H
#include <cmath>
#include <algorithm>
struct vec3
{
    float x, y, z;

    vec3() : x(0), y(0), z(0) {}
    vec3(float a, float b, float c) : x(a), y(b), z(c) {}
    vec3(const vec3& t) = default;

    inline float norm() const { return std::sqrt(x*x+y*y+z*z); }

    vec3& normalize() {
        *this = (*this) * (1.0f / norm());
        return *this;
    }

    inline vec3 operator*(const float t) const {
        return {x*t, y*t, z*t};
    }

    inline vec3 operator-(const vec3& t) const {
        return {x-t.x, y-t.y, z-t.z};
    }

    inline float operator*(const vec3& t) const {
        return x*t.x + y * t.y + z * t.z;
    }

    inline vec3 operator-() const {
        return vec3((*this)) * (-1.0);
    }

    inline vec3 operator+(const vec3& t) const {
        return {x + t.x, y+t.y, z+t.z};
    }

    vec3& operator+=(const vec3& t) {
        x += t.x;
        y += t.y;
        z += t.z;
        return *this;
    }

    inline void crop(float min, float max) {
        x = std::max(min, std::min(x, max));
        y = std::max(min, std::min(y, max));
        z = std::max(min, std::min(z, max));
    }
};

struct vec4
{
    vec4() : x(0), y(0), z(0), w(0) {}
    vec4(float a, float b, float c, float d) : x(a), y(b), z(c), w(d) {}

    float x, y, z, w;
};


#endif //RT_MY_MATH_H
