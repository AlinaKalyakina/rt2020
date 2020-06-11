#ifndef RT_MY_MATH_H
#define RT_MY_MATH_H
#include <cmath>
#include <algorithm>

template<typename T>
struct vec3
{
    T x, y, z;

    vec3<T>() : x(0), y(0), z(0) {}
    vec3<T>(T a, T b, T c) : x(a), y(b), z(c) {}
    vec3<T>(const vec3<T>& t) = default;

    inline float norm() const { return std::sqrt(x*x+y*y+z*z); }

    vec3<T>& normalize() {
        *this = (*this) * (1.0f / norm());
        return *this;
    }

    inline vec3<T> operator*(const T t) const {
        return {x*t, y*t, z*t};
    }

    inline vec3<T> operator-(const vec3<T>& t) const {
        return {x-t.x, y-t.y, z-t.z};
    }

    inline float operator*(const vec3<T>& t) const {
        return x*t.x + y * t.y + z * t.z;
    }

    inline vec3<T> operator-() const {
        return vec3<T>((*this)) * (-1.0);
    }

    inline vec3<T> operator+(const vec3<T>& t) const {
        return {x + t.x, y+t.y, z+t.z};
    }

    vec3<T>& operator+=(const vec3<T>& t) {
        x += t.x;
        y += t.y;
        z += t.z;
        return *this;
    }

    inline void crop(T min, T max) {
        x = std::max(min, std::min(x, max));
        y = std::max(min, std::min(y, max));
        z = std::max(min, std::min(z, max));
    }

    inline T operator[](int i) const {
        switch (i) {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
        }
    }

    inline T& operator[](int i) {
        switch (i) {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
        }
    }
};


typedef vec3<float> vec3f;
typedef vec3<int> vec3i;

struct vec4
{
    vec4() : x(0), y(0), z(0), w(0) {}
    vec4(float a, float b, float c, float d) : x(a), y(b), z(c), w(d) {}

    float x, y, z, w;
};


#endif //RT_MY_MATH_H
