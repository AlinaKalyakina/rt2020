#include <cmath>
#include <algorithm>

#ifndef RT_MY_MATH_H
#define RT_MY_MATH_H

template<typename T>
struct vec3
{
    T x, y, z;

    vec3<T>() : x(0), y(0), z(0) {}
    vec3<T>(T a, T b, T c) : x(a), y(b), z(c) {}
    vec3<T>(const vec3<T>& t) = default;
    explicit vec3<T>(T t):x(t), y(t), z(t) {}

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

    inline vec3<T> operator-(float t) const {
        return {x-t, y-t, z-t};
    }

    inline vec3<T> operator*(const vec3<T>& t) const {
        return {x*t.x, y * t.y, z * t.z};
    }

    inline vec3<T> operator-() const {
        return vec3<T>((*this)) * (-1.0);
    }

    inline vec3<T> operator+(const vec3<T>& t) const {
        return {x + t.x, y+t.y, z+t.z};
    }

    inline vec3<T> operator/(const vec3<T>& t) const {
        return {x / t.x, y/t.y, z/t.z};
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
    inline float operator[](int i) const {
        switch (i) {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
            case 3:
                return w;
        }
    }

    inline float& operator[](int i) {
        switch (i) {
            case 0:
                return x;
            case 1:
                return y;
            case 2:
                return z;
            case 3:
                return w;
        }
    }
    float x, y, z, w;
};


struct vec2f {
    float x, y;

    vec2f() : x(0), y(0) {}
    vec2f(float a, float b) : x(a), y(b) {}

    inline vec2f operator*(float t) const {
        return {x*t, y*t};
    }
    inline float norm() const { return std::sqrt(x*x+y*y); }
    inline vec2f operator-(const vec2f& t) const {
        return {x-t.x, y-t.y};
    }
    inline vec2f operator*(const vec2f& t) const {
        return {x*t.x, y * t.y};
    }


};

template<typename T>
vec3<T> abs(vec3<T> a) {
    return vec3<T>(a.x > 0? a.x:-a.x, a.y > 0? a.y:-a.y,  a.z > 0? a.z:-a.z);
}

template<typename T>
inline T dot(vec3<T> a, vec3<T> b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline float dot(vec2f a, vec2f b) {
    return a.x*b.x + a.y*b.y;
}

template<typename T>
inline T clamp(T a, T min, T max) {
    if (a > max) return max;
    if (a < min) return min;
    return a;
}

template<typename T>
inline T sign(T a) {
    if (a > 0) return 1;
    if (a < 0) return -1;
    return 0;
}

template <typename T>
inline vec3<T> mod(vec3<T> a, float t) {
    return vec3<T>(a.x - t * floor(a.x/t), a.y - t * floor(a.y/t), a.z - t * floor(a.z/t));
}

inline vec3f abs(vec3f a) {
    return vec3f(fabsf(a.x), fabsf(a.y), fabsf(a.z));
}
template<typename T>
inline vec3<T> operator-(float a, const vec3<T>& t) {
    return {a-t.x, a-t.y, a-t.z};
}

template <typename T>
inline vec3<T> operator/(float a, const vec3<T>& t) {
    return {a/t.x, a/t.y, a/t.z};
}

inline vec3f pow(const vec3f &a, float t) {
    return {powf(a.x, t), powf(a.y, t), powf(a.z, t)};
}
#endif //RT_MY_MATH_H
