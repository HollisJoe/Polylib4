/*
 * Polylib - Polygon Management Library
 *
 * Copyright (c) 2010-2011 VCAD System Research Program, RIKEN.
 * All rights reserved.
 *
 * Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef vec2_h
#define vec2_h

#include <iostream>
#include <math.h>

namespace PolylibNS {

// precision
#ifdef _REAL_IS_DOUBLE_
#define REAL_TYPE double
#else
  /** 実数型の指定
   * - デフォルトでは、REAL_TYPE=float
   * - コンパイル時オプション-D_REAL_IS_DOUBLE_を付与することで
   *   REAL_TYPE=doubleになる
   */
#define REAL_TYPE float
#endif

////////////////////////////////////////////////////////////////////////////
///
/// クラス:Vec2<T>
///
////////////////////////////////////////////////////////////////////////////
template<typename T>
class Vec2 {
public:
    Vec2(T v = 0)       { x = y = v; }
    Vec2(T _x, T _y)    { x = _x; y = _y; }
    Vec2(const T v[2])  { x = v[0]; y = v[1]; }
    // use default copy constructor
    // use default operator =()

    Vec2<T>& assign(T _x, T _y) { x = _x; y = _y; return *this; }

    operator       T*()       { return &x; }
    operator const T*() const { return &x; }
          T* ptr()       { return &x; }
    const T* ptr() const { return &x; }
          T& operator [](int i)       { return (&x)[i]; }
    const T& operator [](int i) const { return (&x)[i]; }

    Vec2<T>& operator+=(const Vec2<T>& v) {
        x += v.x; y += v.y;
        return *this;
    }
    Vec2<T>& operator-=(const Vec2<T>& v) {
        x -= v.x; y -= v.y;
        return *this;
    }
    Vec2<T>& operator*=(const Vec2<T>& v) {
        x *= v.x; y *= v.y;
        return *this;
    }
    Vec2<T>& operator/=(const Vec2<T>& v) {
        x /= v.x; y /= v.y;
        return *this;
    }
    Vec2<T>& operator*=(T s) {
        x *= s; y *= s;
        return *this;
    }
    Vec2<T>& operator/=(T s) {
        T inv = 1./s;
        x *= inv; y *= inv;
        return *this;
    }

    Vec2<T> operator+(const Vec2<T>& v) const {
        return Vec2<T>(x + v.x, y + v.y);
    }
    Vec2<T> operator-(const Vec2<T>& v) const {
        return Vec2<T>(x - v.x, y - v.y);
    }
    Vec2<T> operator*(const Vec2<T>& v) const {
        return Vec2<T>(x * v.x, y * v.y);
    }
    Vec2<T> operator/(const Vec2<T>& v) const {
        return Vec2<T>(x / v.x, y / v.y);
    }
    Vec2<T> operator*(T s) const {
        return Vec2<T>(x * s, y * s);
    }
    Vec2<T> operator/(T s) const {
        T inv = 1./s;
        return Vec2<T>(x * inv, y * inv);
    }
    Vec2<T> operator-() const {
        return Vec2<T>(-x, -y);
    }

    bool operator==(const Vec2<T>& v) const {
        return x == v.x && y == v.y;
    }
    bool operator!=(const Vec2<T>& v) const {
        return !(*this == v);
    }

    static Vec2<T> xaxis() { return Vec2<T>(1, 0); }
    static Vec2<T> yaxis() { return Vec2<T>(0, 1); }

    REAL_TYPE lengthSquared() const { return x*x + y*y; }
    REAL_TYPE length() const { return sqrtf(lengthSquared()); }
    Vec2<T>& normalize() {
        REAL_TYPE len = length();
        if (len != 0) {
            return *this /= len;
        }
        else {
            return *this;
        }
    }
    Vec2<T>& normalize(REAL_TYPE* len) {
        *len = length();
        if (*len != 0) {
            return *this /= *len;
        }
        else {
            return *this;
        }
    }
    REAL_TYPE average() const { return (x + y)/(REAL_TYPE)2.0; }

    T x, y;
};

//=========================================================================
// typedef
//=========================================================================
typedef Vec2<unsigned char> Vec2uc;
typedef Vec2<int>           Vec2i;
typedef Vec2<float>         Vec2f;
typedef Vec2<REAL_TYPE>     Vec2r;

//=========================================================================
// inline
//=========================================================================
inline Vec2<REAL_TYPE> operator*(REAL_TYPE s, const Vec2<REAL_TYPE>& v) {
    return Vec2<REAL_TYPE>(s*v.x, s*v.y);
}

inline REAL_TYPE distanceSquared(const Vec2<REAL_TYPE>& a, const Vec2<REAL_TYPE>& b) {
    return (a - b).lengthSquared();
}
inline REAL_TYPE distance(const Vec2<REAL_TYPE>& a, const Vec2<REAL_TYPE>& b) {
    return (a - b).length();
}

template<typename T>
inline std::istream& operator>>(std::istream& is, Vec2<T>& v) {
    return is >> v[0] >> v[1];
}
template<typename T>
inline std::ostream& operator<<(std::ostream& os, const Vec2<T>& v) {
    return os << v[0] << " " << v[1];
}
inline std::istream& operator>>(std::istream& is, Vec2uc& v) {
    int x[2];
    is >> x[0] >> x[1];
    v[0] = x[0]; v[1] = x[1];
    return is;
}
inline std::ostream& operator<<(std::ostream& os, const Vec2uc& v) {
    int x[2];
    x[0] = v[0]; x[1] = v[1];
    os << x[0] << " " << x[1];
    return os;
}

inline bool lessVec2f(const Vec2<REAL_TYPE>& a, const Vec2<REAL_TYPE>& b) {
    return (a.lengthSquared() < b.lengthSquared()) ? true : false;
}

} //namespace PolylibNS

#endif  // vec2_h

