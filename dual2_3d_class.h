#pragma once
#include <type_traits>
#include <iostream>
#include <cmath>
#include <complex>

/*
    Dual2_3d<T> - Mathematical Guide

    This class represents a trivariate function f(x, y, z) and its derivatives up to second order.

    Members:
      val   = f(x, y, z)         // function value
      dx    = df/dx              // first partial derivative w.r.t. x
      dy    = df/dy              // first partial derivative w.r.t. y
      dz    = df/dz              // first partial derivative w.r.t. z
      dxx   = d²f/dx²            // second partial derivative w.r.t. x
	  dxy   = d²f/dxdy           // mixed second partial derivative w.r.t. x and y
	  dxz   = d²f/dxdz           // mixed second partial derivative w.r.t. x and z
      dyy   = d²f/dy²            // second partial derivative w.r.t. y
	  dyz   = d²f/dydz           // mixed second partial derivative w.r.t. y and z
      dzz   = d²f/dz²            // second partial derivative w.r.t. z

    Constructors:
      Dual2_3d()
        - All values and derivatives set to zero.
        - Represents the zero function.

      Dual2_3d(val)
        - Sets f(x, y, z) = val, all derivatives zero.
        - Represents a constant function.

      Dual2_3d(val, dx, dy, dz)
        - Sets f(x, y, z) = val, df/dx = dx, df/dy = dy, df/dz = dz.
        - All second derivatives are zero.
        - Represents a function linear in x, y, z.

      Dual2_3d(val, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz)
        - Sets all values and first/second derivatives explicitly.
        - Represents a function with arbitrary quadratic behavior in x, y, z.

    Note:
      If you omit a derivative in the constructor, it is set to zero by default.
      This is mathematically equivalent to assuming the function does not vary in that direction/order.

      Seeding for Automatic Differentiation:
      - To represent the variable x: Dual2_3d<T>(x0, 1, 0, 0);
      - To represent the variable y: Dual2_3d<T>(y0, 0, 1, 0);
      - To represent the variable z: Dual2_3d<T>(z0, 0, 0, 1);
      - To represent a constant:    Dual2_3d<T>(c0);
      (All second derivatives are zero unless specified.)
*/

// Helper for C++14: common_type_t and conditional_t
template<typename... Ts>
using common_type_t = typename std::common_type<Ts...>::type;
template<bool B, typename T, typename F>
using conditional_t = typename std::conditional<B, T, F>::type;

// Trivariate second-order dual number: tracks value, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz
template<typename T>
struct Dual2_3d {
    T val;   // function value
    T dx;    // df/dx
    T dy;    // df/dy
    T dz;    // df/dz
    T dxx;   // d²f/dx²
    T dxy;   // d²f/dxdy
    T dxz;   // d²f/dxdz
    T dyy;   // d²f/dy²
    T dyz;   // d²f/dydz
    T dzz;   // d²f/dz²

    // Constructors
    Dual2_3d()
        : val(0), dx(0), dy(0), dz(0), dxx(0), dxy(0), dxz(0), dyy(0), dyz(0), dzz(0) {}

    Dual2_3d(const T& v)
        : val(v), dx(0), dy(0), dz(0), dxx(0), dxy(0), dxz(0), dyy(0), dyz(0), dzz(0) {}

    Dual2_3d(const T& v, const T& dx_, const T& dy_, const T& dz_,
        const T& dxx_ = 0, const T& dxy_ = 0, const T& dxz_ = 0,
        const T& dyy_ = 0, const T& dyz_ = 0, const T& dzz_ = 0)
        : val(v), dx(dx_), dy(dy_), dz(dz_),
        dxx(dxx_), dxy(dxy_), dxz(dxz_),
        dyy(dyy_), dyz(dyz_), dzz(dzz_) {}

    // Addition
    Dual2_3d operator+(const Dual2_3d& o) const {
        typedef common_type_t<T, decltype(o.val)> ReturnType;
        return Dual2_3d<ReturnType>(
            val + o.val, dx + o.dx, dy + o.dy, dz + o.dz,
            dxx + o.dxx, dxy + o.dxy, dxz + o.dxz,
            dyy + o.dyy, dyz + o.dyz, dzz + o.dzz
        );
    }

    // Subtraction
    Dual2_3d operator-(const Dual2_3d& o) const {
        typedef common_type_t<T, decltype(o.val)> ReturnType;
        return Dual2_3d<ReturnType>(
            val - o.val, dx - o.dx, dy - o.dy, dz - o.dz,
            dxx - o.dxx, dxy - o.dxy, dxz - o.dxz,
            dyy - o.dyy, dyz - o.dyz, dzz - o.dzz
        );
    }

    // Multiplication
    Dual2_3d operator*(const Dual2_3d& o) const {
        typedef common_type_t<T, decltype(o.val)> ReturnType;
        ReturnType v = val * o.val;
        ReturnType dx_ = dx * o.val + val * o.dx;
        ReturnType dy_ = dy * o.val + val * o.dy;
        ReturnType dz_ = dz * o.val + val * o.dz;
        ReturnType dxx_ = dxx * o.val + 2.0 * dx * o.dx + val * o.dxx;
        ReturnType dxy_ = dxy * o.val + dx * o.dy + dy * o.dx + val * o.dxy;
        ReturnType dxz_ = dxz * o.val + dx * o.dz + dz * o.dx + val * o.dxz;
        ReturnType dyy_ = dyy * o.val + 2.0 * dy * o.dy + val * o.dyy;
        ReturnType dyz_ = dyz * o.val + dy * o.dz + dz * o.dy + val * o.dyz;
        ReturnType dzz_ = dzz * o.val + 2.0 * dz * o.dz + val * o.dzz;
        return Dual2_3d<ReturnType>(v, dx_, dy_, dz_, dxx_, dxy_, dxz_, dyy_, dyz_, dzz_);
    }

    // Division
    Dual2_3d operator/(const Dual2_3d& o) const {
        typedef common_type_t<T, decltype(o.val)> RawType;
        typedef conditional_t<
            std::is_arithmetic<RawType>::value,
            double,
            std::complex<double>
        > ReturnType;
        ReturnType v = static_cast<ReturnType>(val);
        ReturnType dx_ = static_cast<ReturnType>(dx);
        ReturnType dy_ = static_cast<ReturnType>(dy);
        ReturnType dz_ = static_cast<ReturnType>(dz);
        ReturnType dxx_ = static_cast<ReturnType>(dxx);
        ReturnType dxy_ = static_cast<ReturnType>(dxy);
        ReturnType dxz_ = static_cast<ReturnType>(dxz);
        ReturnType dyy_ = static_cast<ReturnType>(dyy);
        ReturnType dyz_ = static_cast<ReturnType>(dyz);
        ReturnType dzz_ = static_cast<ReturnType>(dzz);

        ReturnType ov = static_cast<ReturnType>(o.val);
        ReturnType odx = static_cast<ReturnType>(o.dx);
        ReturnType ody = static_cast<ReturnType>(o.dy);
        ReturnType odz = static_cast<ReturnType>(o.dz);
        ReturnType odxx = static_cast<ReturnType>(o.dxx);
        ReturnType odxy = static_cast<ReturnType>(o.dxy);
        ReturnType odxz = static_cast<ReturnType>(o.dxz);
        ReturnType odyy = static_cast<ReturnType>(o.dyy);
        ReturnType odyz = static_cast<ReturnType>(o.dyz);
        ReturnType odzz = static_cast<ReturnType>(o.dzz);

        ReturnType ov2 = ov * ov;
        ReturnType ov3 = ov2 * ov;

        ReturnType out_val = v / ov;
        ReturnType out_dx = (dx_ * ov - v * odx) / ov2;
        ReturnType out_dy = (dy_ * ov - v * ody) / ov2;
        ReturnType out_dz = (dz_ * ov - v * odz) / ov2;
        ReturnType out_dxx = (dxx_ * ov2 - 2.0 * dx_ * ov * odx - v * ov * odxx + 2.0 * v * odx * odx) / ov3;
        ReturnType out_dxy = (dxy_ * ov2 - dx_ * ov * ody - dy_ * ov * odx - v * ov * odxy + 2.0 * v * odx * ody) / ov3;
        ReturnType out_dxz = (dxz_ * ov2 - dx_ * ov * odz - dz_ * ov * odx - v * ov * odxz + 2.0 * v * odx * odz) / ov3;
        ReturnType out_dyy = (dyy_ * ov2 - 2.0 * dy_ * ov * ody - v * ov * odyy + 2.0 * v * ody * ody) / ov3;
        ReturnType out_dyz = (dyz_ * ov2 - dy_ * ov * odz - dz_ * ov * ody - v * ov * odyz + 2.0 * v * ody * odz) / ov3;
        ReturnType out_dzz = (dzz_ * ov2 - 2.0 * dz_ * ov * odz - v * ov * odzz + 2.0 * v * odz * odz) / ov3;

        return Dual2_3d<ReturnType>(out_val, out_dx, out_dy, out_dz,
            out_dxx, out_dxy, out_dxz, out_dyy, out_dyz, out_dzz);
    }

    // Scalar right-hand side
    template<typename U>
    auto operator+(U o) const -> Dual2_3d<common_type_t<T, U>> {
        return Dual2_3d<common_type_t<T, U>>(val + o, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz);
    }
    template<typename U>
    auto operator-(U o) const -> Dual2_3d<common_type_t<T, U>> {
        return Dual2_3d<common_type_t<T, U>>(val - o, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz);
    }
    template<typename U>
    auto operator*(U o) const -> Dual2_3d<common_type_t<T, U>> {
        return Dual2_3d<common_type_t<T, U>>(val * o, dx * o, dy * o, dz * o, dxx * o, dxy * o, dxz * o, dyy * o, dyz * o, dzz * o);
    }
    template<typename U>
    auto operator/(U o) const -> Dual2_3d<conditional_t<
        std::is_arithmetic<common_type_t<T, U>>::value,
        double,
        std::complex<double>
    >> {
        typedef common_type_t<T, U> RawType;
        typedef conditional_t<
            std::is_arithmetic<RawType>::value,
            double,
            std::complex<double>
        > ReturnType;
        ReturnType o_ = static_cast<ReturnType>(o);
        return Dual2_3d<ReturnType>(
            val / o_,
            dx / o_,
            dy / o_,
            dz / o_,
            dxx / o_,
            dxy / o_,
            dxz / o_,
            dyy / o_,
            dyz / o_,
            dzz / o_
        );
    }

    // Output
    friend std::ostream& operator<<(std::ostream& os, const Dual2_3d& x) {
        os << x.val
            << " (dx=" << x.dx << ", dy=" << x.dy << ", dz=" << x.dz
            << ", dxx=" << x.dxx << ", dxy=" << x.dxy << ", dxz=" << x.dxz
            << ", dyy=" << x.dyy << ", dyz=" << x.dyz << ", dzz=" << x.dzz << ")";
        return os;
    }
};

// --- Dual2_3d <op> scalar (left hand side) ---
template<typename T, typename O>
auto operator+(O lhs, const Dual2_3d<T>& rhs) -> Dual2_3d<common_type_t<O, T>> {
    return Dual2_3d<common_type_t<O, T>>(lhs + rhs.val, rhs.dx, rhs.dy, rhs.dz, rhs.dxx, rhs.dxy, rhs.dxz, rhs.dyy, rhs.dyz, rhs.dzz);
}
template<typename T, typename O>
auto operator-(O lhs, const Dual2_3d<T>& rhs) -> Dual2_3d<common_type_t<O, T>> {
    return Dual2_3d<common_type_t<O, T>>(lhs - rhs.val, -rhs.dx, -rhs.dy, -rhs.dz, -rhs.dxx, -rhs.dxy, -rhs.dxz, -rhs.dyy, -rhs.dyz, -rhs.dzz);
}
template<typename T, typename O>
auto operator*(O lhs, const Dual2_3d<T>& rhs) -> Dual2_3d<common_type_t<O, T>> {
    return Dual2_3d<common_type_t<O, T>>(lhs * rhs.val, lhs * rhs.dx, lhs * rhs.dy, lhs * rhs.dz, lhs * rhs.dxx, lhs * rhs.dxy, lhs * rhs.dxz, lhs * rhs.dyy, lhs * rhs.dyz, lhs * rhs.dzz);
}
template<typename T, typename O>
auto operator/(O lhs, const Dual2_3d<T>& rhs) -> Dual2_3d<conditional_t<
    std::is_arithmetic<common_type_t<O, T>>::value,
    double,
    std::complex<double>
>> {
    typedef common_type_t<O, T> RawType;
    typedef conditional_t<
        std::is_arithmetic<RawType>::value,
        double,
        std::complex<double>
    > ReturnType;
    ReturnType lhs_ = static_cast<ReturnType>(lhs);
    ReturnType v = static_cast<ReturnType>(rhs.val);
    ReturnType dx = static_cast<ReturnType>(rhs.dx);
    ReturnType dy = static_cast<ReturnType>(rhs.dy);
    ReturnType dz = static_cast<ReturnType>(rhs.dz);
    ReturnType dxx = static_cast<ReturnType>(rhs.dxx);
    ReturnType dxy = static_cast<ReturnType>(rhs.dxy);
    ReturnType dxz = static_cast<ReturnType>(rhs.dxz);
    ReturnType dyy = static_cast<ReturnType>(rhs.dyy);
    ReturnType dyz = static_cast<ReturnType>(rhs.dyz);
    ReturnType dzz = static_cast<ReturnType>(rhs.dzz);

    ReturnType v2 = v * v;
    ReturnType v3 = v2 * v;

    ReturnType out_val = lhs_ / v;
    ReturnType out_dx = -lhs_ * dx / v2;
    ReturnType out_dy = -lhs_ * dy / v2;
    ReturnType out_dz = -lhs_ * dz / v2;
    ReturnType out_dxx = (2.0 * lhs_ * dx * dx - lhs_ * v * dxx) / v3;
    ReturnType out_dxy = (2.0 * lhs_ * dx * dy - lhs_ * v * dxy) / v3;
    ReturnType out_dxz = (2.0 * lhs_ * dx * dz - lhs_ * v * dxz) / v3;
    ReturnType out_dyy = (2.0 * lhs_ * dy * dy - lhs_ * v * dyy) / v3;
    ReturnType out_dyz = (2.0 * lhs_ * dy * dz - lhs_ * v * dyz) / v3;
    ReturnType out_dzz = (2.0 * lhs_ * dz * dz - lhs_ * v * dzz) / v3;

    return Dual2_3d<ReturnType>(out_val, out_dx, out_dy, out_dz, out_dxx, out_dxy, out_dxz, out_dyy, out_dyz, out_dzz);
}

// Free functions for common operations
// pow(x, scalar): x^exponent for Dual2_3d
template<typename T, typename U>
auto pow(const Dual2_3d<T>& x, U exponent) -> Dual2_3d<common_type_t<T, U>> {
    using std::pow;
    typedef common_type_t<T, U> ReturnType;

    ReturnType v = pow(x.val, exponent);
    ReturnType dx = exponent * pow(x.val, exponent - 1) * x.dx;
    ReturnType dy = exponent * pow(x.val, exponent - 1) * x.dy;
    ReturnType dz = exponent * pow(x.val, exponent - 1) * x.dz;
    ReturnType dxx = exponent * (exponent - 1) * pow(x.val, exponent - 2) * x.dx * x.dx
                   + exponent * pow(x.val, exponent - 1) * x.dxx;
    ReturnType dxy = exponent * (exponent - 1) * pow(x.val, exponent - 2) * x.dx * x.dy
                   + exponent * pow(x.val, exponent - 1) * x.dxy;
    ReturnType dxz = exponent * (exponent - 1) * pow(x.val, exponent - 2) * x.dx * x.dz
                   + exponent * pow(x.val, exponent - 1) * x.dxz;
    ReturnType dyy = exponent * (exponent - 1) * pow(x.val, exponent - 2) * x.dy * x.dy
                   + exponent * pow(x.val, exponent - 1) * x.dyy;
    ReturnType dyz = exponent * (exponent - 1) * pow(x.val, exponent - 2) * x.dy * x.dz
                   + exponent * pow(x.val, exponent - 1) * x.dyz;
    ReturnType dzz = exponent * (exponent - 1) * pow(x.val, exponent - 2) * x.dz * x.dz
                   + exponent * pow(x.val, exponent - 1) * x.dzz;
    return Dual2_3d<ReturnType>(v, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz);
}

// pow(scalar, x): exponent^x for Dual2_3d
template<typename T, typename U>
auto pow(U scalar, const Dual2_3d<T>& x) -> Dual2_3d<common_type_t<T, U>> {
    using std::pow;
    using std::log;
    typedef common_type_t<T, U> ReturnType;
    ReturnType v = pow(scalar, x.val);
    ReturnType logexponent = log(scalar);
    // First derivatives
    ReturnType dx = v * x.dx * logexponent;
    ReturnType dy = v * x.dy * logexponent;
    ReturnType dz = v * x.dz * logexponent;
    // Second derivatives
    ReturnType dxx = v * (x.dxx * logexponent + x.dx * x.dx * logexponent * logexponent);
    ReturnType dxy = v * (x.dxy * logexponent + x.dx * x.dy * logexponent * logexponent);
    ReturnType dxz = v * (x.dxz * logexponent + x.dx * x.dz * logexponent * logexponent);
    ReturnType dyy = v * (x.dyy * logexponent + x.dy * x.dy * logexponent * logexponent);
    ReturnType dyz = v * (x.dyz * logexponent + x.dy * x.dz * logexponent * logexponent);
    ReturnType dzz = v * (x.dzz * logexponent + x.dz * x.dz * logexponent * logexponent);
    return Dual2_3d<ReturnType>(v, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz);
}

// pow(x, q): x^y for Dual2_3d
template<typename T>
auto pow(const Dual2_3d<T>& x, const Dual2_3d<T>& q) -> Dual2_3d<common_type_t<T, decltype(std::pow(x.val, q.val))>> {
    using std::pow;
    using std::log;
    typedef common_type_t<T, decltype(pow(x.val, q.val))> ReturnType;
    ReturnType v = pow(x.val, q.val);
    ReturnType logx = log(x.val);
    ReturnType inv_x =  1.0 / x.val;
    ReturnType inv_x_sq = inv_x * inv_x;

    // First derivatives
    ReturnType dx = v * (q.dx * logx + q.val * x.dx * inv_x);
    ReturnType dy = v * (q.dy * logx + q.val * x.dy * inv_x);
    ReturnType dz = v * (q.dz * logx + q.val * x.dz * inv_x);

    // Second derivatives
    ReturnType dxx = v * (
        q.val * inv_x_sq * (x.val * x.dxx - x.dx * x.dx) +
        q.dx * x.dx * inv_x + q.dx * x.dx * inv_x + q.dxx * logx +
        (q.val * x.dx * inv_x + q.dx * logx) * (q.val * x.dx * inv_x + q.dx * logx)
        );
    ReturnType dyy = v * (
        q.val * inv_x_sq * (x.val * x.dyy - x.dy * x.dy) +
        q.dy * x.dy * inv_x + q.dy * x.dy * inv_x + q.dyy * logx +
        (q.val * x.dy * inv_x + q.dy * logx) * (q.val * x.dy * inv_x + q.dy * logx)
        );
    ReturnType dzz = v * (
        q.val * inv_x_sq * (x.val * x.dzz - x.dz * x.dz) +
        q.dz * x.dz * inv_x + q.dz * x.dz * inv_x + q.dzz * logx +
        (q.val * x.dz * inv_x + q.dz * logx) * (q.val * x.dz * inv_x + q.dz * logx)
        );
    ReturnType dxy = v * (
        q.val * inv_x_sq * (x.val * x.dxy - x.dx * x.dy) +
        q.dy * x.dx * inv_x + q.dx * x.dy * inv_x + q.dxy * logx +
        (q.val * x.dy * inv_x + q.dy * logx) * (q.val * x.dx * inv_x + q.dx * logx)
    );
    ReturnType dxz = v * (
        q.val * inv_x_sq * (x.val * x.dxz - x.dx * x.dz) +
        q.dz * x.dx * inv_x + q.dx * x.dz * inv_x + q.dxz * logx +
        (q.val * x.dz * inv_x + q.dz * logx) * (q.val * x.dx * inv_x + q.dx * logx)
        );
    ReturnType dyz = v * (
        q.val * inv_x_sq * (x.val * x.dyz - x.dy * x.dz) +
        q.dz * x.dy * inv_x + q.dy * x.dz * inv_x + q.dyz * logx +
        (q.val * x.dz * inv_x + q.dz * logx) * (q.val * x.dy * inv_x + q.dy * logx)
        );

    return Dual2_3d<ReturnType>(v, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz);
}

// Exponential for Dual2_3d
template<typename T>
auto exp(const Dual2_3d<T>& x) -> Dual2_3d<common_type_t<T, decltype(std::exp(x.val))>> {
    typedef common_type_t<T, decltype(std::exp(x.val))> ReturnType;
    ReturnType v = std::exp(x.val);
    ReturnType dx = v * x.dx;
    ReturnType dy = v * x.dy;
    ReturnType dz = v * x.dz;
    ReturnType dxx = v * (x.dxx + x.dx * x.dx);
    ReturnType dxy = v * (x.dxy + x.dx * x.dy);
    ReturnType dxz = v * (x.dxz + x.dx * x.dz);
    ReturnType dyy = v * (x.dyy + x.dy * x.dy);
    ReturnType dyz = v * (x.dyz + x.dy * x.dz);
    ReturnType dzz = v * (x.dzz + x.dz * x.dz);
    return Dual2_3d<ReturnType>(v, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz);
}

// Logarithm for Dual2_3d
template<typename T>
auto log(const Dual2_3d<T>& x) -> Dual2_3d<common_type_t<T, decltype(std::log(x.val))>> {
    typedef common_type_t<T, decltype(std::log(x.val))> ReturnType;
    ReturnType v = std::log(x.val);
    ReturnType inv_x = 1.0 / x.val;
    ReturnType inv_x_sq = inv_x * inv_x;
    ReturnType dx = x.dx * inv_x;
    ReturnType dy = x.dy * inv_x;
    ReturnType dz = x.dz * inv_x;
    ReturnType dxx = (x.dxx * x.val - x.dx * x.dx) * inv_x_sq;
    ReturnType dxy = (x.dxy * x.val - x.dx * x.dy) * inv_x_sq;
    ReturnType dxz = (x.dxz * x.val - x.dx * x.dz) * inv_x_sq;
    ReturnType dyy = (x.dyy * x.val - x.dy * x.dy) * inv_x_sq;
    ReturnType dyz = (x.dyz * x.val - x.dy * x.dz) * inv_x_sq;
    ReturnType dzz = (x.dzz * x.val - x.dz * x.dz) * inv_x_sq;
    return Dual2_3d<ReturnType>(v, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz);
}

// Square root for Dual2_3d
template<typename T>
auto sqrt(const Dual2_3d<T>& x) -> Dual2_3d<common_type_t<T, decltype(x.val)>> {
    typedef common_type_t<T, decltype(x.val)> ReturnType;
    ReturnType v = std::sqrt(x.val);
    ReturnType dx = 0.5 / v * x.dx;
    ReturnType dy = 0.5 / v * x.dy;
    ReturnType dz = 0.5 / v * x.dz;
    ReturnType dxx = 0.5 / v * x.dxx - 0.25 * x.dx * x.dx / (v * v * v);
    ReturnType dxy = 0.5 / v * x.dxy - 0.25 * x.dx * x.dy / (v * v * v);
    ReturnType dxz = 0.5 / v * x.dxz - 0.25 * x.dx * x.dz / (v * v * v);
    ReturnType dyy = 0.5 / v * x.dyy - 0.25 * x.dy * x.dy / (v * v * v);
    ReturnType dyz = 0.5 / v * x.dyz - 0.25 * x.dy * x.dz / (v * v * v);
    ReturnType dzz = 0.5 / v * x.dzz - 0.25 * x.dz * x.dz / (v * v * v);
    return Dual2_3d<ReturnType>(v, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz);
}

// Sine for Dual2_3d
template<typename T>
auto sin(const Dual2_3d<T>& x) -> Dual2_3d<common_type_t<T, decltype(std::sin(x.val))>> {
    typedef common_type_t<T, decltype(std::sin(x.val))> ReturnType;
    ReturnType v = std::sin(x.val);
    ReturnType dx = std::cos(x.val) * x.dx;
    ReturnType dy = std::cos(x.val) * x.dy;
    ReturnType dz = std::cos(x.val) * x.dz;
    ReturnType dxx = -std::sin(x.val) * x.dx * x.dx + std::cos(x.val) * x.dxx;
    ReturnType dxy = -std::sin(x.val) * x.dx * x.dy + std::cos(x.val) * x.dxy;
    ReturnType dxz = -std::sin(x.val) * x.dx * x.dz + std::cos(x.val) * x.dxz;
    ReturnType dyy = -std::sin(x.val) * x.dy * x.dy + std::cos(x.val) * x.dyy;
    ReturnType dyz = -std::sin(x.val) * x.dy * x.dz + std::cos(x.val) * x.dyz;
    ReturnType dzz = -std::sin(x.val) * x.dz * x.dz + std::cos(x.val) * x.dzz;
    return Dual2_3d<ReturnType>(v, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz);
}

// Cosine for Dual2_3d
template<typename T>
auto cos(const Dual2_3d<T>& x) -> Dual2_3d<common_type_t<T, decltype(std::cos(x.val))>> {
    typedef common_type_t<T, decltype(std::cos(x.val))> ReturnType;
    ReturnType v = std::cos(x.val);
    ReturnType dx = -std::sin(x.val) * x.dx;
    ReturnType dy = -std::sin(x.val) * x.dy;
    ReturnType dz = -std::sin(x.val) * x.dz;
    ReturnType dxx = -std::cos(x.val) * x.dx * x.dx - std::sin(x.val) * x.dxx;
    ReturnType dxy = -std::cos(x.val) * x.dx * x.dy - std::sin(x.val) * x.dxy;
    ReturnType dxz = -std::cos(x.val) * x.dx * x.dz - std::sin(x.val) * x.dxz;
    ReturnType dyy = -std::cos(x.val) * x.dy * x.dy - std::sin(x.val) * x.dyy;
    ReturnType dyz = -std::cos(x.val) * x.dy * x.dz - std::sin(x.val) * x.dyz;
    ReturnType dzz = -std::cos(x.val) * x.dz * x.dz - std::sin(x.val) * x.dzz;
    return Dual2_3d<ReturnType>(v, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz);
}

// Tangent for Dual2_3d
template<typename T>
auto tan(const Dual2_3d<T>& x) -> Dual2_3d<common_type_t<T, decltype(std::tan(x.val))>> {
    typedef common_type_t<T, decltype(std::tan(x.val))> ReturnType;
    ReturnType v = std::tan(x.val);
    ReturnType dx = (1.0 / (std::cos(x.val) * std::cos(x.val))) * x.dx;
    ReturnType dy = (1.0 / (std::cos(x.val) * std::cos(x.val))) * x.dy;
    ReturnType dz = (1.0 / (std::cos(x.val) * std::cos(x.val))) * x.dz;
    ReturnType dxx = (2.0 * std::sin(x.val) / (std::cos(x.val) * std::cos(x.val) * std::cos(x.val))) * x.dx * x.dx + (1.0 / (std::cos(x.val) * std::cos(x.val))) * x.dxx;
    ReturnType dxy = (2.0 * std::sin(x.val) / (std::cos(x.val) * std::cos(x.val) * std::cos(x.val))) * x.dx * x.dy + (1.0 / (std::cos(x.val) * std::cos(x.val))) * x.dxy;
    ReturnType dxz = (2.0 * std::sin(x.val) / (std::cos(x.val) * std::cos(x.val) * std::cos(x.val))) * x.dx * x.dz + (1.0 / (std::cos(x.val) * std::cos(x.val))) * x.dxz;
    ReturnType dyy = (2.0 * std::sin(x.val) / (std::cos(x.val) * std::cos(x.val) * std::cos(x.val))) * x.dy * x.dy + (1.0 / (std::cos(x.val) * std::cos(x.val))) * x.dyy;
    ReturnType dyz = (2.0 * std::sin(x.val) / (std::cos(x.val) * std::cos(x.val) * std::cos(x.val))) * x.dy * x.dz + (1.0 / (std::cos(x.val) * std::cos(x.val))) * x.dyz;
    ReturnType dzz = (2.0 * std::sin(x.val) / (std::cos(x.val) * std::cos(x.val) * std::cos(x.val))) * x.dz * x.dz + (1.0 / (std::cos(x.val) * std::cos(x.val))) * x.dzz;
    return Dual2_3d<ReturnType>(v, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz);
}
