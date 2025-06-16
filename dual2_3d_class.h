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
      dx    = ∂f/∂x              // first partial derivative w.r.t. x
      dy    = ∂f/∂y              // first partial derivative w.r.t. y
      dz    = ∂f/∂z              // first partial derivative w.r.t. z
      dxx   = ∂²f/∂x²            // second partial derivative w.r.t. x
	  dxy   = ∂²f/∂x∂y           // mixed second partial derivative w.r.t. x and y
	  dxz   = ∂²f/∂x∂z           // mixed second partial derivative w.r.t. x and z
      dyy   = ∂²f/∂y²            // second partial derivative w.r.t. y
	  dyz   = ∂²f/∂y∂z           // mixed second partial derivative w.r.t. y and z
      dzz   = ∂²f/∂z²            // second partial derivative w.r.t. z

    Constructors:
      Dual2_3d()
        - All values and derivatives set to zero.
        - Represents the zero function.

      Dual2_3d(val)
        - Sets f(x, y, z) = val, all derivatives zero.
        - Represents a constant function.

      Dual2_3d(val, dx, dy, dz)
        - Sets f(x, y, z) = val, ∂f/∂x = dx, ∂f/∂y = dy, ∂f/∂z = dz.
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

// Trivariate second-order dual number: tracks value, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz
template<typename T>
struct Dual2_3d {
    T val;   // function value
    T dx;    // ∂f/∂x
    T dy;    // ∂f/∂y
    T dz;    // ∂f/∂z
    T dxx;   // ∂²f/∂x²
    T dxy;   // ∂²f/∂x∂y
    T dxz;   // ∂²f/∂x∂z
    T dyy;   // ∂²f/∂y²
    T dyz;   // ∂²f/∂y∂z
    T dzz;   // ∂²f/∂z²

    // Constructors
    Dual2_3d()
        : val(0), dx(0), dy(0), dz(0), dxx(0), dxy(0), dxz(0), dyy(0), dyz(0), dzz(0) {
    }

    Dual2_3d(const T& v)
        : val(v), dx(0), dy(0), dz(0), dxx(0), dxy(0), dxz(0), dyy(0), dyz(0), dzz(0) {
    }

	// Constructors with first and second derivatives
    Dual2_3d(const T& v, const T& dx_, const T& dy_, const T& dz_,
        const T& dxx_ = 0, const T& dxy_ = 0, const T& dxz_ = 0,
        const T& dyy_ = 0, const T& dyz_ = 0, const T& dzz_ = 0)
        : val(v), dx(dx_), dy(dy_), dz(dz_),
        dxx(dxx_), dxy(dxy_), dxz(dxz_),
        dyy(dyy_), dyz(dyz_), dzz(dzz_) {
    }

    // Addition
    Dual2_3d operator+(const Dual2_3d& o) const {
        using ReturnType = std::common_type_t<T, decltype(o.val)>;
        return Dual2_3d<ReturnType>(
            val + o.val,
            dx + o.dx,
            dy + o.dy,
            dz + o.dz,
            dxx + o.dxx,
            dxy + o.dxy,
            dxz + o.dxz,
            dyy + o.dyy,
            dyz + o.dyz,
            dzz + o.dzz
        );
    }

    // Subtraction
    Dual2_3d operator-(const Dual2_3d& o) const {
        using ReturnType = std::common_type_t<T, decltype(o.val)>;
        return Dual2_3d<ReturnType>(
            val - o.val,
            dx - o.dx,
            dy - o.dy,
            dz - o.dz,
            dxx - o.dxx,
            dxy - o.dxy,
            dxz - o.dxz,
            dyy - o.dyy,
            dyz - o.dyz,
            dzz - o.dzz
        );
    }

    // Multiplication
    Dual2_3d operator*(const Dual2_3d& o) const {
        using ReturnType = std::common_type_t<T, decltype(o.val)>;
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
        using RawType = std::common_type_t<T, decltype(o.val)>;
        using ReturnType = std::conditional_t<
            std::is_arithmetic_v<RawType>,
            double,
            std::complex<double>
        >;
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
    auto operator+(U o) const {
        using ReturnType = std::common_type_t<T, U>;
        return Dual2_3d<ReturnType>(val + o, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz);
    }
    template<typename U>
    auto operator-(U o) const {
        using ReturnType = std::common_type_t<T, U>;
        return Dual2_3d<ReturnType>(val - o, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz);
    }
    template<typename U>
    auto operator*(U o) const {
        using ReturnType = std::common_type_t<T, U>;
        return Dual2_3d<ReturnType>(val * o, dx * o, dy * o, dz * o, dxx * o, dxy * o, dxz * o, dyy * o, dyz * o, dzz * o);
    }
    template<typename U>
    auto operator/(U o) const {
        using RawType = std::common_type_t<T, U>;
        using ReturnType = std::conditional_t<
            std::is_arithmetic_v<RawType>,
            double,
            std::complex<double>
        >;
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

// Deduction guide
template<typename V, typename DX, typename DY, typename DZ, typename DXX, typename DXY, typename DXZ, typename DYY, typename DYZ, typename DZZ>
Dual2_3d(V, DX, DY, DZ, DXX, DXY, DXZ, DYY, DYZ, DZZ) -> Dual2_3d<std::common_type_t<V, DX, DY, DZ, DXX, DXY, DXZ, DYY, DYZ, DZZ>>;

// --- Dual2_3d <op> scalar (left hand side) ---
template<typename T, typename O>
auto operator+(O lhs, const Dual2_3d<T>& rhs) {
    using ReturnType = std::common_type_t<O, T>;
    return Dual2_3d<ReturnType>(lhs + rhs.val, rhs.dx, rhs.dy, rhs.dz, rhs.dxx, rhs.dxy, rhs.dxz, rhs.dyy, rhs.dyz, rhs.dzz);
}
template<typename T, typename O>
auto operator-(O lhs, const Dual2_3d<T>& rhs) {
    using ReturnType = std::common_type_t<O, T>;
    return Dual2_3d<ReturnType>(lhs - rhs.val, -rhs.dx, -rhs.dy, -rhs.dz, -rhs.dxx, -rhs.dxy, -rhs.dxz, -rhs.dyy, -rhs.dyz, -rhs.dzz);
}
template<typename T, typename O>
auto operator*(O lhs, const Dual2_3d<T>& rhs) {
    using ReturnType = std::common_type_t<O, T>;
    return Dual2_3d<ReturnType>(lhs * rhs.val, lhs * rhs.dx, lhs * rhs.dy, lhs * rhs.dz, lhs * rhs.dxx, lhs * rhs.dxy, lhs * rhs.dxz, lhs * rhs.dyy, lhs * rhs.dyz, lhs * rhs.dzz);
}
template<typename T, typename O>
auto operator/(O lhs, const Dual2_3d<T>& rhs) {
    using RawType = std::common_type_t<O, T>;
    using ReturnType = std::conditional_t<
        std::is_arithmetic_v<RawType>,
        double,
        std::complex<double>
    >;
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
auto pow(const Dual2_3d<T>& x, U exponent) {
    using std::pow;
    using ReturnType = std::common_type_t<T, U>;
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

// pow(x, y): x^y for Dual2_3d
template<typename T>
auto pow(const Dual2_3d<T>& x, const Dual2_3d<T>& y) {
    using std::pow;
    using std::log;
    using ReturnType = std::common_type_t<T, decltype(pow(x.val, y.val))>;
    ReturnType v = pow(x.val, y.val);
    ReturnType logx = log(x.val);

    // First derivatives
    ReturnType dx = v * (y.dx * logx + y.val * x.dx / x.val);
    ReturnType dy = v * (y.dy * logx + y.val * x.dy / x.val);
    ReturnType dz = v * (y.dz * logx + y.val * x.dz / x.val);

    // Second derivatives
    ReturnType dxx = v * (
        y.dxx * logx +
        2.0 * y.dx * x.dx / x.val +
        y.val * (x.dxx * x.val - x.dx * x.dx) / (x.val * x.val)
    );
    ReturnType dyy = v * (
        y.dyy * logx +
        2.0 * y.dy * x.dy / x.val +
        y.val * (x.dyy * x.val - x.dy * x.dy) / (x.val * x.val)
    );
    ReturnType dzz = v * (
        y.dzz * logx +
        2.0 * y.dz * x.dz / x.val +
        y.val * (x.dzz * x.val - x.dz * x.dz) / (x.val * x.val)
    );
    ReturnType dxy = v * (
        y.dxy * logx +
        y.dx * x.dy / x.val +
        y.dy * x.dx / x.val +
        y.val * (x.dxy * x.val - x.dx * x.dy) / (x.val * x.val)
    );
    ReturnType dxz = v * (
        y.dxz * logx +
        y.dx * x.dz / x.val +
        y.dz * x.dx / x.val +
        y.val * (x.dxz * x.val - x.dx * x.dz) / (x.val * x.val)
    );
    ReturnType dyz = v * (
        y.dyz * logx +
        y.dy * x.dz / x.val +
        y.dz * x.dy / x.val +
        y.val * (x.dyz * x.val - x.dy * x.dz) / (x.val * x.val)
    );

    return Dual2_3d<ReturnType>(v, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dzz);
}

template<typename T>
auto exp(const Dual2_3d<T>& x) {
    using ReturnType = std::common_type_t<T, decltype(std::exp(x.val))>;
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

template<typename T>
auto sqrt(const Dual2_3d<T>& x) {
    using ReturnType = std::common_type_t<T, decltype(x.val)>;
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