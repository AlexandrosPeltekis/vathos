#include <cmath>
#include <complex>
#include <array>
#include "dualbase_class.h"

// C++14 compatibility helpers
template<typename... Ts>
using common_type_t = typename std::common_type<Ts...>::type;
template<bool B, typename T, typename F>
using conditional_t = typename std::conditional<B, T, F>::type;

/*
    Dual1<T> - Mathematical Guide

    This class represents a univariate function f(x) and its first derivative.

    Members:
      val   = f(p)        // function value
      d1    = df/dp      // first derivative with respect to a parameter

    Constructors:
      Dual1()
        - Sets value and derivative to zero.
        - Represents the zero function.

      Dual1(val)
        - Sets f(p) = val, derivative zero.
        - Represents a constant function.

      Dual1(val, d1)
        - Sets f(p) = val, df/dp = d1.
        - Represents a function linear in p.

    Note:
      If you omit the derivative in the constructor, it is set to zero by default.
      This is mathematically equivalent to assuming the function does not vary in x.

    Seeding for Automatic Differentiation:
      - To represent the variable p: Dual1<T>(p0, 1);
      - To represent a constant:    Dual1<T>(c0);

    Directional Derivative:
      - If you set d1 = 1 for all variables (e.g., Dual1<T>(p0, 1)), then after evaluating a function f, 
      the d1 member will contain the directional derivative of f in the direction where all variables increase 
      equally (i.e., the sum of partial derivatives with respect to each variable).
	  i.e ., d1 = df/dx + df/dy + df/dz for a function f(x, y, z) evaluated at (x0, y0, z0).
	  and the seeding looks like this:
      Dual1<T>(x0, 1); // for x
      Dual1<T>(y0, 1); // for y
      Dual1<T>(z0, 1); // for z
*/
// First-order Dual number: inherits from DualBase and adds first derivative
template<typename T>
struct Dual1 : public DualBase<T> {
    T d1;

	// Constructors
    Dual1() : DualBase<T>(), d1(1) {}
    Dual1(const T& v) : DualBase<T>(v), d1(1) {}
    Dual1(const T& v, const T& d) : DualBase<T>(v), d1(d) {}

    // Conversion constructor: Dual1<L> -> Dual1<std::complex<T>>
    template<typename L>
    explicit Dual1(const Dual1<L>& other,
        typename std::enable_if<!std::is_same<T, L>::value, int>::type = 0)
        : DualBase<T>(T(other.val)), d1(T(other.d1)) {}

    // Member: Dual1 <op> Dual1 (Right hand side)
    Dual1 operator+(const Dual1& o) const {
        typedef common_type_t<T, decltype(o.val)> ReturnType;
        return Dual1<ReturnType>(this->val + o.val, d1 + o.d1);
    }
    Dual1 operator-(const Dual1& o) const {
        typedef common_type_t<T, decltype(o.val)> ReturnType;
        return Dual1<ReturnType>(this->val - o.val, d1 - o.d1);
    }
    Dual1 operator*(const Dual1& o) const {
        typedef common_type_t<T, decltype(o.val)> ReturnType;
        return Dual1<ReturnType>(this->val * o.val, d1 * o.val + this->val * o.d1);
    }
    Dual1 operator/(const Dual1& o) const {
        typedef common_type_t<decltype(this->val), decltype(o.val)> RawType;
        typedef conditional_t<
            std::is_arithmetic<RawType>::value,
            double,
            std::complex<double>
        > ReturnType;
        ReturnType o_val_ = static_cast<ReturnType>(o.val);
        ReturnType o_d1 = static_cast<ReturnType>(o.d1);
        ReturnType val_ = static_cast<ReturnType>(this->val);
        ReturnType d1_ = static_cast<ReturnType>(d1);
        ReturnType denom = o_val_ * o_val_;
        return Dual1<ReturnType>(
            val_ / o_val_,
            (d1_ * o_val_ - val_ * o_d1) / denom
        );
    }

    // Member: Dual1 <op> scalar (Right hand side)
    template<typename U>
    Dual1<common_type_t<T, U>> operator+(U o) const {
        return Dual1<common_type_t<T, U>>(this->val + o, d1);
    }
    template<typename U>
    Dual1<common_type_t<T, U>> operator-(U o) const {
        return Dual1<common_type_t<T, U>>(this->val - o, d1);
    }
    template<typename U>
    Dual1<common_type_t<T, U>> operator*(U o) const {
        return Dual1<common_type_t<T, U>>(this->val * o, d1 * o);
    }
    template<typename U>
    Dual1<conditional_t<
        std::is_arithmetic<common_type_t<T, U>>::value,
        double,
        std::complex<double>
    >> operator/(U o) const {
        typedef common_type_t<T, U> RawType;
        typedef conditional_t<
            std::is_arithmetic<RawType>::value,
            double,
            std::complex<double>
        > ReturnType;
        ReturnType o_ = static_cast<ReturnType>(o);
        ReturnType val_ = static_cast<ReturnType>(this->val);
        ReturnType d1_ = static_cast<ReturnType>(d1);
        return Dual1<ReturnType>(
            val_ / o_,
            d1_ / o_
        );
    }

    //other member functions
    friend std::ostream& operator<<(std::ostream& os, const Dual1& x) {
        os << x.val << " (d1=" << x.d1 << ")";
        return os;
    }
};

// ---  Dual1 <op> scalar --- (left hand side)
template<typename T, typename O>
Dual1<common_type_t<O, T>> operator+(O lhs, const Dual1<T>& rhs) {
    return Dual1<common_type_t<O, T>>(lhs + rhs.val, rhs.d1);
}
template<typename T, typename O>
Dual1<common_type_t<O, T>> operator-(O lhs, const Dual1<T>& rhs) {
    return Dual1<common_type_t<O, T>>(lhs - rhs.val, -rhs.d1);
}
template<typename T, typename O>
Dual1<common_type_t<O, T>> operator*(O lhs, const Dual1<T>& rhs) {
    return Dual1<common_type_t<O, T>>(lhs * rhs.val, lhs * rhs.d1);
}
template<typename T, typename O>
Dual1<conditional_t<
    std::is_arithmetic<common_type_t<O, T>>::value,
    double,
    std::complex<double>
>> operator/(O lhs, const Dual1<T>& rhs) {
    typedef common_type_t<O, T> RawType;
    typedef conditional_t<
        std::is_arithmetic<RawType>::value,
        double,
        std::complex<double>
    > ReturnType;
    ReturnType lhs_ = static_cast<ReturnType>(lhs);
    ReturnType val_ = static_cast<ReturnType>(rhs.val);
    ReturnType d1_ = static_cast<ReturnType>(rhs.d1);
    ReturnType denom = val_ * val_;

    return Dual1<ReturnType>(
        lhs_ / val_,
        -lhs_ * d1_ / denom
    );
}

// Free functions for Dual1
template<typename T>
Dual1<decltype(std::exp(std::declval<T>()))> exp(const Dual1<T>& x) {
    typedef decltype(std::exp(std::declval<T>())) ReturnType;
    ReturnType exp_val = std::exp(x.val);
    ReturnType d1 = x.d1 * exp_val;
    return Dual1<ReturnType>(exp_val, d1);
}
template<typename T>
Dual1<decltype(std::log(std::declval<T>()))> log(const Dual1<T>& x) {
    typedef decltype(std::log(std::declval<T>())) ReturnType;
    ReturnType log_val = std::log(x.val);
    ReturnType d1 = x.d1 / x.val;
    return Dual1<ReturnType>(log_val, d1);
}
template<typename T, typename U>
Dual1<decltype(std::pow(std::declval<T>(), std::declval<U>()))> pow(const Dual1<T>& x, U exponent) {
    typedef decltype(std::pow(std::declval<T>(), std::declval<U>())) ReturnType;
    ReturnType pow_val = std::pow(x.val, exponent);
    ReturnType d1 = exponent * std::pow(x.val, exponent - 1) * x.d1;
    return Dual1<ReturnType>(pow_val, d1);
}
template<typename T, typename U>
Dual1<decltype(std::pow(std::declval<U>(), std::declval<T>()))> pow(U base, const Dual1<T>& x) {
    typedef decltype(std::pow(std::declval<U>(), std::declval<T>())) ReturnType;
    ReturnType pow_val = std::pow(base, x.val);
    ReturnType d1 = x.d1 * std::log(base) * pow_val;
    return Dual1<ReturnType>(pow_val, d1);
}
template<typename T, typename U>
Dual1<decltype(std::pow(std::declval<T>(), std::declval<U>()))> pow(const Dual1<T>& x, const Dual1<U>& y) {
    typedef decltype(std::pow(std::declval<T>(), std::declval<U>())) ReturnType;
    ReturnType pow_val = std::pow(x.val, y.val);
    ReturnType d1 = x.d1 * y.val * std::pow(x.val, y.val - 1) + y.d1 * std::log(x.val) * pow_val;
    return Dual1<ReturnType>(pow_val, d1);
}
template<typename T>
Dual1<decltype(std::sqrt(std::declval<T>()))> sqrt(const Dual1<T>& x) {
    typedef decltype(std::sqrt(std::declval<T>())) ReturnType;
    ReturnType sqrt_val = std::sqrt(x.val);
    ReturnType d1 = x.d1 / (2 * sqrt_val);
    return Dual1<ReturnType>(sqrt_val, d1);
}
template<typename T>
Dual1<decltype(std::sin(std::declval<T>()))> sin(const Dual1<T>& x) {
    typedef decltype(std::sin(std::declval<T>())) ReturnType;
    ReturnType sin_val = std::sin(x.val);
    ReturnType d1 = x.d1 * std::cos(x.val);
    return Dual1<ReturnType>(sin_val, d1);
}
template<typename T>
Dual1<decltype(std::cos(std::declval<T>()))> cos(const Dual1<T>& x) {
    typedef decltype(std::cos(std::declval<T>())) ReturnType;
    ReturnType cos_val = std::cos(x.val);
    ReturnType d1 = -x.d1 * std::sin(x.val);
    return Dual1<ReturnType>(cos_val, d1);
}
template<typename T>
Dual1<decltype(std::tan(std::declval<T>()))> tan(const Dual1<T>& x) {
    typedef decltype(std::tan(std::declval<T>())) ReturnType;
    ReturnType tan_val = std::tan(x.val);
    ReturnType d1 = x.d1 / (std::cos(x.val) * std::cos(x.val));
    return Dual1<ReturnType>(tan_val, d1);
}
