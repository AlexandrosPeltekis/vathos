#include <cmath>
#include <complex>
#include <array>
#include "dualbase_class.h"

// First-order Dual number: inherits from DualBase and adds first derivative
// --- Dual1<T> member operators for scalar types ---
template<typename T>
struct Dual1 : public DualBase<T> {
    T d1;

	// Constructors
    Dual1() : DualBase<T>(), d1(1) {}
    Dual1(const T& v) : DualBase<T>(v), d1(1) {}
    Dual1(const T& v, const T& d) : DualBase<T>(v), d1(d) {}

    // Conversion constructor: Dual1<L> -> Dual1<std::complex<T>> (i.e. Dual1<std::complex<T>> d = Dual1<L>(a,b); )
    template<typename L, typename = std::enable_if_t<!std::is_same_v<T, L>>>
    explicit Dual1(const Dual1<L>& other)
        : DualBase<T>(T(other.val)), d1(T(other.d1)) {}

    // Member: Dual1 <op> Dual1 (Right hand side)
    // Addition
    Dual1 operator+(const Dual1& o) const {
        using ReturnType = std::common_type_t<T, decltype(o.val)>;
        return Dual1<ReturnType>(this->val + o.val, d1 + o.d1);
    }

    // Subtraction
    Dual1 operator-(const Dual1& o) const {
        using ReturnType = std::common_type_t<T, decltype(o.val)>;
        return Dual1<ReturnType>(this->val - o.val, d1 - o.d1);
    }

    // Multiplication
    Dual1 operator*(const Dual1& o) const {
        using ReturnType = std::common_type_t<T, decltype(o.val)>;
        return Dual1<ReturnType>(this->val * o.val, d1 * o.val + this->val * o.d1);
    }

    // Division
    Dual1 operator/(const Dual1& o) const { // More involved due to potential int division truncation
        using RawType = std::common_type_t<decltype(this->val), decltype(o.val)>;
        using ReturnType = std::conditional_t<
            std::is_arithmetic_v<RawType>,
            double, // Use double for arithmetic types (non-complex) to avoid truncation
            std::complex<double> // Use complex<double> for complex types
        >;
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
    // Addition
    template<typename U>
    auto operator+(U o) const {
        using ReturnType = std::common_type_t<T, U>;
        return Dual1<ReturnType>(this->val + o, d1);
    }
    // Subtraction
    template<typename U>
    auto operator-(U o) const {
        using ReturnType = std::common_type_t<T, U>;
        return Dual1<ReturnType>(this->val - o, d1);
    }
    // Multiplication
    template<typename U>
    auto operator*(U o) const {
        using ReturnType = std::common_type_t<T, U>;
        return Dual1<ReturnType>(this->val * o, d1 * o);
    }
    // Division
    template<typename U>
    auto operator/(U o) const { // More involved due to potential int division truncation
        using RawType = std::common_type_t<T, U>;
        using ReturnType = std::conditional_t<
            std::is_arithmetic_v<RawType>,
            double, // Use double for arithmetic types (non-complex) to avoid truncation
            std::complex<double> // Use complex<double> for complex types
        >;

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

//Deduction Guide for mixing types in Dual1
template<typename V, typename D1>
Dual1(V, D1) -> Dual1<std::common_type_t<V, D1>>;

// ---  Dual1 <op> scalar --- (left hand side)
//Addition
template<typename T, typename O>
auto operator+(O lhs, const Dual1<T>& rhs) {
    using ReturnType = std::common_type_t<O, T>;
    return Dual1<ReturnType>(lhs + rhs.val, rhs.d1);
}

// Subtraction
template<typename T, typename O>
auto operator-(O lhs, const Dual1<T>& rhs) {
    using ReturnType = std::common_type_t<O, T>;
    return Dual1<ReturnType>(lhs - rhs.val, -rhs.d1);
}

// Multiplication
template<typename T, typename O>
auto operator*(O lhs, const Dual1<T>& rhs) {
    using ReturnType = std::common_type_t<O, T>;
    return Dual1<ReturnType>(lhs * rhs.val, lhs * rhs.d1);
}

// Division
template<typename T, typename O>
auto operator/(O lhs, const Dual1<T>& rhs) { // More involved due to potential int division truncation
    using RawType = std::common_type_t<O,T>;
    using ReturnType = std::conditional_t<
        std::is_arithmetic_v<RawType>,
        double, // Use double for arithmetic types (non-complex) to avoid truncation
        std::complex<double> // Use complex<double> for complex types
    >;
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
// Exponential for Dual2_3d
template<typename T>
auto exp(const Dual1<T>& x) {
    using ReturnType = Dual1<decltype(std::exp(x.val))>; // Promote to the type of exp(x.val)
    auto exp_val = std::exp(x.val);
    ReturnType d1 = x.d1 * exp_val;
    return ReturnType(exp_val, d1);
}

// Logarithm for Dual1
template<typename T>
auto log(const Dual1<T>& x) {
    using ReturnType = Dual1<decltype(std::log(x.val))>; // Promote to the type of log(x.val)
    auto log_val = std::log(x.val);
    ReturnType d1 = x.d1 / x.val
    return ReturnType(log_val, d1);
}

// pow(x, scalar): x^exponent for Dual1
template<typename T, typename U>
auto pow(const Dual1<T>& x, U exponent) {
    using ReturnType = Dual1<decltype(std::pow(x.val, exponent))>; // Promote to the type of pow(x.val, exponent)
    auto pow_val = std::pow(x.val, exponent);
    ReturnType d1 = exponent * pow(x.val, exponent - 1) * x.d1
    return ReturnType(pow_val, d1);
}

// pow(scalar, x): exponent^x for Dual1
template<typename T, typename U>
auto pow(U base, const Dual1<T>& x) {
    using ReturnType = Dual1<decltype(std::pow(base, x.val))>; // Promote to the type of pow(base, x.val)
    auto pow_val = std::pow(base, x.val);
    ReturnType d1 = x.d1 * std::log(base) * pow_val
    return ReturnType(pow_val, d1);
}

// pow(x, y): x^y for Dual1
template<typename T, typename U>
auto pow(const Dual1<T>& x, const Dual1<U>& y) {
    using ReturnType = Dual1<decltype(std::pow(x.val, y.val))>; // Promote to the type of pow(x.val, y.val)
    auto pow_val = std::pow(x.val, y.val);
    ReturnType d1 = x.d1 * y.val * std::pow(x.val, y.val - 1) +
        y.d1 * std::log(x.val) * pow_val
    return ReturnType(pow_val, d1);
}

// Square root for Dual1
template<typename T>
auto sqrt(const Dual1<T>& x) {
    using ReturnType = Dual1<decltype(std::sqrt(x.val))>; // Promote to the type of sqrt(x.val)
    auto sqrt_val = std::sqrt(x.val);
    ReturnType d1 = x.d1 / (2 * sqrt_val);
    return ReturnType(sqrt_val, d1);
}
