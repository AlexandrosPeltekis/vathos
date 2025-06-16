#include <cmath>
#include <complex>
#include <array>
#include "dualbase_class.h"


// Second-order Dual number: inherits from DualBase and adds first and second derivatives
template<typename T>
struct Dual2 : public DualBase<T> {
    T d1; // first derivative
    T d2; // second derivative

    Dual2() : DualBase<T>(), d1(1), d2(1) {}
    Dual2(const T& v) : DualBase<T>(v), d1(1), d2(1) {}
    Dual2(const T& v, const T& d1_) : DualBase<T>(v), d1(d1_), d2(1) {}
    Dual2(const T& v, const T& d1_, const T& d2_) : DualBase<T>(v), d1(d1_), d2(d2_) {}

	// Conversion constructor: Dual2<L> -> Dual2<std::complex<T>> (i.e. Dual2<std::complex<T>> d = Dual2<L>(a,b,c); )
	template<typename L>
    explicit Dual2(const Dual2<L>& other)
		requires (!std::is_same_v<T, L>)
	: DualBase<T>(T(other.val)), d1(T(other.d1)), d2(T(other.d2)) {}


    // Member: Dual2 <op> Dual2
    // Addition
    Dual2 operator+(const Dual2& o) const {
        using ReturnType = std::common_type_t<T, decltype(o.val)>;
        return Dual2<ReturnType>(this->val + o.val, d1 + o.d1, d2 + o.d2);
    }
    // Subtraction
    Dual2 operator-(const Dual2& o) const {
        using ReturnType = std::common_type_t<T, decltype(o.val)>;
        return Dual2<ReturnType>(this->val - o.val, d1 - o.d1, d2 - o.d2);
    }
    // Multiplication
    Dual2 operator*(const Dual2& o) const {
        using ReturnType = std::common_type_t<T, decltype(o.val)>;
        ReturnType new_val = this->val * o.val;
        ReturnType new_d1 = d1 * o.val + this->val * o.d1;
        ReturnType new_d2 = d2 * o.val + 2.0 * d1 * o.d1 + this->val * o.d2;
        return Dual2<ReturnType>(new_val, new_d1, new_d2 );

    }
    // Division
    Dual2 operator/(const Dual2& o) const { // More involved due to potential int division truncation
        using RawType = std::common_type_t<decltype(this->val), decltype(o.val)>;
        using ReturnType = std::conditional_t<
            std::is_arithmetic_v<RawType>, 
			double, // Use double for arithmetic types (non-complex) to avoid truncation
			std::complex<double> // Use complex<double> for complex types
        >;
        ReturnType o_val_ = static_cast<ReturnType>(o.val);
        ReturnType o_d1 = static_cast<ReturnType>(o.d1);
        ReturnType o_d2 = static_cast<ReturnType>(o.d2);
        ReturnType val_ = static_cast<ReturnType>(this->val);
        ReturnType d1_ = static_cast<ReturnType>(d1);
        ReturnType d2_ = static_cast<ReturnType>(d2);

        ReturnType denom = o_val_ * o_val_;
        ReturnType new_val = val_ / o_val_;
        ReturnType new_d1 = (d1_ * o_val_ - val_ * o_d1) / denom;
        ReturnType new_d2 = (d2_ * o_val_ * o_val_ - o_val_ * val_ * o_d2 - 2.0 * o_val_ * d1_ * o_d1 - 2.0 * val_ * o_d1 * o_d1) / (denom * o_val_);

        return Dual2<ReturnType>(
            new_val,
            new_d1,
            new_d2
        );
    }

    // Member: Dual2 <op> scalar (Right hand side)
    template<typename U>
    auto operator+(U o) const {
        using ReturnType = std::common_type_t<T, U>;
        return Dual2<ReturnType>(this->val + o, d1, d2);
    }

    template<typename U>
    auto operator-(U o) const {
        using ReturnType = std::common_type_t<T, U>;
        return Dual2<ReturnType>(this->val - o, d1, d2);
    }

    template<typename U>
    auto operator*(U o) const {
        using ReturnType = std::common_type_t<T, U>;
        return Dual2<ReturnType>(this->val * o, d1 * o, d2 * o);
    }

    template<typename U>
    auto operator/(U o) const {
        using RawType = std::common_type_t<T, U>;
        using ReturnType = std::conditional_t<
            std::is_arithmetic_v<RawType>,
            double, // Use double for arithmetic types (non-complex) to avoid truncation
            std::complex<double> // Use complex<double> for complex types
        >;

        ReturnType o_ = static_cast<ReturnType>(o);
        ReturnType val_ = static_cast<ReturnType>(this->val);
        ReturnType d1_ = static_cast<ReturnType>(d1);
		ReturnType d2_ = static_cast<ReturnType>(d2);

        return Dual2<ReturnType>(
            val_ / o_,
            d1_ / o_,
            d2_ / o_
        );
    }

    friend std::ostream& operator<<(std::ostream& os, const Dual2& x) {
        os << x.val << " (d1=" << x.d1 << ", d2=" << x.d2 << ")";
        return os;
    }
};

//Deduction Guide for mixing types in Dual2
template<typename V, typename D1, typename D2>
Dual2(V, D1, D2) -> Dual2<std::common_type_t<V, D1, D2>>;

// ---  Dual2 <op> scalar --- (left hand side)
//Addition
template<typename T, typename O>
auto operator+(O lhs, const Dual2<T>& rhs) {
    using ReturnType = std::common_type_t<O, T>;
	return Dual2<ReturnType>(lhs + rhs.val, rhs.d1, rhs.d2);
}

// Subtraction
template<typename T, typename O>
auto operator-(O lhs, const Dual2<T>& rhs) {
    using ReturnType = std::common_type_t<O, T>;
	return Dual2<ReturnType>(lhs - rhs.val, -rhs.d1, -rhs.d2);
}

// Multiplication
template<typename T, typename O>
auto operator*(O lhs, const Dual2<T>& rhs) {
    using ReturnType = std::common_type_t<O, T>;
	return Dual2<ReturnType>(lhs * rhs.val, lhs * rhs.d1, lhs * rhs.d2);
}

// Division
template<typename T, typename O>
auto operator/(O lhs, const Dual2<T>& rhs) {
    using RawType = std::common_type_t<O, T>;
    using ReturnType = std::conditional_t<
        std::is_arithmetic_v<RawType>,
        double, // Use double for arithmetic types (non-complex) to avoid truncation
        std::complex<double> // Use complex<double> for complex types
    >;
    ReturnType lhs_ = static_cast<ReturnType>(lhs);
    ReturnType val_ = static_cast<ReturnType>(rhs.val);
    ReturnType d1_ = static_cast<ReturnType>(rhs.d1);
	ReturnType d2_ = static_cast<ReturnType>(rhs.d2);
    ReturnType denom = val_ * val_;

    return Dual2<ReturnType>(
        lhs_ / val_,
        -lhs_ * d1_ / denom,
		(2.0 * lhs_ * d1_ * d1_ - lhs_ * val_ * d2_) / (denom * val_)
    );
}

// Free functions for Dual2
template<typename T>
Dual2<T> exp(const Dual2<T>& x) {
    T exp_val = std::exp(x.val);
    return { exp_val, x.d1 * exp_val, (x.d1 * x.d1 + x.d2) * exp_val };
}

template<typename T>
Dual2<std::complex<T>> exp(const Dual2<std::complex<T>>& x) {
    std::complex<T> exp_val = std::exp(x.val);
    return { exp_val, x.d1 * exp_val, (x.d1 * x.d1 + x.d2) * exp_val };
}

template <typename T>
Dual2<T> pow(const Dual2<T>& base, int exponent) {
    if (exponent == 0) {
        return Dual2<T>(1.0); // Any number to the power of 0 is 1
    }

    Dual2<T> result = base;
    for (int i = 1; i < exponent; ++i) {
        result = result * base;
    }
    return result;
}

template <typename T>
Dual2<T> sqrt(const Dual2<T>& x) {
    if (x.val < 0) {
        throw std::domain_error("Cannot compute square root of negative number");
    }
    T sqrt_val = std::sqrt(x.val);
    return { sqrt_val, x.d1 / (2.0 * sqrt_val), x.d2 / (4.0 * sqrt_val * sqrt_val) };
}
