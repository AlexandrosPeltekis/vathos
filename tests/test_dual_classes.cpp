/*
 *  Minimal unit tests for DualBase, Dual1, and Dual2_3d classes.
 *  No external dependencies; uses standard C++ and assert.
 */

#include <iostream>
#include <cmath>
#include <complex>
#include <sstream>
#include <cassert>
#include "../dualbase_class.h"
#include "../dual1_class.h"
#include "../dual2_3d_class.h"

#define EXPECT_EQ(a, b) do { \
    if (!((a) == (b))) { \
        std::cerr << "EXPECT_EQ failed: " << #a << " == " << #b << " (" << (a) << " != " << (b) << ") at " << __FILE__ << ":" << __LINE__ << std::endl; \
        ++failures; \
    } else { ++successes; } \
} while(0)

#define EXPECT_NEAR(a, b, tol) do { \
    if (!(std::abs((a)-(b)) <= (tol))) { \
        std::cerr << "EXPECT_NEAR failed: " << #a << " ~= " << #b << " (" << (a) << " != " << (b) << ") at " << __FILE__ << ":" << __LINE__ << std::endl; \
        ++failures; \
    } else { ++successes; } \
} while(0)

int failures = 0, successes = 0;

void test_dualbase() {
    DualBase<double> d1;
    EXPECT_EQ(d1.val, 0.0);

    DualBase<int> d2(42);
    EXPECT_EQ(d2.val, 42);

    std::ostringstream oss;
    oss << d2;
    EXPECT_EQ(oss.str(), "42");
}

void test_dual1_constructors() {
    Dual1<double> d1;
    EXPECT_EQ(d1.val, 0.0);
    EXPECT_EQ(d1.d1, 1.0);

    Dual1<double> d2(5.0);
    EXPECT_EQ(d2.val, 5.0);
    EXPECT_EQ(d2.d1, 1.0);

    Dual1<double> d3(2.0, 3.0);
    EXPECT_EQ(d3.val, 2.0);
    EXPECT_EQ(d3.d1, 3.0);

    Dual1<int> di(7, 8);
    Dual1<std::complex<double>> dc(di);
    EXPECT_EQ(dc.val, std::complex<double>(7, 0));
    EXPECT_EQ(dc.d1, std::complex<double>(8, 0));
}

void test_dual1_operators() {
    Dual1<double> a(2.0, 3.0);
    Dual1<double> b(5.0, 7.0);

    auto sum = a + b;
    EXPECT_NEAR(sum.val, 7.0, 1e-12);
    EXPECT_NEAR(sum.d1, 10.0, 1e-12);

    auto diff = a - b;
    EXPECT_NEAR(diff.val, -3.0, 1e-12);
    EXPECT_NEAR(diff.d1, -4.0, 1e-12);

    auto prod = a * b;
    EXPECT_NEAR(prod.val, 10.0, 1e-12);
    EXPECT_NEAR(prod.d1, 2.0*7.0 + 3.0*5.0, 1e-12);

    auto quot = a / b;
    EXPECT_NEAR(quot.val, 2.0/5.0, 1e-12);
    EXPECT_NEAR(quot.d1, (3.0*5.0 - 2.0*7.0)/(5.0*5.0), 1e-12);
}

void test_dual1_scalar_ops() {
    Dual1<double> a(2.0, 3.0);

    auto sum = a + 4.0;
    EXPECT_NEAR(sum.val, 6.0, 1e-12);
    EXPECT_NEAR(sum.d1, 3.0, 1e-12);

    auto diff = a - 1.0;
    EXPECT_NEAR(diff.val, 1.0, 1e-12);
    EXPECT_NEAR(diff.d1, 3.0, 1e-12);

    auto prod = a * 2.0;
    EXPECT_NEAR(prod.val, 4.0, 1e-12);
    EXPECT_NEAR(prod.d1, 6.0, 1e-12);

    auto quot = a / 2.0;
    EXPECT_NEAR(quot.val, 1.0, 1e-12);
    EXPECT_NEAR(quot.d1, 1.5, 1e-12);
}

void test_dual1_left_scalar_ops() {
    Dual1<double> a(2.0, 3.0);

    auto sum = 4.0 + a;
    EXPECT_NEAR(sum.val, 6.0, 1e-12);
    EXPECT_NEAR(sum.d1, 3.0, 1e-12);

    auto diff = 5.0 - a;
    EXPECT_NEAR(diff.val, 3.0, 1e-12);
    EXPECT_NEAR(diff.d1, -3.0, 1e-12);

    auto prod = 2.0 * a;
    EXPECT_NEAR(prod.val, 4.0, 1e-12);
    EXPECT_NEAR(prod.d1, 6.0, 1e-12);

    auto quot = 6.0 / a;
    EXPECT_NEAR(quot.val, 3.0, 1e-12);
    EXPECT_NEAR(quot.d1, -4.5, 1e-12);
}

void test_dual1_ostream() {
    Dual1<double> a(2.0, 3.0);
    std::ostringstream oss;
    oss << a;
    EXPECT_EQ(oss.str(), "2 (d1=3)");
}

void test_dual1_exp() {
    Dual1<double> a(1.0, 2.0);
    auto e = exp(a);
    EXPECT_NEAR(e.val, std::exp(1.0), 1e-12);
    EXPECT_NEAR(e.d1, 2.0 * std::exp(1.0), 1e-12);
}

void test_dual2_3d_constructors() {
    Dual2_3d<double> d0;
    EXPECT_EQ(d0.val, 0.0);
    EXPECT_EQ(d0.dx, 0.0);
    EXPECT_EQ(d0.dy, 0.0);
    EXPECT_EQ(d0.dz, 0.0);
    EXPECT_EQ(d0.dxx, 0.0);
    EXPECT_EQ(d0.dxy, 0.0);
    EXPECT_EQ(d0.dxz, 0.0);
    EXPECT_EQ(d0.dyy, 0.0);
    EXPECT_EQ(d0.dyz, 0.0);
    EXPECT_EQ(d0.dzz, 0.0);

    Dual2_3d<double> d1(5.0);
    EXPECT_EQ(d1.val, 5.0);
    EXPECT_EQ(d1.dx, 0.0);

    Dual2_3d<double> d2(1.0, 2.0, 3.0, 4.0);
    EXPECT_EQ(d2.val, 1.0);
    EXPECT_EQ(d2.dx, 2.0);
    EXPECT_EQ(d2.dy, 3.0);
    EXPECT_EQ(d2.dz, 4.0);
    EXPECT_EQ(d2.dxx, 0.0);

    Dual2_3d<double> d3(1,2,3,4,5,6,7,8,9,10);
    EXPECT_EQ(d3.val, 1.0);
    EXPECT_EQ(d3.dx, 2.0);
    EXPECT_EQ(d3.dy, 3.0);
    EXPECT_EQ(d3.dz, 4.0);
    EXPECT_EQ(d3.dxx, 5.0);
    EXPECT_EQ(d3.dxy, 6.0);
    EXPECT_EQ(d3.dxz, 7.0);
    EXPECT_EQ(d3.dyy, 8.0);
    EXPECT_EQ(d3.dyz, 9.0);
    EXPECT_EQ(d3.dzz, 10.0);
}

void test_dual2_3d_operators() {
    Dual2_3d<double> a(1,2,3,4,5,6,7,8,9,10);
    Dual2_3d<double> b(10,9,8,7,6,5,4,3,2,1);

    auto sum = a + b;
    EXPECT_NEAR(sum.val, 11.0, 1e-12);
    EXPECT_NEAR(sum.dx, 11.0, 1e-12);
    EXPECT_NEAR(sum.dxy, 11.0, 1e-12);

    auto diff = a - b;
    EXPECT_NEAR(diff.val, -9.0, 1e-12);
    EXPECT_NEAR(diff.dx, -7.0, 1e-12);
    EXPECT_NEAR(diff.dxy, 1.0, 1e-12);

    auto prod = a * b;
    EXPECT_NEAR(prod.val, 10.0, 1e-12);
    EXPECT_NEAR(prod.dx, 2.0*10.0 + 1.0*9.0, 1e-12);
    EXPECT_NEAR(prod.dxy, 6.0*10.0 + 2.0*8.0 + 3.0*9.0 + 1.0*5.0, 1e-12);
}

void test_dual2_3d_scalar_ops() {
    Dual2_3d<double> a(1,2,3,4,5,6,7,8,9,10);

    auto sum = a + 2.0;
    EXPECT_NEAR(sum.val, 3.0, 1e-12);
    EXPECT_NEAR(sum.dx, 2.0, 1e-12);

    auto diff = a - 1.0;
    EXPECT_NEAR(diff.val, 0.0, 1e-12);
    EXPECT_NEAR(diff.dx, 2.0, 1e-12);

    auto prod = a * 2.0;
    EXPECT_NEAR(prod.val, 2.0, 1e-12);
    EXPECT_NEAR(prod.dx, 4.0, 1e-12);

    auto quot = a / 2.0;
    EXPECT_NEAR(quot.val, 0.5, 1e-12);
    EXPECT_NEAR(quot.dx, 1.0, 1e-12);
}

void test_dual2_3d_left_scalar_ops() {
    Dual2_3d<double> a(1,2,3,4,5,6,7,8,9,10);

    auto sum = 2.0 + a;
    EXPECT_NEAR(sum.val, 3.0, 1e-12);
    EXPECT_NEAR(sum.dx, 2.0, 1e-12);

    auto diff = 5.0 - a;
    EXPECT_NEAR(diff.val, 4.0, 1e-12);
    EXPECT_NEAR(diff.dx, -2.0, 1e-12);

    auto prod = 3.0 * a;
    EXPECT_NEAR(prod.val, 3.0, 1e-12);
    EXPECT_NEAR(prod.dx, 6.0, 1e-12);

    auto quot = 6.0 / a;
    EXPECT_NEAR(quot.val, 6.0, 1e-12);
    EXPECT_NEAR(quot.dx, -12.0, 1e-12);
}

void test_dual2_3d_ostream() {
    Dual2_3d<double> a(1,2,3,4,5,6,7,8,9,10);
    std::ostringstream oss;
    oss << a;
    std::string s = oss.str();
    EXPECT_EQ(s.substr(0, 2), "1 ");
    EXPECT_NEAR(a.dx, 2.0, 1e-12); // just to use a.dx
}

void test_dual2_3d_functions() {
    Dual2_3d<double> a(4,2,0,0,0,0,0,0,0,0);

    auto p = pow(a, 2.0);
    EXPECT_NEAR(p.val, 16.0, 1e-12);
    EXPECT_NEAR(p.dx, 2.0 * 2.0 * 4.0, 1e-12);

    auto e = exp(a);
    EXPECT_NEAR(e.val, std::exp(4.0), 1e-12);
    EXPECT_NEAR(e.dx, 2.0 * std::exp(4.0), 1e-12);

    auto s = sqrt(a);
    EXPECT_NEAR(s.val, 2.0, 1e-12);
    EXPECT_NEAR(s.dx, 0.5 / 2.0 * 2.0, 1e-12);
}

int main() {
    test_dualbase();
    test_dual1_constructors();
    test_dual1_operators();
    test_dual1_scalar_ops();
    test_dual1_left_scalar_ops();
    test_dual1_ostream();
    test_dual1_exp();
    test_dual2_3d_constructors();
    test_dual2_3d_operators();
    test_dual2_3d_scalar_ops();
    test_dual2_3d_left_scalar_ops();
    test_dual2_3d_ostream();
    test_dual2_3d_functions();

    std::cout << "Unit tests: " << successes << " passed, " << failures << " failed." << std::endl;
    return failures == 0 ? 0 : 1;
}