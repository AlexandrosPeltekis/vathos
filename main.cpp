#include <iostream>
#include <cmath>
#include <complex>
#include <array>
#include <vector>
#include "dual1_class.h"
#include "dual2_3d_class.h"
#include <cassert>

/*
 * Test function for Dual1 operators
 * This function tests the addition, subtraction, multiplication, and division of Dual1 numbers
 * with both scalars and other Dual1 numbers.
 */
static void test_dual1_operators(bool verbose = false) {
    double s1 = 2.0;
    std::complex<double> c1(1.0, 2.0);
    Dual1<double> a(2.0, 1.0); // a = 2, da/dx = 1
    Dual1<double> b(3.0, 0.5); // b = 3, db/dx = 0.5

    if (verbose) {
        std::cout << "\nTesting left and right hand side operators with duals:\n";
        std::cout << "a + b: " << (a + b) << "\n";
        std::cout << "a - b: " << (a - b) << "\n";
        std::cout << "a * b: " << (a * b) << "\n";
        std::cout << "a / b: " << (a / b) << "\n";
        std::cout << "\nTesting right hand side operators with dual and scalars:\n";
        std::cout << "a + 2: " << (a + s1) << "\n";
        std::cout << "a - 2: " << (a - s1) << "\n";
        std::cout << "a * 2: " << (a * s1) << "\n";
        std::cout << "a / 2: " << (a / s1) << "\n";
        std::cout << "\nTesting left hand side operators with dual and scalars:\n";
        std::cout << "2 + a: " << (s1 + a) << "\n";
        std::cout << "2 - a: " << (s1 - a) << "\n";
        std::cout << "2 * a: " << (s1 * a) << "\n";
        std::cout << "2 / a: " << (s1 / a) << "\n";
        std::cout << "\nTesting right hand side operators with dual and complex scalars:\n";
        std::cout << "a + c1: " << (a + c1) << "\n";
        std::cout << "a - c1: " << (a - c1) << "\n";
        std::cout << "a * c1: " << (a * c1) << "\n";
        std::cout << "a / c1: " << (a / c1) << "\n";
        std::cout << "\nTesting left hand side operators with dual and complex scalars:\n";
        std::cout << "c1 + a: " << (c1 + a) << "\n";
        std::cout << "c1 - a: " << (c1 - a) << "\n";
        std::cout << "c1 * a: " << (c1 * a) << "\n";
        std::cout << "c1 / a: " << (c1 / a) << "\n";
    }
    else
    {
        // If verbose is false, we can skip printing the results and just perform the operations
        Dual1<double> result1 = a + b;
        Dual1<double> result2 = a - b;
        Dual1<double> result3 = a * b;
        Dual1<double> result4 = a / b;
        Dual1<double> result5 = a + s1;
        Dual1<double> result6 = a - s1;
        Dual1<double> result7 = a * s1;
        Dual1<double> result8 = a / s1;
        Dual1<double> result9 = s1 + a;
        Dual1<double> result10 = s1 - a;
        Dual1<double> result11 = s1 * a;
        Dual1<double> result12 = s1 / a;
        Dual1<std::complex<double>> result13 = a + c1;
        Dual1<std::complex<double>> result14 = a - c1;
        Dual1<std::complex<double>> result15 = a * c1;
        Dual1<std::complex<double>> result16 = a / c1;
        Dual1<std::complex<double>> result17 = c1 + a;
        Dual1<std::complex<double>> result18 = c1 - a;
        Dual1<std::complex<double>> result19 = c1 * a;
        Dual1<std::complex<double>> result20 = c1 / a;

        std::cout << "> Dual1 non-complex functioning...\n";
    }
}

/*
 * Test function for Dual1 complex operators
 * This function tests the addition, subtraction, multiplication, and division of Dual1 complex numbers
 * with both scalars and other Dual1 complex numbers.
 */
static void test_dual1_complex_operators(bool verbose = false) {
    // Initailize complex numbers and the Dual1 complex duals
    std::complex<double> a1(1.0, 1);
    std::complex<double> b1(3.0, 4.0);
    std::complex<double> c1(1, 2.0);
    Dual1<std::complex<double>> a(a1, 1.0); // a = 1 + 2i, da/dx = 1
    Dual1<std::complex<double>> b(b1, 0.5); // b = 3 + 4i, db/dx = 0.5
    if (verbose) {
        // Print the inital values
        std::cout << "Initial values:\n";
        std::cout << "a: " << a << "\n"; // Should print (1 + 2i, 1)
        std::cout << "b: " << b << "\n"; // Should print (3 + 4i, 0.5)
        std::cout << "c1: " << c1 << "\n\n"; // Should print (1 + 2i)


        std::cout << "\nTesting left and right hand side operators with complex duals:\n";
        std::cout << "a + b: " << (a + b) << "\n";
        std::cout << "a - b: " << (a - b) << "\n";
        std::cout << "a * b: " << (a * b) << "\n";
        std::cout << "a / b: " << (a / b) << "\n";
        std::cout << "\nTesting right hand side operators with complex dual and scalars:\n";
        std::cout << "a + 2: " << (a + 2.0) << "\n"; // Should print (4, 1, 0)
        std::cout << "a - 2: " << (a - 2.0) << "\n"; // Should print (0, 1, 0)
        std::cout << "a * 2: " << (a * 2.0) << "\n"; // Should print (4, 2, 0)
        std::cout << "a / 2: " << (a / 2.0) << "\n"; // Should print (1, 0.5, 0)
        std::cout << "\nTesting left hand side operators with complex dual and scalars:\n";
        std::cout << "2 + a: " << (2.0 + a) << "\n"; // Should print (4, 1, 0)
        std::cout << "2 - a: " << (2.0 - a) << "\n"; // Should print (0, -1, 0)
        std::cout << "2 * a: " << (2.0 * a) << "\n"; // Should print (4, 2, 0)
        std::cout << "2 / a: " << (2.0 / a) << "\n"; // Should print (1, -0.5, 0)
        std::cout << "\nTesting right hand side operators with complex dual and complex scalars:\n";
        std::cout << "a + c1: " << (a + c1) << "\n"; // Should print (3, 1, 2)
        std::cout << "a - c1: " << (a - c1) << "\n"; // Should print (1, 1, -2)
        std::cout << "a * c1: " << (a * c1) << "\n"; // Should print (2 + 2i, 1 + 2i, 0)
        std::cout << "a / c1: " << (a / c1) << "\n"; // Should print (0.8 - 0.4i, 0.4 + 0.2i, 0)
        std::cout << "\nTesting left hand side operators with complex dual and complex scalars:\n";
        std::cout << "c1 + a: " << (c1 + a) << "\n"; // Should print (3, 1, 2)
        std::cout << "c1 - a: " << (c1 - a) << "\n"; // Should print (-1, -1, -2)
        std::cout << "c1 * a: " << (c1 * a) << "\n"; // Should print (2 + 2i, 1 + 2i, 0)
        std::cout << "c1 / a: " << (c1 / a) << "\n"; // Should print (0.8 - 0.4i, 0.4 + 0.2i, 0)
    }
    else
    {
        // If verbose is false, we can skip printing the results and just perform the operations
        Dual1<std::complex<double>> result1 = a + b;
        Dual1<std::complex<double>> result2 = a - b;
        Dual1<std::complex<double>> result3 = a * b;
        Dual1<std::complex<double>> result4 = a / b;
        Dual1<std::complex<double>> result5 = a + 2.0;
        Dual1<std::complex<double>> result6 = a - 2.0;
        Dual1<std::complex<double>> result7 = a * 2.0;
        Dual1<std::complex<double>> result8 = a / 2.0;
        Dual1<std::complex<double>> result9 = 2.0 + a;
        Dual1<std::complex<double>> result10 = 2.0 - a;
        Dual1<std::complex<double>> result11 = 2.0 * a;
        Dual1<std::complex<double>> result12 = 2.0 / a;
        Dual1<std::complex<double>> result13 = a + c1;
        Dual1<std::complex<double>> result14 = a - c1;
        Dual1<std::complex<double>> result15 = a * c1;
        Dual1<std::complex<double>> result16 = a / c1;
        Dual1<std::complex<double>> result17 = c1 + a;
        Dual1<std::complex<double>> result18 = c1 - a;
        Dual1<std::complex<double>> result19 = c1 * a;
        Dual1<std::complex<double>> result20 = c1 / a;

        std::cout << "> Dual1 complex functioning...\n";
    }
}

/*
 * Test function for Dual2_3d operators.
 * This function tests the addition, subtraction, multiplication, and division
 * of Dual2_3d objects, both with each other and with scalars (both real and complex).
 * It can be run in verbose mode to print the results of each operation.
 */
static void test_dual2_3d_operators(bool verbose = false) {
    Dual2_3d<double> a(2.0, 1.0, 0.5, 0); // a = 2, da/dx = 1, d^2a/dx^2 = 0.5
    Dual2_3d<double> b(3.0, 0.5, 0.25, 1); // b = 3, db/dx = 0.5, d^2b/dx^2 = 0.25
    std::complex<double> c1(1.0, 2.0);
    if (verbose) {
        std::cout << "\nTesting left and right hand side operators with 3D duals:\n";
        std::cout << "a + b: " << (a + b) << "\n";
        std::cout << "a - b: " << (a - b) << "\n";
        std::cout << "a * b: " << (a * b) << "\n";
        std::cout << "a / b: " << (a / b) << "\n";
        std::cout << "\nTesting right hand side operators with 3D dual and scalars:\n";
        std::cout << "a + 2: " << (a + 2.0) << "\n"; // Should print (4, 1, 0.5)
        std::cout << "a - 2: " << (a - 2.0) << "\n"; // Should print (0, 1, 0.5)
        std::cout << "a * 2: " << (a * 2.0) << "\n"; // Should print (4, 2, 1)
        std::cout << "a / 2: " << (a / 2.0) << "\n"; // Should print (1, 0.5, 0.25)
        std::cout << "\nTesting left hand side operators with 3D dual and scalars:\n";
        std::cout << "2 + a: " << (2.0 + a) << "\n"; // Should print (4, 1, 0.5)
        std::cout << "2 - a: " << (2.0 - a) << "\n"; // Should print (0, -1, -0.5)
        std::cout << "2 * a: " << (2.0 * a) << "\n"; // Should print (4, 2, 1)
        std::cout << "2 / a: " << (2.0 / a) << "\n"; // Should print (1, -0.5, -0.25)
        std::cout << "\nTesting right hand side operators with 3D dual and complex scalars:\n";
        std::cout << "a + c1: " << (a + c1) << "\n"; // Should print (3, 1, 2)
        std::cout << "a - c1: " << (a - c1) << "\n"; // Should print (1, 1, -2)
        std::cout << "a * c1: " << (a * c1) << "\n"; // Should print (2 + 2i, 1 + 2i, 0.5)
        std::cout << "a / c1: " << (a / c1) << "\n"; // Should print (0.8 - 0.4i, 0.4 + 0.2i, -0.25)
        std::cout << "\nTesting left hand side operators with 3D dual and complex scalars:\n";
        std::cout << "c1 + a: " << (c1 + a) << "\n"; // Should print (3, 1, 2)
        std::cout << "c1 - a: " << (c1 - a) << "\n"; // Should print (-1, -1, -2)
        std::cout << "c1 * a: " << (c1 * a) << "\n"; // Should print (2 + 2i, 1 + 2i, 0.5)
        std::cout << "c1 / a: " << (c1 / a) << "\n"; // Should print (0.8 - 0.4i, 0.4 + 0.2i, -0.25);
    }
    else
    {
        // If verbose is false, we can skip printing the results and just perform the operations
        Dual2_3d<double> result1 = a + b;
        Dual2_3d<double> result2 = a - b;
        Dual2_3d<double> result3 = a * b;
        Dual2_3d<double> result4 = a / b;
        Dual2_3d<double> result5 = a + 2.0;
        Dual2_3d<double> result6 = a - 2.0;
        Dual2_3d<double> result7 = a * 2.0;
        Dual2_3d<double> result8 = a / 2.0;
        Dual2_3d<double> result9 = 2.0 + a;
        Dual2_3d<double> result10 = 2.0 - a;
        Dual2_3d<double> result11 = 2.0 * a;
        Dual2_3d<double> result12 = 2.0 / a;
        Dual2_3d<std::complex<double>> result13 = a + c1;
        Dual2_3d<std::complex<double>> result14 = a - c1;
        Dual2_3d<std::complex<double>> result15 = a * c1;
        Dual2_3d<std::complex<double>> result16 = a / c1;
        Dual2_3d<std::complex<double>> result17 = c1 + a;
        Dual2_3d<std::complex<double>> result18 = c1 - a;
        Dual2_3d<std::complex<double>> result19 = c1 * a;
        Dual2_3d<std::complex<double>> result20 = c1 / a;

        std::cout << "> Dual2_3d non-complex functioning...\n";
    }
}

/*
 * Test function for Dual2_3d operators with complex numbers.
 * This function tests the addition, subtraction, multiplication, and division
 * of Dual2_3d objects with both scalars and complex numbers.
 * It can be run with verbose output to see the results of each operation.
 */
static void test_dual2_3d_operators_complex(bool verbose = false) {
    Dual2_3d<std::complex<double>> a(2.0, 1.0, 0.5, 0); // a = 2, da/dx = 1, d^2a/dx^2 = 0.5
    Dual2_3d<std::complex<double>> b(3.0, 0.5, 0.25, 1); // b = 3, db/dx = 0.5, d^2b/dx^2 = 0.25
    std::complex<double> c1(1.0, 2.0);
    if (verbose) {
        std::cout << "\nTesting left and right hand side operators with complex 3D duals:\n";
        std::cout << "a + b: " << (a + b) << "\n";
        std::cout << "a - b: " << (a - b) << "\n";
        std::cout << "a * b: " << (a * b) << "\n";
        std::cout << "a / b: " << (a / b) << "\n";
        std::cout << "\nTesting right hand side operators with complex 3D dual and scalars:\n";
        std::cout << "a + 2: " << (a + 2.0) << "\n"; // Should print (4, 1, 0.5)
        std::cout << "a - 2: " << (a - 2.0) << "\n"; // Should print (0, 1, 0.5)
        std::cout << "a * 2: " << (a * 2.0) << "\n"; // Should print (4, 2, 1)
        std::cout << "a / 2: " << (a / 2.0) << "\n"; // Should print (1, 0.5, 0.25)
        std::cout << "\nTesting left hand side operators with complex 3D dual and scalars:\n";
        std::cout << "2 + a: " << (2.0 + a) << "\n"; // Should print (4, 1, 0.5)
        std::cout << "2 - a: " << (2.0 - a) << "\n"; // Should print (0, -1, -0.5)
        std::cout << "2 * a: " << (2.0 * a) << "\n"; // Should print (4, 2, 1)
        std::cout << "2 / a: " << (2.0 / a) << "\n"; // Should print (1, -0.5, -0.25)
        std::cout << "\nTesting right hand side operators with complex 3D dual and complex scalars:\n";
        std::cout << "a + c1: " << (a + c1) << "\n"; // Should print (3, 1, 2)
        std::cout << "a - c1: " << (a - c1) << "\n"; // Should print (1, 1, -2)
        std::cout << "a * c1: " << (a * c1) << "\n"; // Should print (2 + 2i, 1 + 2i, 0.5)
        std::cout << "a / c1: " << (a / c1) << "\n"; // Should print (0.8 - 0.4i, 0.4 + 0.2i, -0.25)
        std::cout << "\nTesting left hand side operators with complex 3D dual and complex scalars:\n";
        std::cout << "c1 + a: " << (c1 + a) << "\n"; // Should print (3, 1, 2)
        std::cout << "c1 - a: " << (c1 - a) << "\n"; // Should print (-1, -1, -2)
        std::cout << "c1 * a: " << (c1 * a) << "\n"; // Should print (2 + 2i, 1 + 2i, 0.5)
        std::cout << "c1 / a: " << (c1 / a) << "\n"; // Should print (0.8 - 0.4i, 0.4 + 0.2i, -0.25)
    }
    else
    {
        // If verbose is false, we can skip printing the results and just perform the operations
        Dual2_3d<std::complex<double>> result1 = a + b;
        Dual2_3d<std::complex<double>> result2 = a - b;
        Dual2_3d<std::complex<double>> result3 = a * b;
        Dual2_3d<std::complex<double>> result4 = a / b;
        Dual2_3d<std::complex<double>> result5 = a + 2.0;
        Dual2_3d<std::complex<double>> result6 = a - 2.0;
        Dual2_3d<std::complex<double>> result7 = a * 2.0;
        Dual2_3d<std::complex<double>> result8 = a / 2.0;
        Dual2_3d<std::complex<double>> result9 = 2.0 + a;
        Dual2_3d<std::complex<double>> result10 = 2.0 - a;
        Dual2_3d<std::complex<double>> result11 = 2.0 * a;
        Dual2_3d<std::complex<double>> result12 = 2.0 / a;
        Dual2_3d<std::complex<double>> result13 = a + c1;
        Dual2_3d<std::complex<double>> result14 = a - c1;
        Dual2_3d<std::complex<double>> result15 = a * c1;
        Dual2_3d<std::complex<double>> result16 = a / c1;
        Dual2_3d<std::complex<double>> result17 = c1 + a;
        Dual2_3d<std::complex<double>> result18 = c1 - a;
        Dual2_3d<std::complex<double>> result19 = c1 * a;
        Dual2_3d<std::complex<double>> result20 = c1 / a;
        std::cout << "> Dual2_3d complex functioning...\n";
    }
}

/*
 * Test higher order polynomial with dual numbers
 * f(x,y,z) = yx^3 + x^2 + x + 1
 * This function tests the dual number implementation for a higher order polynomial.
 */
static void test_higher_order_polynomial_dual() {
    std::cout << "f(x,y,z) = yx^3 + x^2 + x + 1";
    bool all_match = true;
    for (double i = 0; i <= 5; i += 0.1) {
        for (double k = 0; k <= 5; k += 0.1) {
            // Analytical solution: f(x,y,z) = yx^3 + x^2 + x + 1
            double f_x = k * pow(i, 3) + pow(i, 2) + i + 1.0;
            // First derivatives
            double f_dx = 3 * k * pow(i, 2) + 2 * i + 1.0; // First  with respect to x
            double f_dy = pow(i, 3); // First derivative with respect to y
            double f_dz = 0.0; // First derivative with respect to z (not used in this case)

            // Second derivatives
            double f_dxx = 6 * k * i + 2.0; // Second derivative with respect to x
            double f_dyy = 0.0; // Second derivative with respect to y (not used in this case)
            double f_dzz = 0.0; // Second derivative with respect to z (not used in this case)
            double f_dxy = 3 * pow(i, 2); // Mixed derivative with respect to x and y
            double f_dxz = 0.0; // Mixed derivative with respect to x and z (not used in this case)
            double f_dyz = 0.0; // Mixed derivative with respect to y and z (not used in this case)

            // Create dual numbers
            //Dual1 numbers for dx
            Dual1<double> x_1(i, 1.0); // x = i, dx = 1.0
            Dual1<double> y_1(k, 0.0); // y = k, dy = 0.0

            //Dual1 numbers for dy
            Dual1<double> x_2(i, 0.0); // x = i, dx = 0.0
            Dual1<double> y_2(k, 1.0); // y = k, dy = 1.0

            // Dual2 numbers for x and y
            Dual2_3d<double> x(i, 1.0, 0.0, 0); // x = i, dx = 1.0, dy = 0.0, dz = 0.0
            Dual2_3d<double> y(k, 0.0, 1.0, 0); // y = k, dx = 0.0, dy = 1.0, dz = 0.0

            //Funciton
            // Dual1 version
            // Dual1 version with respect to x
            Dual1<double> x1 = y_1 * pow(x_1, 3) + pow(x_1, 2) + x_1 + 1.0;

            // Dual1 version with respect to y
            Dual1<double> x2 = y_2 * pow(x_2, 3) + pow(x_2, 2) + x_2 + 1.0;

            //Dual2 version
            Dual2_3d<double> x3 = y * pow(x, 3) + pow(x, 2) + x + 1.0;

            auto almost_equal = [](double a, double b, double eps = 1e-9) {
                return std::abs(a - b) < eps;
                };
            // Compare Dual2_3d to analytical
            if (!(almost_equal(x3.val, f_x) &&
                almost_equal(x3.dx, f_dx) &&
                almost_equal(x3.dy, f_dy) &&
                almost_equal(x3.dz, f_dz) &&
                almost_equal(x3.dxx, f_dxx) &&
                almost_equal(x3.dyy, f_dyy) &&
                almost_equal(x3.dzz, f_dzz) &&
                almost_equal(x3.dxy, f_dxy) &&
                almost_equal(x3.dxz, f_dxz) &&
                almost_equal(x3.dyz, f_dyz))) {
                all_match = false;
                std::cout << "x: " << x.val << ", y: " << y.val << "\n";
                std::cout << "  Mismatch in values!\n";
                std::cout << "  Expected f(x,y,z) = " << f_x << ", df/dx = " << f_dx << ", df/dy = " << f_dy << ", d^2f/dx^2 = " << f_dxx << "\n";
                std::cout << "           d^2f/dy^2 = " << f_dyy << ", d^2f/dz^2 = " << f_dzz << ", d^2f/dxy = " << f_dxy << ", d^2f/dxz = " << f_dxz << ", d^2f/dyz = " << f_dyz << "\n";
                std::cout << "        Got f(x,y,z) =      " << x3.val << ", df/dx = " << x3.dx << ", df/dy = " << x3.dy << ", d^2f/dx^2 = " << x3.dxx << "\n";
                std::cout << "           d^2f/dy^2 = " << x3.dyy << ", d^2f/dz^2 = " << x3.dzz << ", d^2f/dxy = " << x3.dxy << ", d^2f/dxz = " << x3.dxz << ", d^2f/dyz = " << x3.dyz << "\n\n";
            }
            // Compare Dual1 to Dual2_3d for value and dx
            if (!(almost_equal(x1.val, x3.val) && almost_equal(x1.d1, x3.dx))) {
                all_match = false;
                std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << "\n";
                std::cout << "  Dual1: val = " << x1.val << ", dx = " << x1.d1 << "\n";
                std::cout << "  Dual2_3d: val = " << x3.val << ", dx = " << x3.dx << "\n\n";
            }
            // Compare Dual1 to Dual2_3d for dy
            if (!(almost_equal(x2.d1, x3.dy))) {
                all_match = false;
                std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << "\n";
                std::cout << "  Dual1: dy = " << x2.d1 << "\n";
                std::cout << "  Dual2_3d: dy = " << x3.dy << "\n\n";
            }
        }
    }
    if (all_match) {
        std::cout << "  ( good )\n";
    }
}

/*
 * Test function for exp(xyz) using dual numbers
 * This function tests the exponential function with respect to three variables x, y, and z.
 * It computes the value and derivatives using both analytical methods and dual numbers.
 */
static void test_exp_func_dual() {
    //exp(xyz)
    std::cout << "f(x,y,z) = exp(xyz)\t     ";
    bool all_match = true;
    for (double i = 0; i <= 5; i += 0.1) {
        for (double j = 0; j <= 5; j += 0.1) {
            for (double k = 0; k <= 5; k += 0.1) {
                // Analytical solution: f(x,y,z) = exp(xyz)
                double f_x = exp(i * j * k);
                // First derivatives
                double f_dx = j * k * exp(i * j * k); // First derivative with respect to x
                double f_dy = i * k * exp(i * j * k); // First derivative with respect to y
                double f_dz = i * j * exp(i * j * k); // First derivative with respect to z
                // Second derivatives
                double f_dxx = j * j * k * k * exp(i * j * k); // Second derivative with respect to x
                double f_dyy = i * i * k * k * exp(i * j * k); // Second derivative with respect to y
                double f_dzz = i * i * j * j * exp(i * j * k); // Second derivative with respect to z
                double f_dxy = (k + i * j * k * k) * exp(i * j * k); // Mixed derivative with respect to x and y
                double f_dxz = (j + i * j * j * k) * exp(i * j * k); // Mixed derivative with respect to x and z
                double f_dyz = (i + i * i * j * k) * exp(i * j * k); // Mixed derivative with respect to y and z


                // Create dual numbers
                // Dual1 numbers for dx
                Dual1<double> x_1(i, 1.0); // x = i, dx = 1.0
                Dual1<double> y_1(j, 0.0); // y = j, dy = 0.0
                Dual1<double> z_1(k, 0.0); // z = k, dz = 0.0

                // Dual1 numbers for dy
                Dual1<double> x_2(i, 0.0); // x = i, dx = 0.0
                Dual1<double> y_2(j, 1.0); // y = j, dy = 1.0
                Dual1<double> z_2(k, 0.0); // z = k, dz = 0.0

                // Dual1 numbers for dz
                Dual1<double> x_3(i, 0.0); // x = i, dx = 0.0
                Dual1<double> y_3(j, 0.0); // y = j, dy = 0.0
                Dual1<double> z_3(k, 1.0); // z = k, dz = 1.0

                // Dual2 numbers for x, y, z
                Dual2_3d<double> x(i, 1.0, 0.0, 0); // x = i, dx = 1.0, dy = 0.0, dz = 0.0
                Dual2_3d<double> y(j, 0.0, 1.0, 0); // y = j, dx = 0.0, dy = 1.0, dz = 0.0
                Dual2_3d<double> z(k, 0.0, 0.0, 1.0); // z = k, dx = 0.0, dy = 0.0, dz = 1.0

                // Function
                // Dual1 version with respect to x
                Dual1<double> x1 = exp(x_1 * y_1 * z_1);

                // Dual1 version with respect to y
                Dual1<double> x2 = exp(x_2 * y_2 * z_2);

                // Dual1 version with respect to z
                Dual1<double> x3 = exp(x_3 * y_3 * z_3);

                // Dual2 version
                Dual2_3d<double> x4 = exp(x * y * z);

                auto almost_equal = [](double a, double b, double rel_eps = 1e-9, double abs_eps = 1e-9) {
                    return std::abs(a - b) <= std::max(rel_eps * std::max(std::abs(a), std::abs(b)), abs_eps);
                    };
                if (!(almost_equal(x4.val, f_x) &&
                    almost_equal(x4.dx, f_dx) &&
                    almost_equal(x4.dy, f_dy) &&
                    almost_equal(x4.dz, f_dz) &&
                    almost_equal(x4.dxx, f_dxx) &&
                    almost_equal(x4.dyy, f_dyy) &&
                    almost_equal(x4.dzz, f_dzz) &&
                    almost_equal(x4.dxy, f_dxy) &&
                    almost_equal(x4.dxz, f_dxz) &&
                    almost_equal(x4.dyz, f_dyz))) {
                    all_match = false;
                    std::cout << "x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Mismatch in values!\n";
                    if (!almost_equal(x4.val, f_x))
                        std::cout << "  f(x,y,z): expected " << f_x << ", actual " << x4.val << "\n";
                    if (!almost_equal(x4.dx, f_dx))
                        std::cout << "  df/dx: expected " << f_dx << ", actual " << x4.dx << "\n";
                    if (!almost_equal(x4.dy, f_dy))
                        std::cout << "  df/dy: expected " << f_dy << ", actual " << x4.dy << "\n";
                    if (!almost_equal(x4.dz, f_dz))
                        std::cout << "  df/dz: expected " << f_dz << ", actual " << x4.dz << "\n";
                    if (!almost_equal(x4.dxx, f_dxx))
                        std::cout << "  d^2f/dx^2: expected " << f_dxx << ", actual " << x4.dxx << "\n";
                    if (!almost_equal(x4.dyy, f_dyy))
                        std::cout << "  d^2f/dy^2: expected " << f_dyy << ", actual " << x4.dyy << "\n";
                    if (!almost_equal(x4.dzz, f_dzz))
                        std::cout << "  d^2f/dz^2: expected " << f_dzz << ", actual " << x4.dzz << "\n";
                    if (!almost_equal(x4.dxy, f_dxy))
                        std::cout << "  d^2f/dxy: expected " << f_dxy << ", actual " << x4.dxy << "\n";
                    if (!almost_equal(x4.dxz, f_dxz))
                        std::cout << "  d^2f/dxz: expected " << f_dxz << ", actual " << x4.dxz << "\n";
                    if (!almost_equal(x4.dyz, f_dyz))
                        std::cout << "  d^2f/dyz: expected " << f_dyz << ", actual " << x4.dyz << "\n";
                    std::cout << "\n";
                    break; // Break out of the loop on first mismatch
                }

                // Compare Dual1 to Dual2_3d for value and dx
                if (!(almost_equal(x1.val, x4.val) && almost_equal(x1.d1, x4.dx))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: val = " << x1.val << ", dx = " << x1.d1 << "\n";
                    std::cout << "  Dual2_3d: val = " << x4.val << ", dx = " << x4.dx << "\n\n";
                }

                // Compare Dual1 to Dual2_3d for dy
                if (!(almost_equal(x2.d1, x4.dy))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: dy = " << x2.d1 << "\n";
                    std::cout << "  Dual2_3d: dy = " << x4.dy << "\n\n";

                }
                // Compare Dual1 to Dual2_3d for dz
                if (!(almost_equal(x3.d1, x4.dz))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: dz = " << x3.d1 << "\n";
                    std::cout << "  Dual2_3d: dz = " << x4.dz << "\n\n";
                }
            }
        }
    }
    if (all_match) {
        std::cout << "  ( good )\n";
    }
}

/*
 * Test function for pow(xyz, 3) using dual numbers.
 * This function tests the power function with a constant base (3) and variables x, y, z.
 * It calculates the value and derivatives using both Dual1 and Dual2_3d representations.
 * The results are compared against analytical solutions to ensure correctness.
 */
static void test_pow_dual_scalar_func_dual() {
    //pow(xyz, 3)
    std::cout << "f(x,y,z) = (xyz)^3\t     ";
    bool all_match = true;

    for (double i = 0; i <= 5; i += 0.1) {
        for (double j = 0; j <= 5; j += 0.1) {
            for (double k = 0; k <= 5; k += 0.1) {
                // Analytical solution: f(x,y,z) = pow(xyz, 3)
                double f_x = pow(i * j * k, 3);
                // First derivatives
                double f_dx = 3 * pow(i, 2) * pow(j, 3) * pow(k, 3); // First derivative with respect to x is 3 * (xyz)^2 * yz
                double f_dy = 3 * pow(i, 3) * pow(j, 2) * pow(k, 3); // First derivative with respect to y is 3 * (xyz)^2 * xz
                double f_dz = 3 * pow(i, 3) * pow(j, 3) * pow(k, 2); // First derivative with respect to z is 3 * (xyz)^2 * xy
                // Second derivatives
                double f_dxx = 6 * i * pow(j, 3) * pow(k, 3); // Second derivative with respect to x is 6 * (xyz)^2 * (yz)^2
                double f_dyy = 6 * pow(i, 3) * j * pow(k, 3); // Second derivative with respect to y is 6 * (xyz)^2 * (xz)^2
                double f_dzz = 6 * pow(i, 3) * pow(j, 3) * k; // Second derivative with respect to z is 6 * (xyz)^2 * (xy)^2
                double f_dxy = 9 * pow(i, 2) * pow(j, 2) * pow(k, 3); // Mixed derivative with respect to x and y is 9 * (xyz)^2 * yz * xz
                double f_dxz = 9 * pow(i, 2) * pow(j, 3) * pow(k, 2); // Mixed derivative with respect to x and z is 9 * (xyz)^2 * yz * xy
                double f_dyz = 9 * pow(i, 3) * pow(j, 2) * pow(k, 2); // Mixed derivative with respect to y and z is 9 * (xyz)^2 * xz * xy

                // Create dual numbers
                // Dual1 numbers for dx
                Dual1<double> x_1(i, 1.0); // x = i, dx = 1.0
                Dual1<double> y_1(j, 0.0); // y = j, dy = 0.0
                Dual1<double> z_1(k, 0.0); // z = k, dz = 0.0

                // Dual1 numbers for dy
                Dual1<double> x_2(i, 0.0); // x = i, dx = 0.0
                Dual1<double> y_2(j, 1.0); // y = j, dy = 1.0
                Dual1<double> z_2(k, 0.0); // z = k, dz = 0.0

                // Dual1 numbers for dz
                Dual1<double> x_3(i, 0.0); // x = i, dx = 0.0
                Dual1<double> y_3(j, 0.0); // y = j, dy = 0.0
                Dual1<double> z_3(k, 1.0); // z = k, dz = 1.0

                // Dual2 numbers for x, y, z
                Dual2_3d<double> x(i, 1.0, 0.0, 0); // x = i, dx = 1.0
                Dual2_3d<double> y(j, 0.0, 1.0, 0); // y = j, dx = 0.0, dy = 1.0, dz = 0.0
                Dual2_3d<double> z(k, 0.0, 0.0, 1.0); // z = k, dx = 0.0, dy = 0.0, dz = 1.0

                // Function
                // Dual1 version with respect to x
                Dual1<double> x1 = pow(x_1 * y_1 * z_1, 3);

                // Dual1 version with respect to y
                Dual1<double> x2 = pow(x_2 * y_2 * z_2, 3);

                // Dual1 version with respect to z
                Dual1<double> x3 = pow(x_3 * y_3 * z_3, 3);

                // Dual2 version
                Dual2_3d<double> x4 = pow(x * y * z, 3);
                auto almost_equal = [](double a, double b, double rel_eps = 1e-9, double abs_eps = 1e-9) {
                    return std::abs(a - b) <= std::max(rel_eps * std::max(std::abs(a), std::abs(b)), abs_eps);
                    };
                if (!(almost_equal(x4.val, f_x) &&
                    almost_equal(x4.dx, f_dx) &&
                    almost_equal(x4.dy, f_dy) &&
                    almost_equal(x4.dz, f_dz) &&
                    almost_equal(x4.dxx, f_dxx) &&
                    almost_equal(x4.dyy, f_dyy) &&
                    almost_equal(x4.dzz, f_dzz) &&
                    almost_equal(x4.dxy, f_dxy) &&
                    almost_equal(x4.dxz, f_dxz) &&
                    almost_equal(x4.dyz, f_dyz))) {
                    all_match = false;
                    std::cout << "x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Mismatch in values!\n";
                    if (!almost_equal(x4.val, f_x))
                        std::cout << "  f(x,y,z): expected " << f_x << ", actual " << x4.val << "\n";
                    if (!almost_equal(x4.dx, f_dx))
                        std::cout << "  df/dx: expected " << f_dx << ", actual " << x4.dx << "\n";
                    if (!almost_equal(x4.dy, f_dy))
                        std::cout << "  df/dy: expected " << f_dy << ", actual " << x4.dy << "\n";
                    if (!almost_equal(x4.dz, f_dz))
                        std::cout << "  df/dz: expected " << f_dz << ", actual " << x4.dz << "\n";
                    if (!almost_equal(x4.dxx, f_dxx))
                        std::cout << "  d^2f/dx^2: expected " << f_dxx << ", actual " << x4.dxx << "\n";
                    if (!almost_equal(x4.dyy, f_dyy))
                        std::cout << "  d^2f/dy^2: expected " << f_dyy << ", actual " << x4.dyy << "\n";
                    if (!almost_equal(x4.dzz, f_dzz))
                        std::cout << "  d^2f/dz^2: expected " << f_dzz << ", actual " << x4.dzz << "\n";
                    if (!almost_equal(x4.dxy, f_dxy))
                        std::cout << "  d^2f/dxy: expected " << f_dxy << ", actual " << x4.dxy << "\n";
                    if (!almost_equal(x4.dxz, f_dxz))
                        std::cout << "  d^2f/dxz: expected " << f_dxz << ", actual " << x4.dxz << "\n";
                    if (!almost_equal(x4.dyz, f_dyz))
                        std::cout << "  d^2f/dyz: expected " << f_dyz << ", actual " << x4.dyz << "\n";
                    std::cout << "\n";
                }

                // Compare Dual1 to Dual2_3d for value and dx
                if (!(almost_equal(x1.val, x4.val) && almost_equal(x1.d1, x4.dx))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: val = " << x1.val << ", dx = " << x1.d1 << "\n";
                    std::cout << "  Dual2_3d: val = " << x4.val << ", dx = " << x4.dx << "\n\n";
                }
                // Compare Dual1 to Dual2_3d for dy
                if (!(almost_equal(x2.d1, x4.dy))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: dy = " << x2.d1 << "\n";
                    std::cout << "  Dual2_3d: dy = " << x4.dy << "\n\n";
                }
                // Compare Dual1 to Dual2_3d for dz
                if (!(almost_equal(x3.d1, x4.dz))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: dz = " << x3.d1 << "\n";
                    std::cout << "  Dual2_3d: dz = " << x4.dz << "\n\n";
                }
            }
        }
    }
    if (all_match) {
        std::cout << "  ( good )\n";
    }
}

/*
 * Test function for pow(3, xyz) using dual numbers.
 * This function tests the power function with a constant base (3) and variables x, y, z.
 * It calculates the value and derivatives using both Dual1 and Dual2_3d representations.
 * The results are compared against analytical solutions to ensure correctness.
 */
static void test_pow_scalar_dual_func_dual() {
    //pow(3,xyz)
    std::cout << "f(x,y,z) = 3^(xyz)\t      ";
    bool all_match = true;
    for (double i = 0; i <= 5; i += 0.1) {
        for (double j = 0; j <= 5; j += 0.1) {
            for (double k = 0; k <= 5; k += 0.1) {
                // Analytical solution: f(x,y,z) = pow(xyz, 3)
                double f_x = pow(3, i * j * k);
                // First derivatives
                double f_dx = log(3) * pow(3, i * j * k) * j * k; // First derivative with respect to x is log(3) * 3^(xyz) * yz
                double f_dy = log(3) * pow(3, i * j * k) * i * k; // First derivative with respect to y is log(3) * 3^(xyz) * xz
                double f_dz = log(3) * pow(3, i * j * k) * i * j; // First derivative with respect to z is log(3) * 3^(xyz) * xy
                // Second derivatives
                double f_dxx = log(3) * log(3) * pow(3, i * j * k) * pow(j, 2) * pow(k, 2);// Second derivative with respect to x is log(3)^2 * 3^(xyz) * (yz)^2
                double f_dyy = log(3) * log(3) * pow(3, i * j * k) * pow(i, 2) * pow(k, 2); // Second derivative with respect to y is log(3)^2 * 3^(xyz) * (xz)^2
                double f_dzz = log(3) * log(3) * pow(3, i * j * k) * pow(i, 2) * pow(j, 2); // Second derivative with respect to z is log(3)^2 * 3^(xyz) * (xy)^2
                double f_dxy = log(3) * k * (log(3) * i * j * k + 1) * pow(3, i * j * k); // Mixed derivative with respect to x and y is log(3)^2 * 3^(xyz) * yz * xz
                double f_dxz = log(3) * j * (log(3) * i * j * k + 1) * pow(3, i * j * k); // Mixed derivative with respect to x and z is log(3)^2 * 3^(xyz) * yz * xy
                double f_dyz = log(3) * i * (log(3) * i * j * k + 1) * pow(3, i * j * k); // Mixed derivative with respect to y and z is log(3)^2 * 3^(xyz) * xz * xy

                // Create dual numbers
                //Dual1 numbers for dx
                Dual1<double> x_1(i, 1.0); // x = i, dx = 1.0
                Dual1<double> y_1(j, 0.0); // y = j, dy = 0.0
                Dual1<double> z_1(k, 0.0); // z = k, dz = 0.0

                //Dual1 numbers for dy
                Dual1<double> x_2(i, 0.0); // x = i, dx = 0.0
                Dual1<double> y_2(j, 1.0); // y = j, dy = 1.0
                Dual1<double> z_2(k, 0.0); // z = k, dz = 0.0

                //Dual1 numbers for dz
                Dual1<double> x_3(i, 0.0); // x = i, dx = 0.0
                Dual1<double> y_3(j, 0.0); // y = j, dy = 0.0
                Dual1<double> z_3(k, 1.0); // z = k, dz = 1.0

                // Dual2 numbers for x, y, z
                Dual2_3d<double> x(i, 1.0, 0.0, 0); // x = i, dx = 1.0
                Dual2_3d<double> y(j, 0.0, 1.0, 0); // y = j, dx = 0.0, dy = 1.0, dz = 0.0
                Dual2_3d<double> z(k, 0.0, 0.0, 1.0); // z = k, dx = 0.0, dy = 0.0, dz = 1.0

                // Function
                // Dual1 version with respect to x
                Dual1<double> x1 = pow(3, x_1 * y_1 * z_1);

                // Dual1 version with respect to y
                Dual1<double> x2 = pow(3, x_2 * y_2 * z_2);

                // Dual1 version with respect to z
                Dual1<double> x3 = pow(3, x_3 * y_3 * z_3);

                // Dual2 version
                Dual2_3d<double> x4 = pow(3, x * y * z);
                auto almost_equal = [](double a, double b, double rel_eps = 1e-9, double abs_eps = 1e-9) {
                    return std::abs(a - b) <= std::max(rel_eps * std::max(std::abs(a), std::abs(b)), abs_eps);
                    };
                if (!(almost_equal(x4.val, f_x) &&
                    almost_equal(x4.dx, f_dx) &&
                    almost_equal(x4.dy, f_dy) &&
                    almost_equal(x4.dz, f_dz) &&
                    almost_equal(x4.dxx, f_dxx) &&
                    almost_equal(x4.dyy, f_dyy) &&
                    almost_equal(x4.dzz, f_dzz) &&
                    almost_equal(x4.dxy, f_dxy) &&
                    almost_equal(x4.dxz, f_dxz) &&
                    almost_equal(x4.dyz, f_dyz))) {
                    all_match = false;
                    std::cout << "\nx: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Mismatch in values!\n";
                    if (!almost_equal(x4.val, f_x))
                        std::cout << "  f(x,y,z): expected " << f_x << ", actual " << x4.val << "\n";
                    if (!almost_equal(x4.dx, f_dx))
                        std::cout << "  df/dx: expected " << f_dx << ", actual " << x4.dx << "\n";
                    if (!almost_equal(x4.dy, f_dy))
                        std::cout << "  df/dy: expected " << f_dy << ", actual " << x4.dy << "\n";
                    if (!almost_equal(x4.dz, f_dz))
                        std::cout << "  df/dz: expected " << f_dz << ", actual " << x4.dz << "\n";
                    if (!almost_equal(x4.dxx, f_dxx))
                        std::cout << "  d^2f/dx^2: expected " << f_dxx << ", actual " << x4.dxx << "\n";
                    if (!almost_equal(x4.dyy, f_dyy))
                        std::cout << "  d^2f/dy^2: expected " << f_dyy << ", actual " << x4.dyy << "\n";
                    if (!almost_equal(x4.dzz, f_dzz))
                        std::cout << "  d^2f/dz^2: expected " << f_dzz << ", actual " << x4.dzz << "\n";
                    if (!almost_equal(x4.dxy, f_dxy))
                        std::cout << "  d^2f/dxy: expected " << f_dxy << ", actual " << x4.dxy << "\n";
                    if (!almost_equal(x4.dxz, f_dxz))
                        std::cout << "  d^2f/dxz: expected " << f_dxz << ", actual " << x4.dxz << "\n";
                    if (!almost_equal(x4.dyz, f_dyz))
                        std::cout << "  d^2f/dyz: expected " << f_dyz << ", actual " << x4.dyz << "\n";
                    std::cout << "\n";
                }

                // Compare Dual1 to Dual2_3d for value and dx
                if (!(almost_equal(x1.val, x4.val) && almost_equal(x1.d1, x4.dx))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: val = " << x1.val << ", dx = " << x1.d1 << "\n";
                    std::cout << "  Dual2_3d: val = " << x4.val << ", dx = " << x4.dx << "\n\n";
                }
                // Compare Dual1 to Dual2_3d for dy
                if (!(almost_equal(x2.d1, x4.dy))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: dy = " << x2.d1 << "\n";
                    std::cout << "  Dual2_3d: dy = " << x4.dy << "\n\n";
                }
                // Compare Dual1 to Dual2_3d for dz
                if (!(almost_equal(x3.d1, x4.dz))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: dz = " << x3.d1 << "\n";
                    std::cout << "  Dual2_3d: dz = " << x4.dz << "\n\n";
                }
            }
        }

    }
    if (all_match) {
        std::cout << " ( good )\n";
    }
}

/*
 * Test function for pow with dual numbers and a function of dual numbers
 * f(x,y,z) = (x + y + z)^(x)
 */
static void test_pow_dual_dual_func_dual() {
    std::cout << "f(x,y,z) = (x + y + z)^(x)   ";
    bool all_match = true;
    for (double i = 2; i <= 5; i += 0.1) {
        for (double j = 1; j <= 5; j += 0.1) {
            for (double k = 1; k <= 5; k += 0.1) {
                // Analytical solution: f(x,y,z) = (x + y + z)^(x)
                double s = i + j + k;
                double f_x = pow(s, i);
                double log_s = std::log(s);
                // First derivatives
                double f_dx = f_x * (log_s + i / s);
                double f_dy = i * pow(s, i - 1);
                double f_dz = i * pow(s, i - 1);
                // Second derivatives
                double f_dxx = f_x * ((log_s + (i / s)) * (log_s + (i / s)) + ((2.0 / s) - (i / (s * s))));
                double f_dyy = i * (i - 1) * pow(s, i - 2);
                double f_dzz = i * (i - 1) * pow(s, i - 2);
                double f_dxy = i * pow(s, i - 1) * (log_s + (i - 1) / s) + pow(s, i - 1);
                double f_dxz = i * pow(s, i - 1) * (log_s + (i - 1) / s) + pow(s, i - 1);
                double f_dyz = i * (i - 1) * pow(s, i - 2);

                // Create dual numbers
                // Dual1 numbers for dx
                Dual1<double> x_1(i, 1.0); // x = i, dx = 1.0
                Dual1<double> y_1(j, 0.0); // y = j, dy = 0.0
                Dual1<double> z_1(k, 0.0); // z = k, dz = 0.0

                // Dual1 numbers for dy
                Dual1<double> x_2(i, 0.0); // x = i, dx = 0.0
                Dual1<double> y_2(j, 1.0); // y = j, dy = 1.0
                Dual1<double> z_2(k, 0.0); // z = k, dz = 0.0

                // Dual1 numbers for dz
                Dual1<double> x_3(i, 0.0); // x = i, dx = 0.0
                Dual1<double> y_3(j, 0.0); // y = j, dy = 0.0
                Dual1<double> z_3(k, 1.0); // z = k, dz = 1.0

                // Dual2 numbers for x, y, z
                Dual2_3d<double> x(i, 1.0, 0.0, 0); // x = i, dx = 1.0
                Dual2_3d<double> y(j, 0.0, 1.0, 0); // y = j, dx = 0.0, dy = 1.0, dz = 0.0
                Dual2_3d<double> z(k, 0.0, 0.0, 1.0); // z = k, dx = 0.0, dy = 0.0, dz = 1.0

                // Function
                // Dual1 version with respect to x
                Dual1<double> x1 = pow(x_1 + y_1 + z_1, x_1);

                // Dual1 version with respect to y
                Dual1<double> x2 = pow(x_2 + y_2 + z_2, x_2);

                // Dual1 version with respect to z
                Dual1<double> x3 = pow(x_3 + y_3 + z_3, x_3);

                // Dual2 version
                Dual2_3d<double> x4 = pow(x + y + z, x);

                auto almost_equal = [](double a, double b, double rel_eps = 1e-9, double abs_eps = 1e-9) {
                    return std::abs(a - b) <= std::max(rel_eps * std::max(std::abs(a), std::abs(b)), abs_eps);
                    };
                if (!(almost_equal(x4.val, f_x) &&
                    almost_equal(x4.dx, f_dx) &&
                    almost_equal(x4.dy, f_dy) &&
                    almost_equal(x4.dz, f_dz) &&
                    almost_equal(x4.dxx, f_dxx) &&
                    almost_equal(x4.dyy, f_dyy) &&
                    almost_equal(x4.dzz, f_dzz) &&
                    almost_equal(x4.dxy, f_dxy) &&
                    almost_equal(x4.dxz, f_dxz) &&
                    almost_equal(x4.dyz, f_dyz))) {
                    all_match = false;
                    std::cout << "x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Mismatch in values!\n";
                    if (!almost_equal(x4.val, f_x))
                        std::cout << "  f(x,y,z): expected " << f_x << ", actual " << x4.val << "\n";
                    if (!almost_equal(x4.dx, f_dx))
                        std::cout << "  df/dx: expected " << f_dx << ", actual " << x4.dx << "\n";
                    if (!almost_equal(x4.dy, f_dy))
                        std::cout << "  df/dy: expected " << f_dy << ", actual " << x4.dy << "\n";
                    if (!almost_equal(x4.dz, f_dz))
                        std::cout << "  df/dz: expected " << f_dz << ", actual " << x4.dz << "\n";
                    if (!almost_equal(x4.dxx, f_dxx))
                        std::cout << "  d^2f/dx^2: expected " << f_dxx << ", actual " << x4.dxx << "\n";
                    if (!almost_equal(x4.dyy, f_dyy))
                        std::cout << "  d^2f/dy^2: expected " << f_dyy << ", actual " << x4.dyy << "\n";
                    if (!almost_equal(x4.dzz, f_dzz))
                        std::cout << "  d^2f/dz^2: expected " << f_dzz << ", actual " << x4.dzz << "\n";
                    if (!almost_equal(x4.dxy, f_dxy))
                        std::cout << "  d^2f/dxy: expected " << f_dxy << ", actual " << x4.dxy << "\n";
                    if (!almost_equal(x4.dxz, f_dxz))
                        std::cout << "  d^2f/dxz: expected " << f_dxz << ", actual " << x4.dxz << "\n";
                    if (!almost_equal(x4.dyz, f_dyz))
                        std::cout << "  d^2f/dyz: expected " << f_dyz << ", actual " << x4.dyz << "\n";
                    std::cout << "\n";

                }

                // Compare Dual1 to Dual2_3d for value and dx
                if (!(almost_equal(x1.val, x4.val) && almost_equal(x1.d1, x4.dx))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: val = " << x1.val << ", dx = " << x1.d1 << "\n";
                    std::cout << "  Dual2_3d: val = " << x4.val << ", dx = " << x4.dx << "\n\n";
                }
                // Compare Dual1 to Dual2_3d for dy
                if (!(almost_equal(x2.d1, x4.dy))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: dy = " << x2.d1 << "\n";
                    std::cout << "  Dual2_3d: dy = " << x4.dy << "\n\n";
                }
                // Compare Dual1 to Dual2_3d for dz
                if (!(almost_equal(x3.d1, x4.dz))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: dz = " << x3.d1 << "\n";
                    std::cout << "  Dual2_3d: dz = " << x4.dz << "\n\n";
                }
            }


        }
    }
    if (all_match)
        std::cout << "  ( good )\n";
}

/*
 * Test function for log function with dual numbers
 * This function tests the log function with three variables (x, y, z)
 * using both Dual1 and Dual2_3d representations.
 */
static void test_log_func_dual() {
    std::cout << "f(x,y,z) = log((x+y+z)       ";
    bool all_match = true;
    for (double i = 0.1; i <= 5; i += 0.1) { // Start from 0.1 to avoid log(0)
        for (double j = 0.1; j <= 5; j += 0.1) {
            for (double k = 0.1; k <= 5; k += 0.1) {
                // Analytical solution: f(x,y,z) = log(xyz)
                double f_x = log(i + j + k);
                // First derivatives
                double f_dx = 1 / (i + j + k); // First derivative with respect to x
                double f_dy = 1 / (i + j + k); // First derivative with respect to y
                double f_dz = 1 / (i + j + k); // First derivative with respect to z
                // Second derivatives
                double f_dxx = -1 / ((i + j + k) * (i + j + k)); // Second derivative with respect to x
                double f_dyy = -1 / ((i + j + k) * (i + j + k)); // Second derivative with respect to y
                double f_dzz = -1 / ((i + j + k) * (i + j + k)); // Second derivative with respect to z
                double f_dxy = -1 / ((i + j + k) * (i + j + k)); // Mixed derivative with respect to x and y
                double f_dxz = -1 / ((i + j + k) * (i + j + k)); // Mixed derivative with respect to x and z
                double f_dyz = -1 / ((i + j + k) * (i + j + k)); // Mixed derivative with respect to y and z

                // Create dual numbers
                // Dual1 numbers for dx
                Dual1<double> x_1(i, 1.0); // x = i, dx = 1.0
                Dual1<double> y_1(j, 0.0); // y = j, dy = 0.0
                Dual1<double> z_1(k, 0.0); // z = k, dz = 0.0

                // Dual1 numbers for dy
                Dual1<double> x_2(i, 0.0); // x = i, dx = 0.0
                Dual1<double> y_2(j, 1.0); // y = j, dy = 1.0
                Dual1<double> z_2(k, 0.0); // z = k, dz = 0.0

                // Dual1 numbers for dz
                Dual1<double> x_3(i, 0.0); // x = i, dx = 0.0
                Dual1<double> y_3(j, 0.0); // y = j, dy = 0.0
                Dual1<double> z_3(k, 1.0); // z = k, dz = 1.0

                // Dual2 numbers for x, y, z
                Dual2_3d<double> x(i, 1.0, 0.0, 0); // x = i, dx = 1.0
                Dual2_3d<double> y(j, 0.0, 1.0, 0); // y = j, dx = 0.0, dy = 1.0, dz = 0.0
                Dual2_3d<double> z(k, 0.0, 0.0, 1.0); // z

                // Function
                // Dual1 version with respect to x
                Dual1<double> x1 = log(x_1 + y_1 + z_1);

                // Dual1 version with respect to y
                Dual1<double> x2 = log(x_2 + y_2 + z_2);

                // Dual1 version with respect to z
                Dual1<double> x3 = log(x_3 + y_3 + z_3);

                // Dual2 version
                Dual2_3d<double> x4 = log(x + y + z);
                auto almost_equal = [](double a, double b, double rel_eps = 1e-9, double abs_eps = 1e-9) {
                    return std::abs(a - b) <= std::max(rel_eps * std::max(std::abs(a), std::abs(b)), abs_eps);
                    };
                if (!(almost_equal(x4.val, f_x) &&
                    almost_equal(x4.dx, f_dx) &&
                    almost_equal(x4.dy, f_dy) &&
                    almost_equal(x4.dz, f_dz) &&
                    almost_equal(x4.dxx, f_dxx) &&
                    almost_equal(x4.dyy, f_dyy) &&
                    almost_equal(x4.dzz, f_dzz) &&
                    almost_equal(x4.dxy, f_dxy) &&
                    almost_equal(x4.dxz, f_dxz) &&
                    almost_equal(x4.dyz, f_dyz))) {
                    all_match = false;
                    std::cout << "x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Mismatch in values!\n";
                    if (!almost_equal(x4.val, f_x))
                        std::cout << "  f(x,y,z): expected " << f_x << ", actual " << x4.val << "\n";
                    if (!almost_equal(x4.dx, f_dx))
                        std::cout << "  df/dx: expected " << f_dx << ", actual " << x4.dx << "\n";
                    if (!almost_equal(x4.dy, f_dy))
                        std::cout << "  df/dy: expected " << f_dy << ", actual " << x4.dy << "\n";
                    if (!almost_equal(x4.dz, f_dz))
                        std::cout << "  df/dz: expected " << f_dz << ", actual " << x4.dz << "\n";
                    if (!almost_equal(x4.dxx, f_dxx))
                        std::cout << "  d^2f/dx^2: expected " << f_dxx << ", actual " << x4.dxx << "\n";
                    if (!almost_equal(x4.dyy, f_dyy))
                        std::cout << "  d^2f/dy^2: expected " << f_dyy << ", actual " << x4.dyy << "\n";
                    if (!almost_equal(x4.dzz, f_dzz))
                        std::cout << "  d^2f/dz^2: expected " << f_dzz << ", actual " << x4.dzz << "\n";
                    if (!almost_equal(x4.dxy, f_dxy))
                        std::cout << "  d^2f/dxy: expected " << f_dxy << ", actual " << x4.dxy << "\n";
                    if (!almost_equal(x4.dxz, f_dxz))
                        std::cout << "  d^2f/dxz: expected " << f_dxz << ", actual " << x4.dxz << "\n";
                    if (!almost_equal(x4.dyz, f_dyz))
                        std::cout << "  d^2f/dyz: expected " << f_dyz << ", actual " << x4.dyz << "\n";
                    std::cout << "\n";
                }

                // Compare Dual1 to Dual2_3d for value and dx
                if (!(almost_equal(x1.val, x4.val) && almost_equal(x1.d1, x4.dx))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: val = " << x1.val << ", dx = " << x1.d1 << "\n";
                    std::cout << "  Dual2_3d: val = " << x4.val << ", dx = " << x4.dx << "\n\n";
                }

                // Compare Dual1 to Dual2_3d for dy
                if (!(almost_equal(x2.d1, x4.dy))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: dy = " << x2.d1 << "\n";
                    std::cout << "  Dual2_3d: dy = " << x4.dy << "\n\n";
                }
                // Compare Dual1 to Dual2_3d for dz
                if (!(almost_equal(x3.d1, x4.dz))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: dz = " << x3.d1 << "\n";
                    std::cout << "  Dual2_3d: dz = " << x4.dz << "\n\n";
                }

            }
        }
    }
    if (all_match) {
        std::cout << "  ( good )\n";
    }
}

/*
 * Test function for sqrt with three variables using dual numbers.
 * Function f(x,y,z) = sqrt(x + y + z)
 * Compares the results of Analytical derivatives with Dual1 and Dual2_3d implementations.
 * Checks the value, first derivatives, second derivatives, and mixed derivatives.
 * The function iterates through a range of values for x, y, and z.
 */
static void test_sqrt_func_dual() {
    std::cout << "f(x,y,z) = sqrt(x + y + z)   ";
    bool all_match = true;
    for (double i = 0.1; i <= 5; i += 0.1) { // Start from 0.1 to avoid sqrt(0)
        for (double j = 0.1; j <= 5; j += 0.1) {
            for (double k = 0.1; k <= 5; k += 0.1) {
                // Analytical solution: f(x,y,z) = sqrt(xyz)
                double f_x = sqrt(i + j + k);
                // First derivatives
                double f_dx = 0.5 / sqrt(i + j + k); // First derivative with respect to x
                double f_dy = 0.5 / sqrt(i + j + k); // First derivative with respect to y
                double f_dz = 0.5 / sqrt(i + j + k); // First derivative with respect to z
                // Second derivatives
                double f_dxx = -0.25 / ((i + j + k) * sqrt(i + j + k)); // Second derivative with respect to x
                double f_dyy = -0.25 / ((i + j + k) * sqrt(i + j + k)); // Second derivative with respect to y
                double f_dzz = -0.25 / ((i + j + k) * sqrt(i + j + k)); // Second derivative with respect to z
                double f_dxy = -0.25 / ((i + j + k) * sqrt(i + j + k)); // Mixed derivative with respect to x and y
                double f_dxz = -0.25 / ((i + j + k) * sqrt(i + j + k)); // Mixed derivative with respect to x and z
                double f_dyz = -0.25 / ((i + j + k) * sqrt(i + j + k)); // Mixed derivative with respect to y and z

                // Create dual numbers
                //Dual1 numbers for dx
                Dual1<double> x_1(i, 1.0); // x = i, dx = 1.0
                Dual1<double> y_1(j, 0.0); // y = j, dy = 0.0
                Dual1<double> z_1(k, 0.0); // z = k, dz = 0.0

                //Dual1 numbers for dy
                Dual1<double> x_2(i, 0.0); // x = i, dx = 0.0
                Dual1<double> y_2(j, 1.0); // y = j, dy = 1.0
                Dual1<double> z_2(k, 0.0); // z = k, dz = 0.0

                //Dual1 numbers for dz
                Dual1<double> x_3(i, 0.0); // x = i, dx = 0.0
                Dual1<double> y_3(j, 0.0); // y = j, dy = 0.0
                Dual1<double> z_3(k, 1.0); // z = k, dz = 1.0

                // Dual2 numbers for x, y, z
                Dual2_3d<double> x(i, 1.0, 0.0, 0); // x = i, dx = 1.0
                Dual2_3d<double> y(j, 0.0, 1.0, 0); // y = j, dx = 0.0, dy = 1.0, dz = 0.0
                Dual2_3d<double> z(k, 0.0, 0.0, 1.0); // z = k, dx = 0.0, dy = 0.0, dz = 1.0

                // Function
                // Dual1 version with respect to x
                Dual1<double> x1 = sqrt(x_1 + y_1 + z_1);

                // Dual1 version with respect to y
                Dual1<double> x2 = sqrt(x_2 + y_2 + z_2);

                // Dual1 version with respect to z
                Dual1<double> x3 = sqrt(x_3 + y_3 + z_3);

                // Dual2 version
                Dual2_3d<double> x4 = sqrt(x + y + z);
                auto almost_equal = [](double a, double b, double rel_eps = 1e-9, double abs_eps = 1e-9) {
                    return std::abs(a - b) <= std::max(rel_eps * std::max(std::abs(a), std::abs(b)), abs_eps);
                    };
                if (!(almost_equal(x4.val, f_x) &&
                    almost_equal(x4.dx, f_dx) &&
                    almost_equal(x4.dy, f_dy) &&
                    almost_equal(x4.dz, f_dz) &&
                    almost_equal(x4.dxx, f_dxx) &&
                    almost_equal(x4.dyy, f_dyy) &&
                    almost_equal(x4.dzz, f_dzz) &&
                    almost_equal(x4.dxy, f_dxy) &&
                    almost_equal(x4.dxz, f_dxz) &&
                    almost_equal(x4.dyz, f_dyz))) {
                    all_match = false;
                    std::cout << "x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Mismatch in values!\n";
                    if (!almost_equal(x4.val, f_x))
                        std::cout << "  f(x,y,z): expected " << f_x << ", actual " << x4.val << "\n";
                    if (!almost_equal(x4.dx, f_dx))
                        std::cout << "  df/dx: expected " << f_dx << ", actual " << x4.dx << "\n";
                    if (!almost_equal(x4.dy, f_dy))
                        std::cout << "  df/dy: expected " << f_dy << ", actual " << x4.dy << "\n";
                    if (!almost_equal(x4.dz, f_dz))
                        std::cout << "  df/dz: expected " << f_dz << ", actual " << x4.dz << "\n";
                    if (!almost_equal(x4.dxx, f_dxx))
                        std::cout << "  d^2f/dx^2: expected " << f_dxx << ", actual " << x4.dxx << "\n";
                    if (!almost_equal(x4.dyy, f_dyy))
                        std::cout << "  d^2f/dy^2: expected " << f_dyy << ", actual " << x4.dyy << "\n";
                    if (!almost_equal(x4.dzz, f_dzz))
                        std::cout << "  d^2f/dz^2: expected " << f_dzz << ", actual " << x4.dzz << "\n";
                    if (!almost_equal(x4.dxy, f_dxy))
                        std::cout << "  d^2f/dxy: expected " << f_dxy << ", actual " << x4.dxy << "\n";
                    if (!almost_equal(x4.dxz, f_dxz))
                        std::cout << "  d^2f/dxz: expected " << f_dxz << ", actual " << x4.dxz << "\n";
                    if (!almost_equal(x4.dyz, f_dyz))
                        std::cout << "  d^2f/dyz: expected " << f_dyz << ", actual " << x4.dyz << "\n";
                    std::cout << "\n";
                }

                // Compare Dual1 to Dual2_3d for value and dx
                if (!(almost_equal(x1.val, x4.val) && almost_equal(x1.d1, x4.dx))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: val = " << x1.val << ", dx = " << x1.d1 << "\n";
                    std::cout << "  Dual2_3d: val = " << x4.val << ", dx = " << x4.dx << "\n\n";
                }
                // Compare Dual1 to Dual2_3d for dy
                if (!(almost_equal(x2.d1, x4.dy))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: dy = " << x2.d1 << "\n";
                    std::cout << "  Dual2_3d: dy = " << x4.dy << "\n\n";
                }
                // Compare Dual1 to Dual2_3d for dz
                if (!(almost_equal(x3.d1, x4.dz))) {
                    all_match = false;
                    std::cout << "Dual1 vs Dual2_3d mismatch at x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Dual1: dz = " << x3.d1 << "\n";
                    std::cout << "  Dual2_3d: dz = " << x4.dz << "\n\n";
                }
            }
        }
    }
    if (all_match) {
        std::cout << "  ( good )\n";
    }
}

/*
 * Test the chain rule and symmetry of derivatives in Dual1 and Dual2_3d classes.
 * The function f(x, y, z) = sin(x * y + z^2) is used to verify the correctness of derivatives.
 */
static void test_chain_rule_and_symmetry() {
    double x0 = 1.2, y0 = -0.7, z0 = 0.5;

    // Create dual numbers
    // Dual1 numbers for dx
    Dual1<double> x1(x0, 1.0); // x = x0, dx = 1.0
    Dual1<double> y1(y0, 0.0); // y = y0, dy = 0.0
    Dual1<double> z1(z0, 0.0); // z = z0, dz = 0.0

    // Dual1 numbers for dy
    Dual1<double> x2(x0, 0.0); // x = x0, dx = 0.0
    Dual1<double> y2(y0, 1.0); // y = y0, dy = 1.0
    Dual1<double> z2(z0, 0.0); // z = z0, dz = 0.0

    // Dual1 numbers for dz
    Dual1<double> x3(x0, 0.0); // x = x0, dx = 0.0
    Dual1<double> y3(y0, 0.0); // y = y0, dy = 0.0
    Dual1<double> z3(z0, 1.0); // z = z0, dz = 1.0

    // Dual2_3d numbers for x, y, z
    Dual2_3d<double> x(x0, 1, 0, 0);
    Dual2_3d<double> y(y0, 0, 1, 0);
    Dual2_3d<double> z(z0, 0, 0, 1);

    // f(x, y, z) = sin(x * y + z^2)
    std::cout << "f(x, y, z) = sin(x * y + z^2)  ";
    // Dual1 version with respect to x
    Dual1<double> f1 = sin(x1 * y1 + pow(z1, 2));

    // Dual1 version with respect to y
    Dual1<double> f2 = sin(x2 * y2 + pow(z2, 2));

    // Dual1 version with respect to z
    Dual1<double> f3 = sin(x3 * y3 + pow(z3, 2));

    // Dual2_3d version
    Dual2_3d<double> f = sin(x * y + pow(z, 2));

    // Analytical derivatives
    double s = x0 * y0 + z0 * z0;
    double fx = cos(s) * y0;
    double fy = cos(s) * x0;
    double fz = cos(s) * 2 * z0;
    double fxx = -sin(s) * y0 * y0;
    double fyy = -sin(s) * x0 * x0;
    double fzz = -sin(s) * 4 * z0 * z0 + cos(s) * 2;
    double fxy = cos(s) - sin(s) * x0 * y0;
    double fxz = -sin(s) * y0 * 2 * z0;
    double fyz = -sin(s) * x0 * 2 * z0;

    auto almost_equal = [](double a, double b, double eps = 1e-9) {
        return std::abs(a - b) < eps;
        };

    assert(almost_equal(f.dx, fx));
    assert(almost_equal(f.dy, fy));
    assert(almost_equal(f.dz, fz));
    assert(almost_equal(f.dxx, fxx));
    assert(almost_equal(f.dyy, fyy));
    assert(almost_equal(f.dzz, fzz));
    assert(almost_equal(f.dxy, fxy));
    assert(almost_equal(f.dxz, fxz));
    assert(almost_equal(f.dyz, fyz));

    // Symmetry
    assert(almost_equal(f.dxy, f.dxy));
    assert(almost_equal(f.dxz, f.dxz));
    assert(almost_equal(f.dyz, f.dyz));

    // Dual1 vs Dual2_3d
    assert(almost_equal(f1.val, f.val) && almost_equal(f1.d1, f.dx));
    assert(almost_equal(f2.d1, f.dy));
    assert(almost_equal(f3.d1, f.dz));

    // Check chain rule
    double chain_rule_x = cos(s) * y0;
    double chain_rule_y = cos(s) * x0;
    double chain_rule_z = cos(s) * 2 * z0;
    assert(almost_equal(f1.d1, chain_rule_x));
    assert(almost_equal(f2.d1, chain_rule_y));
    assert(almost_equal(f3.d1, chain_rule_z));

    std::cout << "( good )\n";
}

//Test super hard function with all operators and functions
template<typename T>
T all_ops_func(const T& x, const T& y, const T& z) {
    // Use all arithmetic operators, exp, pow (all variants), log, sin, cos, tan
    T a = x + y - z;
    T b = x * y / (z + 1.0);
    T c = exp(a);
    T d = pow(b, 2.5);
    T e = pow(2.0, x - y);
    T f = pow(x + y + z, z);
    T g = log(x + y + z + 1.0);
    T h = sin(x) + cos(y) - tan(z);
    T result = ((a + b) * c - d / (e + 1.0)) + f * g + h;
    return result;
}

//Test super hard function with all operators and functions
static void test_all_combined() {
    // This function combines addition, multiplication, exponentiation, and more.
    for (double i = 0.1; i <= 5; i += 0.1) {
        for (double j = 0.1; j <= 5; j += 0.1) {
            for (double k = 0.1; k <= 5; k += 0.1) {
                // Create dual numbers
                // Dual1 numbers
                // Dual1 numbers for dx
                Dual1<double> x(i, 1.0);
                Dual1<double> y(j, 0.0);
                Dual1<double> z(k, 0.0);
                // Dual1 numbers for dy
                Dual1<double> x2(i, 0.0);
                Dual1<double> y2(j, 1.0);
                Dual1<double> z2(k, 0.0);
                // Dual1 numbers for dz
                Dual1<double> x3(i, 0.0);
                Dual1<double> y3(j, 0.0);
                Dual1<double> z3(k, 1.0);

                // Dual2_3d numbers
                Dual2_3d<double> x4(i, 1.0, 0.0, 0); // x = i, dx = 1.0
                Dual2_3d<double> y4(j, 0.0, 1.0, 0); // y = j, dx = 0.0, dy = 1.0, dz = 0.0
                Dual2_3d<double> z4(k, 0.0, 0.0, 1.0); // z = k, dx = 0.0, dy = 0.0, dz = 1.0

                // Function

                Dual1<double> f_1 = all_ops_func(x, y, z);
                Dual1<double> f_2 = all_ops_func(x2, y2, z2);
                Dual1<double> f_3 = all_ops_func(x3, y3, z3);
                Dual2_3d<double> f_4 = all_ops_func(x4, y4, z4);
                // Check values
                auto almost_equal = [](double a, double b, double rel_eps = 1e-9, double abs_eps = 1e-9) {
                    return std::abs(a - b) <= std::max(rel_eps * std::max(std::abs(a), std::abs(b)), abs_eps);
                    };
                if (!(almost_equal(f_4.val, f_1.val) &&
                    almost_equal(f_4.dx, f_1.d1) &&
                    almost_equal(f_4.dy, f_2.d1) &&
                    almost_equal(f_4.dz, f_3.d1))) {
                    std::cout << "Mismatch at x: " << i << ", y: " << j << ", z: " << k << "\n";
                    std::cout << "  Dual1: val = " << f_1.val << ", dx = " << f_1.d1 << ", dy = " << f_2.d1 << ", dz = " << f_3.d1 << "\n";
                    std::cout << "  Dual2_3d: val = " << f_4.val << ", dx = " << f_4.dx << ", dy = " << f_4.dy << ", dz = " << f_4.dz << "\n\n";
                }
            }
        }
    }
    std::cout << "( good )\n";
}

// Main function to run all tests
int main() {
    std::cout << "Testing Dual1 and Dual2_3d classes with various operators and functions:\n";

    std::cout << "\nTesting Dual1 class:\n";
    // Test the Dual1 class against operators
    test_dual1_operators();

    // Test the Dual1 class with complex numbers
    test_dual1_complex_operators();

    std::cout << "\nTesting Dual2_3d class:\n";
    // Test the Dual2_3d class with non-complex numbers
    test_dual2_3d_operators();

    // Test the Dual2_3d class with complex numbers
    test_dual2_3d_operators_complex();

    // Example for second-order Dual2_3d with higher order polynomial (yx^3 + x^2 + x + 1)
    std::cout << "\n\n";
    std::cout << "Testing both classes with funtions:\n";

    test_higher_order_polynomial_dual();

    // Example for exp function with Dual2_3d
    test_exp_func_dual();

    //Example for pow function with Dual2_3d
    test_pow_dual_scalar_func_dual();
    test_pow_scalar_dual_func_dual();
    test_pow_dual_dual_func_dual();

    // Example for log function with Dual2_3d
    test_log_func_dual();

    // Example for sqrt function with Dual2_3d
    test_sqrt_func_dual();

    // Test chain rule and symmetry properties
    std::cout << "\nTesting chain rule and symmetry properties:\n";
    test_chain_rule_and_symmetry();

    // Test all combined operators and functions
    std::cout << "\nTesting all operators and functions in one function  ";
    test_all_combined();

    std::cout << "\nAll tests completed.\n";

    return 0;
}
