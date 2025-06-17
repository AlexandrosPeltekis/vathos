#include <iostream>
#include <cmath>
#include <complex>
#include <array>
#include <vector>
#include "dual1_class.h"
#include "dual2_3d_class.h"

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

static void test_higher_order_polynomial_dual2_3d () {
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

            Dual2_3d<double> x(i, 1.0, 0.0, 0); // x = i, dx = 1.0, dy = 0.0, dz = 0.0
            Dual2_3d<double> y(k, 0.0, 1.0, 0); // y = k, dx = 0.0, dy = 1.0, dz = 0.0

            //Funciton
            Dual2_3d<double> x3 = y * pow(x, 3) + pow(x, 2) + x + 1.0;

            auto almost_equal = [](double a, double b, double eps = 1e-9) {
                return std::abs(a - b) < eps;
                };
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

        }
    }
    if (all_match) {
        std::cout << "  ( good )\n";
    }
}

static void test_exp_func_dual2_3d()  {
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
                double f_dxx = j * k * j * k * exp(i * j * k); // Second derivative with respect to x
                double f_dyy = i * k * i * k * exp(i * j * k); // Second derivative with respect to y
                double f_dzz = i * j * i * j * exp(i * j * k); // Second derivative with respect to z
                double f_dxy = k * exp(i * j * k); // Mixed derivative with respect to x and y
                double f_dxz = j * exp(i * j * k); // Mixed derivative with respect to x and z
                double f_dyz = i * exp(i * j * k); // Mixed derivative with respect to y and z
                Dual2_3d<double> x(i, 1.0, 0.0, 0); // x = i, dx = 1.0, dy = 0.0, dz = 0.0
                Dual2_3d<double> y(j, 0.0, 1.0, 0); // y = j, dx = 0.0, dy = 1.0, dz = 0.0
                Dual2_3d<double> z(k, 0.0, 0.0, 1.0); // z = k, dx = 0.0, dy = 0.0, dz = 1.0

                // Function
                Dual2_3d<double> x3 = exp(x * y * z);
                auto almost_equal = [](double a, double b, double rel_eps = 1e-9, double abs_eps = 1e-9) {
                    return std::abs(a - b) <= std::max(rel_eps * std::max(std::abs(a), std::abs(b)), abs_eps);
                    };
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
                    std::cout << "x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Mismatch in values!\n";
                    if (!almost_equal(x3.val, f_x))
                        std::cout << "  f(x,y,z): expected " << f_x << ", actual " << x3.val << "\n";
                    if (!almost_equal(x3.dx, f_dx))
                        std::cout << "  df/dx: expected " << f_dx << ", actual " << x3.dx << "\n";
                    if (!almost_equal(x3.dy, f_dy))
                        std::cout << "  df/dy: expected " << f_dy << ", actual " << x3.dy << "\n";
                    if (!almost_equal(x3.dz, f_dz))
                        std::cout << "  df/dz: expected " << f_dz << ", actual " << x3.dz << "\n";
                    if (!almost_equal(x3.dxx, f_dxx))
                        std::cout << "  d^2f/dx^2: expected " << f_dxx << ", actual " << x3.dxx << "\n";
                    if (!almost_equal(x3.dyy, f_dyy))
                        std::cout << "  d^2f/dy^2: expected " << f_dyy << ", actual " << x3.dyy << "\n";
                    if (!almost_equal(x3.dzz, f_dzz))
                        std::cout << "  d^2f/dz^2: expected " << f_dzz << ", actual " << x3.dzz << "\n";
                    if (!almost_equal(x3.dxy, f_dxy))
                        std::cout << "  d^2f/dxy: expected " << f_dxy << ", actual " << x3.dxy << "\n";
                    if (!almost_equal(x3.dxz, f_dxz))
                        std::cout << "  d^2f/dxz: expected " << f_dxz << ", actual " << x3.dxz << "\n";
                    if (!almost_equal(x3.dyz, f_dyz))
                        std::cout << "  d^2f/dyz: expected " << f_dyz << ", actual " << x3.dyz << "\n";
                    std::cout << "\n";
                }
                break;
            }
        }
    }
	if (all_match) {
        std::cout << "  ( good )\n";
    }
}

static void test_pow_dual_scalar_func_dual2_3d() {
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
                Dual2_3d<double> x(i, 1.0, 0.0, 0); // x = i, dx = 1.0
                Dual2_3d<double> y(j, 0.0, 1.0, 0); // y = j, dx = 0.0, dy = 1.0, dz = 0.0
                Dual2_3d<double> z(k, 0.0, 0.0, 1.0); // z = k, dx = 0.0, dy = 0.0, dz = 1.0
                // Function
                Dual2_3d<double> x3 = pow(x * y * z, 3);
                auto almost_equal = [](double a, double b, double rel_eps = 1e-9, double abs_eps = 1e-9) {
                    return std::abs(a - b) <= std::max(rel_eps * std::max(std::abs(a), std::abs(b)), abs_eps);
                    };
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
                    std::cout << "x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Mismatch in values!\n";
                    if (!almost_equal(x3.val, f_x))
                        std::cout << "  f(x,y,z): expected " << f_x << ", actual " << x3.val << "\n";
                    if (!almost_equal(x3.dx, f_dx))
                        std::cout << "  df/dx: expected " << f_dx << ", actual " << x3.dx << "\n";
                    if (!almost_equal(x3.dy, f_dy))
                        std::cout << "  df/dy: expected " << f_dy << ", actual " << x3.dy << "\n";
                    if (!almost_equal(x3.dz, f_dz))
                        std::cout << "  df/dz: expected " << f_dz << ", actual " << x3.dz << "\n";
                    if (!almost_equal(x3.dxx, f_dxx))
                        std::cout << "  d^2f/dx^2: expected " << f_dxx << ", actual " << x3.dxx << "\n";
                    if (!almost_equal(x3.dyy, f_dyy))
                        std::cout << "  d^2f/dy^2: expected " << f_dyy << ", actual " << x3.dyy << "\n";
                    if (!almost_equal(x3.dzz, f_dzz))
                        std::cout << "  d^2f/dz^2: expected " << f_dzz << ", actual " << x3.dzz << "\n";
                    if (!almost_equal(x3.dxy, f_dxy))
                        std::cout << "  d^2f/dxy: expected " << f_dxy << ", actual " << x3.dxy << "\n";
                    if (!almost_equal(x3.dxz, f_dxz))
                        std::cout << "  d^2f/dxz: expected " << f_dxz << ", actual " << x3.dxz << "\n";
                    if (!almost_equal(x3.dyz, f_dyz))
                        std::cout << "  d^2f/dyz: expected " << f_dyz << ", actual " << x3.dyz << "\n";
                    std::cout << "\n";
                }
            }
        }
    }
    if (all_match) {
        std::cout << "  ( good )\n";
	}
}

static void test_pow_scalar_dual_func_dual2_3d() {
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
				double f_dxy = log(3) * k *(log(3) * i * j * k + 1) * pow(3, i * j * k); // Mixed derivative with respect to x and y is log(3)^2 * 3^(xyz) * yz * xz
				double f_dxz = log(3) * j * (log(3) * i * j * k + 1) * pow(3, i * j * k); // Mixed derivative with respect to x and z is log(3)^2 * 3^(xyz) * yz * xy
				double f_dyz = log(3) * i * (log(3) * i * j * k + 1) * pow(3, i * j * k); // Mixed derivative with respect to y and z is log(3)^2 * 3^(xyz) * xz * xy

                Dual2_3d<double> x(i, 1.0, 0.0, 0); // x = i, dx = 1.0
                Dual2_3d<double> y(j, 0.0, 1.0, 0); // y = j, dx = 0.0, dy = 1.0, dz = 0.0
                Dual2_3d<double> z(k, 0.0, 0.0, 1.0); // z = k, dx = 0.0, dy = 0.0, dz = 1.0
                // Function
                Dual2_3d<double> x3 = pow(3, x * y * z);
                auto almost_equal = [](double a, double b, double rel_eps = 1e-9, double abs_eps = 1e-9) {
                    return std::abs(a - b) <= std::max(rel_eps * std::max(std::abs(a), std::abs(b)), abs_eps);
                    };
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
                    std::cout << "\nx: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Mismatch in values!\n";
                    if (!almost_equal(x3.val, f_x))
                        std::cout << "  f(x,y,z): expected " << f_x << ", actual " << x3.val << "\n";
                    if (!almost_equal(x3.dx, f_dx))
                        std::cout << "  df/dx: expected " << f_dx << ", actual " << x3.dx << "\n";
					if (!almost_equal(x3.dy, f_dy))
						std::cout << "  df/dy: expected " << f_dy << ", actual " << x3.dy << "\n";
					if (!almost_equal(x3.dz, f_dz))
						std::cout << "  df/dz: expected " << f_dz << ", actual " << x3.dz << "\n";
					if (!almost_equal(x3.dxx, f_dxx))
						std::cout << "  d^2f/dx^2: expected " << f_dxx << ", actual " << x3.dxx << "\n";
					if (!almost_equal(x3.dyy, f_dyy))
						std::cout << "  d^2f/dy^2: expected " << f_dyy << ", actual " << x3.dyy << "\n";
					if (!almost_equal(x3.dzz, f_dzz))
						std::cout << "  d^2f/dz^2: expected " << f_dzz << ", actual " << x3.dzz << "\n";
					if (!almost_equal(x3.dxy, f_dxy))
						std::cout << "  d^2f/dxy: expected " << f_dxy << ", actual " << x3.dxy << "\n";
					if (!almost_equal(x3.dxz, f_dxz))
						std::cout << "  d^2f/dxz: expected " << f_dxz << ", actual " << x3.dxz << "\n";
                    if (!almost_equal(x3.dyz, f_dyz))
						std::cout << "  d^2f/dyz: expected " << f_dyz << ", actual " << x3.dyz << "\n";
                    std::cout << "\n";
                }
			}
        }
        
    }
    if (all_match) {
		std::cout << " ( good )\n";
    }
}

static void test_pow_dual_dual_func_dual2_3d() {
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
                double f_dxx = f_x * ((log_s + (i / s)) * (log_s + (i / s)) + ((2.0 / s) - (i / (s * s) )) );
                double f_dyy = i * (i - 1) * pow(s, i - 2);
                double f_dzz = i * (i - 1) * pow(s, i - 2);
                double f_dxy = i * pow(s, i - 1) * (log_s + (i - 1) / s) + pow(s, i - 1);
                double f_dxz = i * pow(s, i - 1) * (log_s + (i - 1) / s) + pow(s, i - 1);
                double f_dyz = i * (i - 1) * pow(s, i - 2);
                Dual2_3d<double> x(i, 1.0, 0.0, 0); // x = i, dx = 1.0
                Dual2_3d<double> y(j, 0.0, 1.0, 0); // y = j, dx = 0.0, dy = 1.0, dz = 0.0
                Dual2_3d<double> z(k, 0.0, 0.0, 1.0); // z = k, dx = 0.0, dy = 0.0, dz = 1.0

                // Function
                Dual2_3d<double> x3 = pow(x + y + z, x);
				Dual2_3d<double> autoCorrect = exp(x * log(x + y + z)); // Using the property a^b = exp(b * log(a))

                auto almost_equal = [](double a, double b, double rel_eps = 1e-9, double abs_eps = 1e-9) {
                    return std::abs(a - b) <= std::max(rel_eps * std::max(std::abs(a), std::abs(b)), abs_eps);
                    };
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
                    std::cout << "x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Mismatch in values!\n";
                    if (!almost_equal(x3.val, f_x))
						std::cout << "  f(x,y,z): expected " << f_x << ", actual " << x3.val << "\n";
                    if (!almost_equal(x3.dx, f_dx))
						std::cout << "  df/dx: expected " << f_dx << ", actual " << x3.dx << "\n";
					if (!almost_equal(x3.dy, f_dy))
						std::cout << "  df/dy: expected " << f_dy << ", actual " << x3.dy << "\n";
					if (!almost_equal(x3.dz, f_dz))
						std::cout << "  df/dz: expected " << f_dz << ", actual " << x3.dz << "\n";
					if (!almost_equal(x3.dxx, f_dxx))
						std::cout << "  d^2f/dx^2: expected " << f_dxx << ", actual " << x3.dxx << "\n";
					if (!almost_equal(x3.dyy, f_dyy))
						std::cout << "  d^2f/dy^2: expected " << f_dyy << ", actual " << x3.dyy << "\n";
					if (!almost_equal(x3.dzz, f_dzz))
						std::cout << "  d^2f/dz^2: expected " << f_dzz << ", actual " << x3.dzz << "\n";
					if (!almost_equal(x3.dxy, f_dxy))
						std::cout << "  d^2f/dxy: expected " << f_dxy << ", actual " << x3.dxy << "\n";
					if (!almost_equal(x3.dxz, f_dxz))
						std::cout << "  d^2f/dxz: expected " << f_dxz << ", actual " << x3.dxz << "\n";
					if (!almost_equal(x3.dyz, f_dyz))
                        std::cout << "  d^2f/dyz: expected " << f_dyz << ", actual " << x3.dyz << "\n";
					std::cout << "\n";
					
                }
            }
        }
    }
    if (all_match)
        std::cout << "  ( good )\n";
}

static void test_log_func_dual2_3d() {
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
                Dual2_3d<double> x(i, 1.0, 0.0, 0); // x = i, dx = 1.0
                Dual2_3d<double> y(j, 0.0, 1.0, 0); // y = j, dx = 0.0, dy = 1.0, dz = 0.0
                Dual2_3d<double> z(k, 0.0, 0.0, 1.0); // z

                // Function
                Dual2_3d<double> x3 = log(x + y + z);
                auto almost_equal = [](double a, double b, double rel_eps = 1e-9, double abs_eps = 1e-9) {
                    return std::abs(a - b) <= std::max(rel_eps * std::max(std::abs(a), std::abs(b)), abs_eps);
                    };
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
                    std::cout << "x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Mismatch in values!\n";
                    if (!almost_equal(x3.val, f_x))
                        std::cout << "  f(x,y,z): expected " << f_x << ", actual " << x3.val << "\n";
                    if (!almost_equal(x3.dx, f_dx))
                        std::cout << "  df/dx: expected " << f_dx << ", actual " << x3.dx << "\n";
                    if (!almost_equal(x3.dy, f_dy))
                        std::cout << "  df/dy: expected " << f_dy << ", actual " << x3.dy << "\n";
                    if (!almost_equal(x3.dz, f_dz))
                        std::cout << "  df/dz: expected " << f_dz << ", actual " << x3.dz << "\n";
                    if (!almost_equal(x3.dxx, f_dxx))
                        std::cout << "  d^2f/dx^2: expected " << f_dxx << ", actual " << x3.dxx << "\n";
                    if (!almost_equal(x3.dyy, f_dyy))
                        std::cout << "  d^2f/dy^2: expected " << f_dyy << ", actual " << x3.dyy << "\n";
                    if (!almost_equal(x3.dzz, f_dzz))
                        std::cout << "  d^2f/dz^2: expected " << f_dzz << ", actual " << x3.dzz << "\n";
                    if (!almost_equal(x3.dxy, f_dxy))
                        std::cout << "  d^2f/dxy: expected " << f_dxy << ", actual " << x3.dxy << "\n";
                    if (!almost_equal(x3.dxz, f_dxz))
                        std::cout << "  d^2f/dxz: expected " << f_dxz << ", actual " << x3.dxz << "\n";
                    if (!almost_equal(x3.dyz, f_dyz))
                        std::cout << "  d^2f/dyz: expected " << f_dyz << ", actual " << x3.dyz << "\n";
                    std::cout << "\n";
                }
            }
        }
    }
    if (all_match) {
        std::cout << "  ( good )\n";
    }
}

static void test_sqrt_func_dual2_3d() {
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
                Dual2_3d<double> x(i, 1.0, 0.0, 0); // x = i, dx = 1.0
                Dual2_3d<double> y(j, 0.0, 1.0, 0); // y = j, dx = 0.0, dy = 1.0, dz = 0.0

                Dual2_3d<double> z(k, 0.0, 0.0, 1.0); // z = k, dx = 0.0, dy = 0.0, dz = 1.0
                // Function
                Dual2_3d<double> x3 = sqrt(x + y + z);
                auto almost_equal = [](double a, double b, double rel_eps = 1e-9, double abs_eps = 1e-9) {
                    return std::abs(a - b) <= std::max(rel_eps * std::max(std::abs(a), std::abs(b)), abs_eps);
                    };
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
                    std::cout << "x: " << x.val << ", y: " << y.val << ", z: " << z.val << "\n";
                    std::cout << "  Mismatch in values!\n";
                    if (!almost_equal(x3.val, f_x))
                        std::cout << "  f(x,y,z): expected " << f_x << ", actual " << x3.val << "\n";
                    if (!almost_equal(x3.dx, f_dx))
                        std::cout << "  df/dx: expected " << f_dx << ", actual " << x3.dx << "\n";
                    if (!almost_equal(x3.dy, f_dy))
                        std::cout << "  df/dy: expected " << f_dy << ", actual " << x3.dy << "\n";
                    if (!almost_equal(x3.dz, f_dz))
                        std::cout << "  df/dz: expected " << f_dz << ", actual " << x3.dz << "\n";
                    if (!almost_equal(x3.dxx, f_dxx))
                        std::cout << "  d^2f/dx^2: expected " << f_dxx << ", actual " << x3.dxx << "\n";
                    if (!almost_equal(x3.dyy, f_dyy))
                        std::cout << "  d^2f/dy^2: expected " << f_dyy << ", actual " << x3.dyy << "\n";
                    if (!almost_equal(x3.dzz, f_dzz))
                        std::cout << "  d^2f/dz^2: expected " << f_dzz << ", actual " << x3.dzz << "\n";
                    if (!almost_equal(x3.dxy, f_dxy))
                        std::cout << "  d^2f/dxy: expected " << f_dxy << ", actual " << x3.dxy << "\n";
                    if (!almost_equal(x3.dxz, f_dxz))
                        std::cout << "  d^2f/dxz: expected " << f_dxz << ", actual " << x3.dxz << "\n";
                    if (!almost_equal(x3.dyz, f_dyz))
                        std::cout << "  d^2f/dyz: expected " << f_dyz << ", actual " << x3.dyz << "\n";
                    std::cout << "\n";
                }
            }
        }
    }
    if (all_match) {
        std::cout << "  ( good )\n";
    }
}

int main() {
	std::cout << "Testing Dual1 and Dual2_3d classes with various operators and functions:\n";
	
    // Test the Dual1 class against operators
    test_dual1_operators();

    // Test the Dual1 class with complex numbers
    test_dual1_complex_operators();

    // Test the Dual2_3d class with non-complex numbers
    test_dual2_3d_operators();

    // Test the Dual2_3d class with complex numbers
    test_dual2_3d_operators_complex();
    
    // Example for second-order Dual2_3d with higher order polynomial (yx^3 + x^2 + x + 1)
    std::cout << "\n\n";
    test_higher_order_polynomial_dual2_3d();

	// Example for exp function with Dual2_3d
	test_exp_func_dual2_3d();

	//Example for pow function with Dual2_3d
	test_pow_dual_scalar_func_dual2_3d();
    test_pow_scalar_dual_func_dual2_3d();
	test_pow_dual_dual_func_dual2_3d();

	// Example for log function with Dual2_3d
    test_log_func_dual2_3d();  

	// Example for sqrt function with Dual2_3d
    test_sqrt_func_dual2_3d();
	std::cout << "\nAll tests completed.\n";

	return 0;
}
