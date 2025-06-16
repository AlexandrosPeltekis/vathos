#include <iostream>
#include <cmath>
#include <complex>
#include <array>
#include <vector>
#include "dual1_class.h"
#include "dual2_3d_class.h"

template<typename T>
std::array<T, 3> compute_k(const std::array<double, 3>& O, const std::array<T, 3>& H, double sign) {
    std::array<T, 3> k;
    std::array<T, 3> cp = {
        O[1] * H[2] - O[2] * H[1],
        O[2] * H[0] - O[0] * H[2],
        O[0] * H[1] - O[1] * H[0]
    };
    for (int i = 0; i < 3; ++i)
        k[i] = sign * 0.5 * cp[i];
    return k;
}

template<typename T>
using Vec3 = std::array<T, 3>;

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
        Dual1<std::complex<double> > result13 = a + c1;
        Dual1<std::complex<double> > result14 = a - c1;
        Dual1<std::complex<double> > result15 = a * c1;
        Dual1<std::complex<double> > result16 = a / c1;
        Dual1<std::complex<double> > result17 = c1 + a;
        Dual1<std::complex<double> > result18 = c1 - a;
        Dual1<std::complex<double> > result19 = c1 * a;
        Dual1<std::complex<double> > result20 = c1 / a;

        std::cout << "Completed Dual1 non-complex operations successfully.\n";
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
        Dual1<std::complex<double> > result1 = a + b;
        Dual1<std::complex<double> > result2 = a - b;
        Dual1<std::complex<double> > result3 = a * b;
        Dual1<std::complex<double> > result4 = a / b;
        Dual1<std::complex<double> > result5 = a + 2.0;
        Dual1<std::complex<double> > result6 = a - 2.0;
        Dual1<std::complex<double> > result7 = a * 2.0;
        Dual1<std::complex<double> > result8 = a / 2.0;
        Dual1<std::complex<double> > result9 = 2.0 + a;
        Dual1<std::complex<double> > result10 = 2.0 - a;
        Dual1<std::complex<double> > result11 = 2.0 * a;
        Dual1<std::complex<double> > result12 = 2.0 / a;
        Dual1<std::complex<double> > result13 = a + c1;
        Dual1<std::complex<double> > result14 = a - c1;
        Dual1<std::complex<double> > result15 = a * c1;
        Dual1<std::complex<double> > result16 = a / c1;
        Dual1<std::complex<double> > result17 = c1 + a;
        Dual1<std::complex<double> > result18 = c1 - a;
        Dual1<std::complex<double> > result19 = c1 * a;
        Dual1<std::complex<double> > result20 = c1 / a;

        std::cout << "Completed Dual1 complex operations successfully.\n";
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
        Dual2_3d<std::complex<double> > result13 = a + c1;
        Dual2_3d<std::complex<double> > result14 = a - c1;
        Dual2_3d<std::complex<double> > result15 = a * c1;
        Dual2_3d<std::complex<double> > result16 = a / c1;
        Dual2_3d<std::complex<double> > result17 = c1 + a;
        Dual2_3d<std::complex<double> > result18 = c1 - a;
        Dual2_3d<std::complex<double> > result19 = c1 * a;
        Dual2_3d<std::complex<double> > result20 = c1 / a;

        std::cout << "Completed Dual2_3d operations successfully.\n";
    }
}

static void test_dual2_3d_operators_complex(bool verbose = false) {
    Dual2_3d<std::complex<double> > a(std::complex<double>(2.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(0.5, 0.0), std::complex<double>(0.0, 0.0));
    Dual2_3d<std::complex<double> > b(std::complex<double>(3.0, 0.0), std::complex<double>(0.5, 0.0), std::complex<double>(0.25, 0.0), std::complex<double>(1.0, 0.0));
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
        Dual2_3d<std::complex<double> > result1 = a + b;
        Dual2_3d<std::complex<double> > result2 = a - b;
        Dual2_3d<std::complex<double> > result3 = a * b;
        Dual2_3d<std::complex<double> > result4 = a / b;
        Dual2_3d<std::complex<double> > result5 = a + 2.0;
        Dual2_3d<std::complex<double> > result6 = a - 2.0;
        Dual2_3d<std::complex<double> > result7 = a * 2.0;
        Dual2_3d<std::complex<double> > result8 = a / 2.0;
        Dual2_3d<std::complex<double> > result9 = 2.0 + a;
        Dual2_3d<std::complex<double> > result10 = 2.0 - a;
        Dual2_3d<std::complex<double> > result11 = 2.0 * a;
        Dual2_3d<std::complex<double> > result12 = 2.0 / a;
        Dual2_3d<std::complex<double> > result13 = a + c1;
        Dual2_3d<std::complex<double> > result14 = a - c1;
        Dual2_3d<std::complex<double> > result15 = a * c1;
        Dual2_3d<std::complex<double> > result16 = a / c1;
        Dual2_3d<std::complex<double> > result17 = c1 + a;
        Dual2_3d<std::complex<double> > result18 = c1 - a;
        Dual2_3d<std::complex<double> > result19 = c1 * a;
        Dual2_3d<std::complex<double> > result20 = c1 / a;
        std::cout << "Completed Dual2_3D complex operations successfully.\n";
    }
}

static void test_higher_order_polynomial_dual2_3d () {
    std::cout << "Second-order Dual2_3d example:\n";
    std::cout << "x: f(x) = yx^3 + x^2 + x + 1\n";
    bool all_match = true;
    for (double i = 0; i <= 5; i += 0.1) {
        for (double k = 0; k <= 5; k += 0.1) {
            // Analytical solution: f(x) = yx^3 + x^2 + x + 1
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
                std::cout << "  Expected f(x) = " << f_x << ", df/dx = " << f_dx << ", df/dy = " << f_dy << ", d^2f/dx^2 = " << f_dxx << "\n";
                std::cout << "           d^2f/dy^2 = " << f_dyy << ", d^2f/dz^2 = " << f_dzz << ", d^2f/dxy = " << f_dxy << ", d^2f/dxz = " << f_dxz << ", d^2f/dyz = " << f_dyz << "\n";
                std::cout << "        Got f(x) =      " << x3.val << ", df/dx = " << x3.dx << ", df/dy = " << x3.dy << ", d^2f/dx^2 = " << x3.dxx << "\n";
                std::cout << "           d^2f/dy^2 = " << x3.dyy << ", d^2f/dz^2 = " << x3.dzz << ", d^2f/dxy = " << x3.dxy << ", d^2f/dxz = " << x3.dxz << ", d^2f/dyz = " << x3.dyz << "\n\n";
            }
        }
    }
    if (all_match) {
        std::cout << "All values match the analytical solutions!!\n";
    }
}

int main() {
	// Example for second-order Dual2_3d with higher order polynomial (yx^3 + x^2 + x + 1)
    test_higher_order_polynomial_dual2_3d();

    // Test the Dual1 class against operators
    std::cout << "\nTesting non-complex Dual1 operators:\n";
    test_dual1_operators();

    // Test the Dual1 class with complex numbers
    std::cout << "\nTesting complex Dual1 operators:\n";
    test_dual1_complex_operators();

    // Test the Dual2_3d class with non-complex numbers
    std::cout << "\nTesting non-complex Dual2_3d operators:\n";
    test_dual2_3d_operators();
    
    // Test the Dual2_3d class with complex numbers
    std::cout << "\nTesting complex Dual2_3d operators:\n";
    test_dual2_3d_operators_complex();

	return 0;
}
