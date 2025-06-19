#include <iostream>
#include <cmath>
#include <memory>
#include <vector>
#include <functional>
#include "reverseVar.h"
#include "dual1_class.h"
#include "dual2_3d_class.h"


// The function to differentiate
template<typename Vec>
auto my_func(const Vec& vars) -> decltype(vars[0] * vars[0] * vars[1] + sin(vars[0])) {
    return vars[0] * vars[0] * vars[1] + sin(vars[0]);
}

// Another function to differentiate with many fucntions like log, trig,pow, exp etc.
template<typename Vec>
auto my_func2(const Vec& vars) -> decltype(log(vars[0]) + pow(vars[1], 2.0) + exp(vars[0] * vars[1])) {
    return log(vars[0]) + pow(vars[1], 2.0) + exp(vars[0] * vars[1]);
}



// Function to compute the Hessian using reverse-on-forward AD
template<typename T, typename Func>
void compute_hessian_with_reverse_on_forward(const std::vector<T>& x0, std::vector<std::vector<T>>& hessian, Func func) {
    for (auto& row : hessian)
        std::fill(row.begin(), row.end(), 0.0);
    int N = x0.size();
    for (int i = 0; i < N; ++i) {
        std::vector<rVar<T>> x_rev(N);
        for (int j = 0; j < N; ++j)
            x_rev[j] = rVar<T>(x0[j]);
        auto grad_i = [&](const std::vector<rVar<T>>& x) {
            std::vector<Dual1<rVar<T>>> x_fwd(N);
            for (int k = 0; k < N; ++k)
                x_fwd[k] = Dual1<rVar<T>>(x[k], (k == i) ? rVar<T>(1.0) : rVar<T>(0.0));
            return func(x_fwd).d1;
            };
        rVar<T> g = grad_i(x_rev);
        g.backward();
        for (int j = 0; j < N; ++j)
            hessian[j][i] = x_rev[j].grad();
        for (int j = 0; j < N; ++j)
            x_rev[j].zero_grad();
    }
}

int main() {
    constexpr int N = 2;
    std::vector<double> x0 = { 1.0, 2.0 };
    std::vector<std::vector<double>> hessian(N, std::vector<double>(N, 0.0));
    
	// Compute Hessian using reverse-on-forward AD
	std::cout << "Computing Hessian for a function...\n";
    compute_hessian_with_reverse_on_forward<double>(x0, hessian, my_func<std::vector<Dual1<rVar<double>>>>);

    // Print Hessian
	std::cout << "Hessian from reverse-on-forward AD:\n";
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j)
            std::cout << hessian[i][j] << " ";
        std::cout << "\n";
    }

	//Dual2_3d to compare to hessian
	Dual2_3d<double> x(1.0, 1.0, 0.0, 0.0, 0.0); // x = 1.0, dx/dx = 1.0
	Dual2_3d<double> y(2.0, 0.0, 1.0, 0.0, 0.0); // y = 2.0, dy/dx = 1.0

	Dual2_3d<double> f = my_func<std::vector<Dual2_3d<double>>>({ x, y });

	//print the hessian from Dual2_3d
	std::cout << "\nHessian from Forward mode (Dual2_3d):\n";
    std::vector<std::vector<double>> hessian_dual2 = {
        {f.dxx, f.dxy},
        {f.dxy, f.dyy}
    };
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j)
            std::cout << hessian_dual2[i][j] << " ";
        std::cout << "\n";
	}

	// Try another function
	std::cout << "\n\nComputing Hessian for another function...\n";
    std::vector<std::vector<double>> hessian2(N, std::vector<double>(N, 0.0));
    compute_hessian_with_reverse_on_forward<double>(x0, hessian2, my_func2<std::vector<Dual1<rVar<double>>>>);
    // Print Hessian for the second function
    std::cout << "Hessian from reverse-on-forward AD for second function:\n";
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j)
            std::cout << hessian2[i][j] << " ";
        std::cout << "\n";
	}

    // Dual2_3d example for the second function
    Dual2_3d<double> x2(1.0, 1.0, 0.0, 0.0, 0.0); // x = 1.0, dx/dx = 1.0
    Dual2_3d<double> y2(2.0, 0.0, 1.0, 0.0, 0.0); // y = 2.0, dy/dx = 1.0
    Dual2_3d<double> f2 = my_func2<std::vector<Dual2_3d<double>>>({ x2, y2 });
    //print the hessian from Dual2_3d for the second function
    std::cout << "\nHessian from Forward mode (Dual2_3d) for second function:\n";
    std::vector<std::vector<double>> hessian_dual2_2 = {
        {f2.dxx, f2.dxy},
        {f2.dxy, f2.dyy}
    };
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j)
            std::cout << hessian_dual2_2[i][j] << " ";
        std::cout << "\n";
	}
    
	return 0;
}

// Example usage of reverse mode AD with rVar
//int main() {
 //   rVar x(1.0);
 //   rVar y(2.0);
 //   //rVar f = exp(x * y + sin(x)) + pow(x, 2.0) / y;
	////rVar f = sin(x * y) + x * cos(y) + 2 * y;
	////rVar f = pow(2, x*y);
	//rVar f = x * y;
 //   f.backward();
 //   std::cout << "f = " << f.val() << "\n";
	//std::cout << "df/dx = " << x.grad() << "\n"; //  
	//std::cout << "df/dy = " << y.grad() << "\n"; // 
	//return 0;									    
//}
