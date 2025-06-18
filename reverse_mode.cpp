#include <iostream>
#include <cmath>
#include <memory>
#include <vector>
#include <functional>
#include "reverseVar.h"

int main() {
    rVar x(1.0);
    rVar y(2);

    //rVar f = exp(x * y + sin(x)) + pow(x, 2.0) / y;
	//rVar f = sin(x * y) + x * cos(y) + 2 * y;
	//rVar f = pow(2, x*y);
	rVar f = x * y;

    f.backward();

    std::cout << "f = " << f.val() << "\n";
	std::cout << "df/dx = " << x.grad() << "\n"; //  
	std::cout << "df/dy = " << y.grad() << "\n"; // 
	return 0;									    
}
