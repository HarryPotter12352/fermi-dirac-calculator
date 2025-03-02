#include <iostream>
#include <complex>

std::complex<double> polylog(std::complex<double> s, std::complex<double> z){
    std::complex<double> result(0,0);
    int n = 1;
    while (n < 1000){
        std::complex<double> term(0,0);
        term = pow(z,n)/pow(n,s);
        if(std::isnan(real(term))||std::isnan(imag(term))||std::isnan(real(result))||std::isnan(imag(result))){
            std::cout << "Encountered NaN during computation of polylogarithm. Terminating computation of polylogarithm" << std::endl;
            break;
        }
        else if (std::abs(term) < 1e-10){
            std::cout << "Encountred convergence barrier during computation of polylogarithm. Terminating computation of polylogarithm";
        }
        else{
            result+=term;
            n++;
        }
        return result;
    }
}

std::complex<double> gamma_integrand(std::complex<double> t, std::complex<double> z) {
    if (std::abs(t) < 1e-10) {
        std::cout << "Invalid value of t. Terminating computation" << std::endl;
        return std::complex<double>(0, 0);
    }
    return std::pow(t, z - 1.0) * std::exp(-t);
}

std::complex<double> gamma(std::complex<double> z, double lower_limit = 0, double upper_limit = 1e5, double steps = 1e5) {
    int num_steps = static_cast<int>(steps);
    double step_size = (upper_limit - lower_limit) / num_steps;
    std::complex<double> integral(0, 0);
    integral += 0.5 * (gamma_integrand(std::complex<double>(lower_limit, 0), z) + gamma_integrand(std::complex<double>(upper_limit, 0), z));

    for (int i = 1; i < num_steps; i++) {
        double t = lower_limit + i * step_size;
        integral += gamma_integrand(std::complex<double>(t, 0), z);
    }

    integral *= step_size;
    
    return integral;
}


#include <iostream>
#include <complex>
#include <cmath>

std::complex<double> digamma_integrand(std::complex<double> x, std::complex<double> z) {
    if (std::abs(std::real(x) - 1.0) < 1e-10 || std::abs(std::real(x)) < 1e-10) {
        std::cout << "Invalid value of x. Skipping computation for x = " << x << std::endl;
        return std::complex<double>(0, 0);
    }
    return (1.0 - std::pow(x, z - 1.0)) / (1.0 - x);
}

double euler_mascheroni = 0.57721; // Euler-Mascheroni constant (corrected value)

std::complex<double> digamma(std::complex<double> z, double steps = 1e6) {
    double lower_limit = 0;
    double upper_limit = 1;
    std::complex<double> integral(0, 0);
    int num_steps = static_cast<int>(steps);
    double step_size = (upper_limit - lower_limit) / num_steps;
    
    integral += 0.5 * (digamma_integrand(std::complex<double>(lower_limit, 0), z) + digamma_integrand(std::complex<double>(upper_limit, 0), z));

    for (int i = 1; i < num_steps; i++) {
        double x_real = lower_limit + i * step_size;
        std::complex<double> x(x_real, 0);
        integral += digamma_integrand(x, z);
    }
    
    integral *= step_size;
    return integral - std::complex<double>(euler_mascheroni, 0);
}
