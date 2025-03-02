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


std::complex<double> gamma_integrand(std::complex<double> t, std::complex<double> z){
    if (std::abs(t) < 1e-10){
        std::cout << "Invalid value of t. Terminating computation";
        return std::complex<double>(0,0);
    }
    return std::pow(t,z-1.0)*std::exp(-t);
}

std::complex<double> gamma(std::complex<double> z, double lower_limit = 0, double upper_limit = 1e5, double steps = 1e-5){
    double step_size = (upper_limit - lower_limit)/steps;
    std::complex<double> integral(0,0);
    integral += 0.5 *(gamma_integrand(lower_limit, z) + gamma_integrand(upper_limit, z));

    for (int i = 1; i < steps; i++){
        double t = lower_limit + i * step_size;
        integral+= gamma_integrand(t,z);
    }
    return integral;
}