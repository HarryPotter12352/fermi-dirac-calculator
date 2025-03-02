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