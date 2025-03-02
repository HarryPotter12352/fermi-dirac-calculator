# fermi-dirac-calculator

A script written in C++ to compute the various Fermi-Dirac quantities for a given species.

While writing these scripts, I try to use my own functions as much as I can, all the way to basic numerical integration. I resort to using the MATLAB engine for symbolic integration.

# Definitions of functions

## Polylogarithm

The polylogarithm is defined by a power series in a complex parameter $z$, given by

$$
Li_s(z) = \sum_{n=1}^{\infty}\frac{z^n}{n^s}
$$

The polylogarithm is defined for all complex orders $s$ and all complex arguments where $|z| < 1$.

## Gamma function

The gamma function is a continuation of the factorial function to complex numbers. For any complex number $z$ with $\mathcal{R}(z) > 0$, the gamma function is defined as

$$
\Gamma(z) = \int_{0}^{\infty}t^{z-1}e^{-t}dt
$$
