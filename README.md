# fermi-dirac-calculator

A script written in C++ to compute the various Fermi-Dirac quantities for a given species.

While writing these scripts, I try to use my own functions as much as I can, all the way to basic numerical integration. I resort to using the MATLAB engine for symbolic integration.
I operate almost solely in the complex domain to ensure compliance with the initial definitions of the functions.

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

## Digamma function

The digamma function appears multiple times in the relativistic expansion of Fermi-Dirac integrals. It is commonly expressed as

$$
\psi(z) = \frac{d}{dz}\ln \Gamma(z) = \frac{\Gamma'(z)}{\Gamma(z)}
$$

I use the integral representation to allow for easier computation.

$$
\psi(z) = -\gamma + \int_{0}^{1}\frac{1-x^{z-1}}{1-x}dx
$$

where $\gamma$ represents the Euler-Mashceroni constant, $\gamma \approx 0.57721$

## Modified Bessel function of the second kind

Modified Bessel functions of the second kind are often used in non-relativistic expansions of Fermi-dirac integrals. With some simplification relevant to the work,

$$
K_2(z) \simeq \sqrt{\frac{\pi}{2z}}\sum_{n=0}^{\infty}\frac{\Gamma(5/2 + n)}{\Gamma(5/2 - n)n! (2z)^n}
$$
