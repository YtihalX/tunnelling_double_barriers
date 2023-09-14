#include <complex>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
using namespace ::std;
int main(int argc, char *argv[]) {

  double V_0, kappa;
  sscanf(argv[1], "%lf", &kappa);
  sscanf(argv[2], "%lf", &V_0);
  cout << -4.0960000000000009e-13 * V_0 /
              exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * conj(pow(V_0 * kappa, (1.0 / 2.0)))) /
              exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
              (V_0 - V_0 * kappa) *
              exp(-2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                  conj(pow(V_0 - V_0 * kappa, (1.0 / 2.0)))) *
              conj(V_0) * conj(kappa) /
              (-8.0000000000000002e-08 *
                   exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
                   exp(-2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 - V_0 * kappa, (1.0 / 2.0)))) *
                   pow(conj(V_0), 2.0) * pow(conj(kappa), 2.0) +
               4.0000000000000001e-08 *
                   exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
                   exp(-4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 - V_0 * kappa, (1.0 / 2.0)))) *
                   pow(conj(V_0) * conj(kappa) - conj(V_0), 2.0) +
               -8.0000000000000002e-08 *
                   exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
                   conj(V_0) *
                   exp(-4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 - V_0 * kappa, (1.0 / 2.0)))) *
                   conj(kappa) * (conj(V_0) * conj(kappa) - conj(V_0)) +
               8.0000000000000002e-08 *
                   exp(-2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 - V_0 * kappa, (1.0 / 2.0)))) *
                   pow(conj(V_0) * conj(kappa) - conj(V_0), 2.0) +
               std::complex<double>(0.0, 1.6000000000000000e-07) *
                   exp(-4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 - V_0 * kappa, (1.0 / 2.0)))) *
                   conj(pow(V_0 - V_0 * kappa, (3.0 / 2.0))) *
                   conj(pow(V_0 * kappa, (1.0 / 2.0))) +
               std::complex<double>(0.0, -1.6000000000000000e-07) *
                   exp(-4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 - V_0 * kappa, (1.0 / 2.0)))) *
                   conj(pow(V_0 * kappa, (3.0 / 2.0))) *
                   conj(pow(V_0 - V_0 * kappa, (1.0 / 2.0))) +
               -2.4000000000000003e-07 * conj(V_0) *
                   exp(-4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 - V_0 * kappa, (1.0 / 2.0)))) *
                   conj(kappa) * (conj(V_0) * conj(kappa) - conj(V_0)) +
               4.0000000000000001e-08 *
                   exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
                   pow(conj(V_0), 2.0) *
                   exp(-4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 - V_0 * kappa, (1.0 / 2.0)))) *
                   pow(conj(kappa), 2.0) +
               -4.0000000000000001e-08 *
                   exp(-4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 - V_0 * kappa, (1.0 / 2.0)))) *
                   pow(conj(V_0) * conj(kappa) - conj(V_0), 2.0) +
               -8.0000000000000002e-08 *
                   exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
                   exp(-2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 - V_0 * kappa, (1.0 / 2.0)))) *
                   pow(conj(V_0) * conj(kappa) - conj(V_0), 2.0) +
               -8.0000000000000002e-08 *
                   exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
                   conj(V_0) * conj(kappa) *
                   (conj(V_0) * conj(kappa) - conj(V_0)) +
               -1.6000000000000000e-07 *
                   exp(-2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 - V_0 * kappa, (1.0 / 2.0)))) *
                   conj(V_0) * conj(kappa) *
                   (conj(V_0) * conj(kappa) - conj(V_0)) +
               4.0000000000000001e-08 *
                   exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
                   pow(conj(V_0) * conj(kappa) - conj(V_0), 2.0) +
               std::complex<double>(0.0, 1.6000000000000000e-07) *
                   conj(pow(V_0 * kappa, (3.0 / 2.0))) *
                   conj(pow(V_0 - V_0 * kappa, (1.0 / 2.0))) +
               -4.0000000000000001e-08 * pow(conj(V_0), 2.0) *
                   pow(conj(kappa), 2.0) +
               -4.0000000000000001e-08 *
                   pow(conj(V_0) * conj(kappa) - conj(V_0), 2.0) +
               8.0000000000000002e-08 *
                   exp(-2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 - V_0 * kappa, (1.0 / 2.0)))) *
                   pow(conj(V_0), 2.0) * pow(conj(kappa), 2.0) +
               4.0000000000000001e-08 *
                   exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
                   pow(conj(V_0), 2.0) * pow(conj(kappa), 2.0) +
               -2.4000000000000003e-07 * conj(V_0) * conj(kappa) *
                   (conj(V_0) * conj(kappa) - conj(V_0)) +
               1.6000000000000000e-07 *
                   exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
                   exp(-2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 - V_0 * kappa, (1.0 / 2.0)))) *
                   conj(V_0) * conj(kappa) *
                   (conj(V_0) * conj(kappa) - conj(V_0)) +
               std::complex<double>(0.0, -1.6000000000000000e-07) *
                   conj(pow(V_0 - V_0 * kappa, (3.0 / 2.0))) *
                   conj(pow(V_0 * kappa, (1.0 / 2.0))) +
               -4.0000000000000001e-08 * pow(conj(V_0), 2.0) *
                   exp(-4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                       conj(pow(V_0 - V_0 * kappa, (1.0 / 2.0)))) *
                   pow(conj(kappa), 2.0)) /
              (-4.0000000000000001e-08 * pow(V_0 - V_0 * kappa, 2.0) +
               -8.0000000000000002e-08 * (V_0 * V_0) *
                   exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
                   (kappa * kappa) *
                   exp(-2.0000000000000000e-02 *
                       pow(V_0 - V_0 * kappa, (1.0 / 2.0)) *
                       pow(2.0, (1.0 / 2.0))) +
               4.0000000000000001e-08 *
                   exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
                   pow(V_0 - V_0 * kappa, 2.0) +
               -4.0000000000000001e-08 * (V_0 * V_0) * (kappa * kappa) +
               -4.0000000000000001e-08 * (V_0 * V_0) *
                   exp(-4.0000000000000001e-02 *
                       pow(V_0 - V_0 * kappa, (1.0 / 2.0)) *
                       pow(2.0, (1.0 / 2.0))) *
                   (kappa * kappa) +
               8.0000000000000002e-08 * V_0 *
                   exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
                   (V_0 - V_0 * kappa) *
                   exp(-4.0000000000000001e-02 *
                       pow(V_0 - V_0 * kappa, (1.0 / 2.0)) *
                       pow(2.0, (1.0 / 2.0))) *
                   kappa +
               8.0000000000000002e-08 * V_0 *
                   exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
                   (V_0 - V_0 * kappa) * kappa +
               4.0000000000000001e-08 *
                   exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
                   pow(V_0 - V_0 * kappa, 2.0) *
                   exp(-4.0000000000000001e-02 *
                       pow(V_0 - V_0 * kappa, (1.0 / 2.0)) *
                       pow(2.0, (1.0 / 2.0))) +
               -1.6000000000000000e-07 * V_0 *
                   exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
                   (V_0 - V_0 * kappa) * kappa *
                   exp(-2.0000000000000000e-02 *
                       pow(V_0 - V_0 * kappa, (1.0 / 2.0)) *
                       pow(2.0, (1.0 / 2.0))) +
               -4.0000000000000001e-08 * pow(V_0 - V_0 * kappa, 2.0) *
                   exp(-4.0000000000000001e-02 *
                       pow(V_0 - V_0 * kappa, (1.0 / 2.0)) *
                       pow(2.0, (1.0 / 2.0))) +
               std::complex<double>(0.0, -1.6000000000000000e-07) *
                   pow(V_0 - V_0 * kappa, (3.0 / 2.0)) *
                   pow(V_0 * kappa, (1.0 / 2.0)) *
                   exp(-4.0000000000000001e-02 *
                       pow(V_0 - V_0 * kappa, (1.0 / 2.0)) *
                       pow(2.0, (1.0 / 2.0))) +
               2.4000000000000003e-07 * V_0 * (V_0 - V_0 * kappa) *
                   exp(-4.0000000000000001e-02 *
                       pow(V_0 - V_0 * kappa, (1.0 / 2.0)) *
                       pow(2.0, (1.0 / 2.0))) *
                   kappa +
               std::complex<double>(0.0, 1.6000000000000000e-07) *
                   pow(V_0 - V_0 * kappa, (3.0 / 2.0)) *
                   pow(V_0 * kappa, (1.0 / 2.0)) +
               4.0000000000000001e-08 * (V_0 * V_0) *
                   exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
                   exp(-4.0000000000000001e-02 *
                       pow(V_0 - V_0 * kappa, (1.0 / 2.0)) *
                       pow(2.0, (1.0 / 2.0))) *
                   (kappa * kappa) +
               1.6000000000000000e-07 * V_0 * (V_0 - V_0 * kappa) * kappa *
                   exp(-2.0000000000000000e-02 *
                       pow(V_0 - V_0 * kappa, (1.0 / 2.0)) *
                       pow(2.0, (1.0 / 2.0))) +
               2.4000000000000003e-07 * V_0 * (V_0 - V_0 * kappa) * kappa +
               8.0000000000000002e-08 * (V_0 * V_0) * (kappa * kappa) *
                   exp(-2.0000000000000000e-02 *
                       pow(V_0 - V_0 * kappa, (1.0 / 2.0)) *
                       pow(2.0, (1.0 / 2.0))) +
               std::complex<double>(0.0, -1.6000000000000000e-07) *
                   pow(V_0 - V_0 * kappa, (1.0 / 2.0)) *
                   pow(V_0 * kappa, (3.0 / 2.0)) +
               8.0000000000000002e-08 * pow(V_0 - V_0 * kappa, 2.0) *
                   exp(-2.0000000000000000e-02 *
                       pow(V_0 - V_0 * kappa, (1.0 / 2.0)) *
                       pow(2.0, (1.0 / 2.0))) +
               4.0000000000000001e-08 * (V_0 * V_0) *
                   exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
                   (kappa * kappa) +
               std::complex<double>(0.0, 1.6000000000000000e-07) *
                   pow(V_0 - V_0 * kappa, (1.0 / 2.0)) *
                   pow(V_0 * kappa, (3.0 / 2.0)) *
                   exp(-4.0000000000000001e-02 *
                       pow(V_0 - V_0 * kappa, (1.0 / 2.0)) *
                       pow(2.0, (1.0 / 2.0))) +
               -8.0000000000000002e-08 *
                   exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                       pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
                   pow(V_0 - V_0 * kappa, 2.0) *
                   exp(-2.0000000000000000e-02 *
                       pow(V_0 - V_0 * kappa, (1.0 / 2.0)) *
                       pow(2.0, (1.0 / 2.0)))) *
              kappa *
              exp(-2.0000000000000000e-02 *
                  pow(V_0 - V_0 * kappa, (1.0 / 2.0)) * pow(2.0, (1.0 / 2.0))) *
              (conj(V_0) * conj(kappa) - conj(V_0))

       << '\n';
}
