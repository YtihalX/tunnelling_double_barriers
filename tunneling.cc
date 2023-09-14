#include <ginac/ex.h>
#include <ginac/ginac.h>
#include <ginac/numeric.h>
#include <ginac/operators.h>
#include <ginac/print.h>
#include <iostream>
#include <string>
#include <vector>

using namespace ::GiNaC;
using namespace ::std;

#define NUM_V0 512
#define NUM_KAPPA 512

int main() {
  symbol A[4], B[4];
  for (int i = 0; i < 4; i++) {
    A[i] = symbol("A" + to_string(i));
    B[i] = symbol("B" + to_string(i));
  }
  symbol k0("k0"), k1("k1");
  symbol V_0("V_0"), a("a"), b("b"), m("m"), E("E");
  symbol x("x"), kappa("kappa");

  ex psi1 = exp(k0 * x * I) + A[0] * exp(-k0 * x * I);
  ex psi2 = B[0] * exp(k1 * x) + B[1] * exp(-k1 * x);
  ex psi3 = A[1] * exp(k0 * x * I) + A[2] * exp(-k0 * x * I);
  ex psi4 = B[2] * exp(k1 * x) + B[3] * exp(-k1 * x);
  ex psi5 = A[3] * exp(k0 * x * I);

  lst eqs, vars;
  eqs.append(psi1.subs(x == 0) == psi2.subs(x == 0));
  eqs.append(psi1.diff(x).subs(x == 0) == psi2.diff(x).subs(x == 0));
  eqs.append(psi2.subs(x == a) == psi3.subs(x == a));
  eqs.append(psi2.diff(x).subs(x == a) == psi3.diff(x).subs(x == a));
  eqs.append(psi3.subs(x == a + a) == psi4.subs(x == a + a));
  eqs.append(psi3.diff(x).subs(x == a + a) == psi4.diff(x).subs(x == a + a));
  eqs.append(psi4.subs(x == 2 * a + a) == psi5.subs(x == 2 * a + a));
  eqs.append(psi4.diff(x).subs(x == 2 * a + a) ==
             psi5.diff(x).subs(x == 2 * a + a));

  for (int i = 0; i < 4; i++) {
    vars.append(A[i]);
    vars.append(B[i]);
  }
  ex sol = lsolve(eqs, vars);
  ex R = expand(sol[6].rhs());
  R = R.subs(
      lst{k0 == sqrt(2 * m * E) * 0.01, k1 == sqrt(2 * m * (V_0 - E)) * 0.01});
  R = R.subs(lst{E == kappa * V_0, m == 1, a == 1, b == 1});
  ex prob = R * conjugate(R);
  cout << csrc_double << prob << '\n';

  return 0;
  // cout<<conjugate(R.subs(lst{a == 1, b == 1, m == 1}))*R<<'\n';
}
