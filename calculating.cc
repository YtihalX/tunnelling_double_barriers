#include <complex>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

using namespace ::std;

#define NUM_THREADS 12
#define NUM_V0 256
#define NUM_KAPPA 196608

void *routine(void *arg);
complex<double> f(double kappa, double V_0);

char res[NUM_V0 * (NUM_KAPPA + 1) * 27];
double value[NUM_KAPPA * NUM_V0];
struct parg {
  double *arr;
  int thread_id;
};

int main() {

  pthread_t threads[NUM_THREADS];
  parg args[NUM_THREADS];

  for (int i = 0; i < NUM_THREADS; i++) {

    args[i].arr = value;
    args[i].thread_id = i;

    if (pthread_create(threads + i, NULL, routine, (void *)(args + i)) != 0) {
      perror("pthread_create");
      exit(EXIT_FAILURE);
    }
  }

  for (int i = 0; i < NUM_THREADS; i++) {
    if (pthread_join(threads[i], nullptr) != 0) {
      perror("pthread_join");
      exit(EXIT_FAILURE);
    }
  }

  int len = 0;
  for (int i = 0; i < NUM_V0; i++) {
    for (int j = 0; j < NUM_KAPPA; j++) {
      len += sprintf(res + len, "%f\t%f\n",
                     1.0172526041666666e-5 + 2. / NUM_KAPPA * j,
                     value[i * NUM_KAPPA + j]);
    }
    len += sprintf(res + len, "\n\n");
  }
  FILE *fptr = fopen("prob.dat", "w");
  fprintf(fptr, "%s", res);
  fclose(fptr);

  fptr = popen("gnuplot -persistent", "w");
  fprintf(fptr, "set xlabel \"kappa\"\n");
  fprintf(fptr, "set ylabel \"probability\"\n");
  fprintf(fptr, "set term gif animate size 1280, 960 delay 0\n");
  fprintf(fptr, "set output \"animated_plot.gif\"\n");
  fprintf(fptr, "do for [i=0:%d] {\n", NUM_V0 - 1);
  fprintf(fptr, "    plot \"prob.dat\" index i using 1:2 with lines title "
                "\"probability of passing through\"\n");
  fprintf(fptr, "}\n");
  fflush(fptr);
  pclose(fptr);
  return 0;
}

void *routine(void *arg) {
  parg warg = *(parg *)arg;
  int t = (NUM_KAPPA + NUM_THREADS - 1) / NUM_THREADS;
  for (int i = 0; i < t; i++) {
    int column = i * NUM_THREADS + warg.thread_id;
    if (column >= NUM_KAPPA)
      continue;
    double kappa = 1.0172526041666666e-5 + 2. / NUM_KAPPA * column;
    for (int row = 0; row < NUM_V0; row++) {
      // printf("i: %d, j: %d\n", row, column);
      double V_0 = 3000 + 1320000. / NUM_V0 * row;
      warg.arr[row * NUM_KAPPA + column] = f(kappa, V_0).real();
    }
  }
  return NULL;
}

complex<double> f(double kappa, double V_0) {
  return 4.0960000000000009e-13 * V_0 /
         exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
             pow(2.0, (1.0 / 2.0)) * conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
         conj(V_0) *
         exp(2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
             conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                      (1.0 / 2.0)))) /
         (-1.6000000000000000e-07 *
              exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
              conj(V_0) *
              exp(2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                  conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                           (1.0 / 2.0)))) *
              conj(kappa) * (conj(V_0) - conj(V_0) * conj(kappa)) +
          -8.0000000000000002e-08 *
              exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
              pow(conj(V_0), 2.0) *
              exp(2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                  conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                           (1.0 / 2.0)))) *
              pow(conj(kappa), 2.0) +
          -4.0000000000000001e-08 * pow(conj(V_0), 2.0) *
              exp(4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                  conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                           (1.0 / 2.0)))) *
              pow(conj(kappa), 2.0) +
          -4.0000000000000001e-08 *
              pow(conj(V_0) - conj(V_0) * conj(kappa), 2.0) +
          2.4000000000000003e-07 * conj(V_0) *
              exp(4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                  conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                           (1.0 / 2.0)))) *
              conj(kappa) * (conj(V_0) - conj(V_0) * conj(kappa)) +
          8.0000000000000002e-08 *
              exp(2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                  conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                           (1.0 / 2.0)))) *
              pow(conj(V_0) - conj(V_0) * conj(kappa), 2.0) +
          -4.0000000000000001e-08 *
              exp(4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                  conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                           (1.0 / 2.0)))) *
              pow(conj(V_0) - conj(V_0) * conj(kappa), 2.0) +
          std::complex<double>(0.0, 1.6000000000000000e-07) *
              conj(pow(V_0 * kappa, (3.0 / 2.0))) *
              exp(4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                  conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                           (1.0 / 2.0)))) *
              conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                       (1.0 / 2.0))) +
          4.0000000000000001e-08 *
              exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
              pow(conj(V_0), 2.0) * pow(conj(kappa), 2.0) +
          8.0000000000000002e-08 * pow(conj(V_0), 2.0) *
              exp(2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                  conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                           (1.0 / 2.0)))) *
              pow(conj(kappa), 2.0) +
          4.0000000000000001e-08 *
              exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
              pow(conj(V_0), 2.0) *
              exp(4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                  conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                           (1.0 / 2.0)))) *
              pow(conj(kappa), 2.0) +
          std::complex<double>(0.0, -1.6000000000000000e-07) *
              conj(pow(V_0 * kappa, (3.0 / 2.0))) *
              conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                       (1.0 / 2.0))) +
          2.4000000000000003e-07 * conj(V_0) * conj(kappa) *
              (conj(V_0) - conj(V_0) * conj(kappa)) +
          -8.0000000000000002e-08 *
              exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
              exp(2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                  conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                           (1.0 / 2.0)))) *
              pow(conj(V_0) - conj(V_0) * conj(kappa), 2.0) +
          4.0000000000000001e-08 *
              exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
              pow(conj(V_0) - conj(V_0) * conj(kappa), 2.0) +
          std::complex<double>(0.0, -1.6000000000000000e-07) *
              conj(pow(V_0 * kappa, (1.0 / 2.0))) *
              exp(4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                  conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                           (1.0 / 2.0)))) *
              conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                       (3.0 / 2.0))) +
          1.6000000000000000e-07 * conj(V_0) *
              exp(2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                  conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                           (1.0 / 2.0)))) *
              conj(kappa) * (conj(V_0) - conj(V_0) * conj(kappa)) +
          8.0000000000000002e-08 *
              exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
              conj(V_0) * conj(kappa) * (conj(V_0) - conj(V_0) * conj(kappa)) +
          -4.0000000000000001e-08 * pow(conj(V_0), 2.0) *
              pow(conj(kappa), 2.0) +
          4.0000000000000001e-08 *
              exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
              exp(4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                  conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                           (1.0 / 2.0)))) *
              pow(conj(V_0) - conj(V_0) * conj(kappa), 2.0) +
          8.0000000000000002e-08 *
              exp(std::complex<double>(0.0, -2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * conj(pow(V_0 * kappa, (1.0 / 2.0)))) *
              conj(V_0) *
              exp(4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                  conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                           (1.0 / 2.0)))) *
              conj(kappa) * (conj(V_0) - conj(V_0) * conj(kappa)) +
          std::complex<double>(0.0, 1.6000000000000000e-07) *
              conj(pow(V_0 * kappa, (1.0 / 2.0))) *
              conj(pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                       (3.0 / 2.0)))) /
         exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
             pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
         (std::complex<double>(V_0 - V_0 * kappa, 0.0)) * conj(kappa) *
         (conj(V_0) - conj(V_0) * conj(kappa)) /
         (8.0000000000000002e-08 * V_0 *
              exp(4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                  pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                      (1.0 / 2.0))) *
              exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
              (std::complex<double>(V_0 - V_0 * kappa, 0.0)) * kappa +
          4.0000000000000001e-08 *
              exp(4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                  pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                      (1.0 / 2.0))) *
              exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
              pow(std::complex<double>(V_0 - V_0 * kappa, 0.0), 2.0) +
          8.0000000000000002e-08 * (V_0 * V_0) * (kappa * kappa) *
              exp(2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                  pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                      (1.0 / 2.0))) +
          4.0000000000000001e-08 *
              exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
              pow(std::complex<double>(V_0 - V_0 * kappa, 0.0), 2.0) +
          2.4000000000000003e-07 * V_0 *
              (std::complex<double>(V_0 - V_0 * kappa, 0.0)) * kappa +
          -4.0000000000000001e-08 *
              pow(std::complex<double>(V_0 - V_0 * kappa, 0.0), 2.0) +
          8.0000000000000002e-08 *
              pow(std::complex<double>(V_0 - V_0 * kappa, 0.0), 2.0) *
              exp(2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                  pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                      (1.0 / 2.0))) +
          -4.0000000000000001e-08 * (V_0 * V_0) *
              exp(4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                  pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                      (1.0 / 2.0))) *
              (kappa * kappa) +
          -4.0000000000000001e-08 * (V_0 * V_0) * (kappa * kappa) +
          2.4000000000000003e-07 * V_0 *
              exp(4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                  pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                      (1.0 / 2.0))) *
              (std::complex<double>(V_0 - V_0 * kappa, 0.0)) * kappa +
          std::complex<double>(0.0, -1.6000000000000000e-07) *
              pow(std::complex<double>(V_0 - V_0 * kappa, 0.0), (3.0 / 2.0)) *
              pow(V_0 * kappa, (1.0 / 2.0)) +
          8.0000000000000002e-08 * V_0 *
              exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
              (std::complex<double>(V_0 - V_0 * kappa, 0.0)) * kappa +
          4.0000000000000001e-08 * (V_0 * V_0) *
              exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
              (kappa * kappa) +
          1.6000000000000000e-07 * V_0 *
              (std::complex<double>(V_0 - V_0 * kappa, 0.0)) * kappa *
              exp(2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                  pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                      (1.0 / 2.0))) +
          4.0000000000000001e-08 * (V_0 * V_0) *
              exp(4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                  pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                      (1.0 / 2.0))) *
              exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
              (kappa * kappa) +
          -8.0000000000000002e-08 *
              exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
              pow(std::complex<double>(V_0 - V_0 * kappa, 0.0), 2.0) *
              exp(2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                  pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                      (1.0 / 2.0))) +
          std::complex<double>(0.0, 1.6000000000000000e-07) *
              exp(4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                  pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                      (1.0 / 2.0))) *
              pow(std::complex<double>(V_0 - V_0 * kappa, 0.0), (3.0 / 2.0)) *
              pow(V_0 * kappa, (1.0 / 2.0)) +
          std::complex<double>(0.0, 1.6000000000000000e-07) *
              pow(std::complex<double>(V_0 - V_0 * kappa, 0.0), (1.0 / 2.0)) *
              pow(V_0 * kappa, (3.0 / 2.0)) +
          -4.0000000000000001e-08 *
              exp(4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                  pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                      (1.0 / 2.0))) *
              pow(std::complex<double>(V_0 - V_0 * kappa, 0.0), 2.0) +
          std::complex<double>(0.0, -1.6000000000000000e-07) *
              exp(4.0000000000000001e-02 * pow(2.0, (1.0 / 2.0)) *
                  pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                      (1.0 / 2.0))) *
              pow(std::complex<double>(V_0 - V_0 * kappa, 0.0), (1.0 / 2.0)) *
              pow(V_0 * kappa, (3.0 / 2.0)) +
          -8.0000000000000002e-08 * (V_0 * V_0) *
              exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
              (kappa * kappa) *
              exp(2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                  pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                      (1.0 / 2.0))) +
          -1.6000000000000000e-07 * V_0 *
              exp(std::complex<double>(0.0, 2.0000000000000000e-02) *
                  pow(2.0, (1.0 / 2.0)) * pow(V_0 * kappa, (1.0 / 2.0))) *
              (std::complex<double>(V_0 - V_0 * kappa, 0.0)) * kappa *
              exp(2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
                  pow(std::complex<double>(V_0 - V_0 * kappa, 0.0),
                      (1.0 / 2.0)))) *
         kappa *
         exp(2.0000000000000000e-02 * pow(2.0, (1.0 / 2.0)) *
             pow(std::complex<double>(V_0 - V_0 * kappa, 0.0), (1.0 / 2.0)));
}
