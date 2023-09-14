#include <cuComplex.h>
#include <stdio.h>
#include <stdlib.h>
#define numk 128
#define numv 1024

__global__ void cal(double *prob);

__device__ cuDoubleComplex point(double kappa, double V_0);

__device__ cuDoubleComplex cuConjugate(cuDoubleComplex x) {
  return make_cuDoubleComplex(cuCreal(x), -cuCimag((x)));
}

__device__ cuDoubleComplex cuConjugate(double x) {
  return make_cuDoubleComplex(x, 0);
}

__global__ void test(float kappa, float V_0) {
  auto pt = point(kappa, V_0);
  printf("here\n");
  printf("%f\t%f\n", cuCreal(pt), cuCimag(pt));
}

int main() {
  // const int numk = 100;
  // const int numv = 1000;
  double *d_prob;
  double *prob = (double *)malloc(numk * numv * sizeof(double));
  cudaMalloc((void **)&d_prob, numk * numv * sizeof(double));
  cal<<<128, 256>>>(d_prob);
  cudaDeviceSynchronize();
  auto err = cudaGetLastError();
  if (err != cudaSuccess) {
    printf("error: %s\n", cudaGetErrorString(err));
    exit(err);
  }
  cudaMemcpy(prob, d_prob, numk * numv * sizeof(double),
             cudaMemcpyDeviceToHost);
  // test<<<1, 1>>>(0.2341, 170);
  // cudaDeviceSynchronize();
  char res[7500001];
  int len = 0;
  for (int j = 0; j < numv; j++) {
    for (int i = 0; i < numk; i++) {
      len += sprintf(res + len, "%f\t%f\n", 0.001 + 1. / 128. * i,
                     prob[j * numk + i]);
    }
    len += sprintf(res + len, "\n\n");
  }
  FILE *fptr = fopen("prob.dat", "w");
  fprintf(fptr, "%s", res);
  fclose(fptr);
  return 0;
}

__global__ void cal(double *prob) {
  // printf("here\n");
  int i = blockIdx.x;
  for (int offset = 0; offset < 4; offset++) {
    int j = threadIdx.x * 4 + offset;
    if (i >= numk || j >= numv) {
      return;
    }
    double kappa = 0.001 + (1. / 128.) * i;
    double V_0 = 80 + 0.1 * j;
    auto pt = point(kappa, V_0);
    // printf("%f %f\n", pt.x, pt.y);
    // printf("%f\n", x);
    // printf("%f %f\n", kappa, V_0);
    prob[j * gridDim.x + i] = pt.x;
  }
  return;
}

__device__ cuDoubleComplex operator*(const cuDoubleComplex &z,
                                     const cuDoubleComplex &r) {
  return cuCmul(z, r);
}

__device__ cuDoubleComplex operator*(const cuDoubleComplex &z,
                                     const double &r) {
  return make_cuDoubleComplex(cuCreal(z) * r, cuCimag(z) * r);
}

__device__ cuDoubleComplex operator*(const double &r,
                                     const cuDoubleComplex &z) {
  return z * r;
}

__device__ cuDoubleComplex operator+(const cuDoubleComplex &z,
                                     const cuDoubleComplex &r) {
  return cuCadd(z, r);
}

__device__ cuDoubleComplex operator+(const cuDoubleComplex &z,
                                     const double &r) {
  return cuCadd(z, make_cuDoubleComplex(r, 0));
}

__device__ cuDoubleComplex operator+(const double &r,
                                     const cuDoubleComplex &z) {
  return z + r;
}

__device__ cuDoubleComplex operator-(const cuDoubleComplex &lhs, const cuDoubleComplex &rhs) {
    return make_cuDoubleComplex(lhs.x - rhs.x, lhs.y - rhs.y);
}

__device__ cuDoubleComplex operator/(const cuDoubleComplex &z,
                                     const double &a) {
  return z * (1 / a);
}
__device__ cuDoubleComplex operator/(const double &a,
                                     const cuDoubleComplex &z) {
  return a * cuConjugate(z) /
         (cuCreal(z) * cuCreal(z) - cuCimag(z) * cuCimag(z));
}

__device__ cuDoubleComplex operator/(const cuDoubleComplex &z,
                                     const cuDoubleComplex &a) {
  return z * cuConjugate(a) /
         (cuCreal(a) * cuCreal(a) - cuCimag(a) * cuCimag(a));
}

__device__ cuDoubleComplex exp(const cuDoubleComplex z) {
  double a = cuCreal(z);
  double b = cuCimag(z);
  return exp(a) * make_cuDoubleComplex(cos(b), sin(b));
}

__device__ cuDoubleComplex point(double kappa, double V_0) {
  // printf("here\n");
return -4.0960000000000009e-13*(kappa)/exp(make_cuDoubleComplex(0.0,2.0000000000000000e-02)*pow(2.0,(1.0/2.0))*pow(kappa*V_0,(1.0/2.0)))*exp(-2.0000000000000000e-02*pow(2.0,(1.0/2.0))*pow( V_0-kappa*V_0,(1.0/2.0)))*( (kappa)*(V_0)-(V_0))*kappa*exp(-2.0000000000000000e-02*(pow( V_0-kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))/( 2.4000000000000003e-07*(kappa)*( (kappa)*(V_0)-(V_0))*(V_0)+1.6000000000000000e-07*(kappa)*( (kappa)*(V_0)-(V_0))*exp(-2.0000000000000000e-02*(pow( V_0-kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*(V_0)+8.0000000000000002e-08*(kappa)*( (kappa)*(V_0)-(V_0))*exp(make_cuDoubleComplex(0.0,-2.0000000000000000e-02)*(pow(kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*(V_0)+4.0000000000000001e-08*pow((kappa),2.0)*exp(-4.0000000000000001e-02*(pow( V_0-kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*pow((V_0),2.0)+-1.6000000000000000e-07*(kappa)*( (kappa)*(V_0)-(V_0))*exp(-2.0000000000000000e-02*(pow( V_0-kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*exp(make_cuDoubleComplex(0.0,-2.0000000000000000e-02)*(pow(kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*(V_0)+8.0000000000000002e-08*pow((kappa),2.0)*exp(-2.0000000000000000e-02*(pow( V_0-kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*exp(make_cuDoubleComplex(0.0,-2.0000000000000000e-02)*(pow(kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*pow((V_0),2.0)+-8.0000000000000002e-08*pow((kappa),2.0)*exp(-2.0000000000000000e-02*(pow( V_0-kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*pow((V_0),2.0)+-4.0000000000000001e-08*pow((kappa),2.0)*exp(make_cuDoubleComplex(0.0,-2.0000000000000000e-02)*(pow(kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*pow((V_0),2.0)+-4.0000000000000001e-08*exp(-4.0000000000000001e-02*(pow( V_0-kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*pow( (kappa)*(V_0)-(V_0),2.0)*exp(make_cuDoubleComplex(0.0,-2.0000000000000000e-02)*(pow(kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))+make_cuDoubleComplex(0.0,-1.6000000000000000e-07)*(pow(kappa*V_0,(3.0/2.0)))*(pow( V_0-kappa*V_0,(1.0/2.0)))+make_cuDoubleComplex(0.0,-1.6000000000000000e-07)*(pow(kappa*V_0,(1.0/2.0)))*exp(-4.0000000000000001e-02*(pow( V_0-kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*(pow( V_0-kappa*V_0,(3.0/2.0)))+-4.0000000000000001e-08*pow( (kappa)*(V_0)-(V_0),2.0)*exp(make_cuDoubleComplex(0.0,-2.0000000000000000e-02)*(pow(kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))+4.0000000000000001e-08*exp(-4.0000000000000001e-02*(pow( V_0-kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*pow( (kappa)*(V_0)-(V_0),2.0)+8.0000000000000002e-08*(kappa)*exp(-4.0000000000000001e-02*(pow( V_0-kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*( (kappa)*(V_0)-(V_0))*exp(make_cuDoubleComplex(0.0,-2.0000000000000000e-02)*(pow(kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*(V_0)+8.0000000000000002e-08*pow( (kappa)*(V_0)-(V_0),2.0)*exp(-2.0000000000000000e-02*(pow( V_0-kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*exp(make_cuDoubleComplex(0.0,-2.0000000000000000e-02)*(pow(kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))+4.0000000000000001e-08*pow( (kappa)*(V_0)-(V_0),2.0)+-4.0000000000000001e-08*pow((kappa),2.0)*exp(-4.0000000000000001e-02*(pow( V_0-kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*exp(make_cuDoubleComplex(0.0,-2.0000000000000000e-02)*(pow(kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*pow((V_0),2.0)+2.4000000000000003e-07*(kappa)*exp(-4.0000000000000001e-02*(pow( V_0-kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*( (kappa)*(V_0)-(V_0))*(V_0)+4.0000000000000001e-08*pow((kappa),2.0)*pow((V_0),2.0)+make_cuDoubleComplex(0.0,1.6000000000000000e-07)*(pow(kappa*V_0,(1.0/2.0)))*(pow( V_0-kappa*V_0,(3.0/2.0)))+-8.0000000000000002e-08*pow( (kappa)*(V_0)-(V_0),2.0)*exp(-2.0000000000000000e-02*(pow( V_0-kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))+make_cuDoubleComplex(0.0,1.6000000000000000e-07)*(pow(kappa*V_0,(3.0/2.0)))*exp(-4.0000000000000001e-02*(pow( V_0-kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*(pow( V_0-kappa*V_0,(1.0/2.0))))*V_0/exp(make_cuDoubleComplex(0.0,-2.0000000000000000e-02)*(pow(kappa*V_0,(1.0/2.0)))*pow(2.0,(1.0/2.0)))*(V_0)/( make_cuDoubleComplex(0.0,-1.6000000000000000e-07)*pow( V_0-kappa*V_0,(3.0/2.0))*pow(kappa*V_0,(1.0/2.0))+-2.4000000000000003e-07*exp(-4.0000000000000001e-02*pow(2.0,(1.0/2.0))*pow( V_0-kappa*V_0,(1.0/2.0)))*kappa*V_0*( V_0-kappa*V_0)+-8.0000000000000002e-08*exp(make_cuDoubleComplex(0.0,2.0000000000000000e-02)*pow(2.0,(1.0/2.0))*pow(kappa*V_0,(1.0/2.0)))*kappa*V_0*( V_0-kappa*V_0)+4.0000000000000001e-08*exp(-4.0000000000000001e-02*pow(2.0,(1.0/2.0))*pow( V_0-kappa*V_0,(1.0/2.0)))*pow( V_0-kappa*V_0,2.0)+make_cuDoubleComplex(0.0,1.6000000000000000e-07)*exp(-4.0000000000000001e-02*pow(2.0,(1.0/2.0))*pow( V_0-kappa*V_0,(1.0/2.0)))*pow( V_0-kappa*V_0,(3.0/2.0))*pow(kappa*V_0,(1.0/2.0))+make_cuDoubleComplex(0.0,1.6000000000000000e-07)*pow( V_0-kappa*V_0,(1.0/2.0))*pow(kappa*V_0,(3.0/2.0))+4.0000000000000001e-08*exp(-4.0000000000000001e-02*pow(2.0,(1.0/2.0))*pow( V_0-kappa*V_0,(1.0/2.0)))*(kappa*kappa)*(V_0*V_0)+-4.0000000000000001e-08*exp(make_cuDoubleComplex(0.0,2.0000000000000000e-02)*pow(2.0,(1.0/2.0))*pow(kappa*V_0,(1.0/2.0)))*pow( V_0-kappa*V_0,2.0)+-2.4000000000000003e-07*kappa*V_0*( V_0-kappa*V_0)+1.6000000000000000e-07*exp(make_cuDoubleComplex(0.0,2.0000000000000000e-02)*pow(2.0,(1.0/2.0))*pow(kappa*V_0,(1.0/2.0)))*exp(-2.0000000000000000e-02*pow(2.0,(1.0/2.0))*pow( V_0-kappa*V_0,(1.0/2.0)))*kappa*V_0*( V_0-kappa*V_0)+-8.0000000000000002e-08*exp(-2.0000000000000000e-02*pow(2.0,(1.0/2.0))*pow( V_0-kappa*V_0,(1.0/2.0)))*pow( V_0-kappa*V_0,2.0)+make_cuDoubleComplex(0.0,-1.6000000000000000e-07)*exp(-4.0000000000000001e-02*pow(2.0,(1.0/2.0))*pow( V_0-kappa*V_0,(1.0/2.0)))*pow( V_0-kappa*V_0,(1.0/2.0))*pow(kappa*V_0,(3.0/2.0))+-4.0000000000000001e-08*exp(-4.0000000000000001e-02*pow(2.0,(1.0/2.0))*pow( V_0-kappa*V_0,(1.0/2.0)))*exp(make_cuDoubleComplex(0.0,2.0000000000000000e-02)*pow(2.0,(1.0/2.0))*pow(kappa*V_0,(1.0/2.0)))*pow( V_0-kappa*V_0,2.0)+4.0000000000000001e-08*(kappa*kappa)*(V_0*V_0)+-4.0000000000000001e-08*exp(-4.0000000000000001e-02*pow(2.0,(1.0/2.0))*pow( V_0-kappa*V_0,(1.0/2.0)))*exp(make_cuDoubleComplex(0.0,2.0000000000000000e-02)*pow(2.0,(1.0/2.0))*pow(kappa*V_0,(1.0/2.0)))*(kappa*kappa)*(V_0*V_0)+-8.0000000000000002e-08*exp(-4.0000000000000001e-02*pow(2.0,(1.0/2.0))*pow( V_0-kappa*V_0,(1.0/2.0)))*exp(make_cuDoubleComplex(0.0,2.0000000000000000e-02)*pow(2.0,(1.0/2.0))*pow(kappa*V_0,(1.0/2.0)))*kappa*V_0*( V_0-kappa*V_0)+-8.0000000000000002e-08*exp(-2.0000000000000000e-02*pow(2.0,(1.0/2.0))*pow( V_0-kappa*V_0,(1.0/2.0)))*(kappa*kappa)*(V_0*V_0)+4.0000000000000001e-08*pow( V_0-kappa*V_0,2.0)+8.0000000000000002e-08*exp(make_cuDoubleComplex(0.0,2.0000000000000000e-02)*pow(2.0,(1.0/2.0))*pow(kappa*V_0,(1.0/2.0)))*exp(-2.0000000000000000e-02*pow(2.0,(1.0/2.0))*pow( V_0-kappa*V_0,(1.0/2.0)))*pow( V_0-kappa*V_0,2.0)+-4.0000000000000001e-08*exp(make_cuDoubleComplex(0.0,2.0000000000000000e-02)*pow(2.0,(1.0/2.0))*pow(kappa*V_0,(1.0/2.0)))*(kappa*kappa)*(V_0*V_0)+8.0000000000000002e-08*exp(make_cuDoubleComplex(0.0,2.0000000000000000e-02)*pow(2.0,(1.0/2.0))*pow(kappa*V_0,(1.0/2.0)))*exp(-2.0000000000000000e-02*pow(2.0,(1.0/2.0))*pow( V_0-kappa*V_0,(1.0/2.0)))*(kappa*kappa)*(V_0*V_0)+-1.6000000000000000e-07*exp(-2.0000000000000000e-02*pow(2.0,(1.0/2.0))*pow( V_0-kappa*V_0,(1.0/2.0)))*kappa*V_0*( V_0-kappa*V_0))*( V_0-kappa*V_0)
;
}
