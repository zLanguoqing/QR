#ifndef _QR_H_
#define _QR_H_
#define MIN (1e-30)
int QR(double a[], int m, int n, double q[]);
void Hessenberg(double a[], int n);
int HessenbergQR(double a[],int n,double u[],double v[],double eps,int Iteration);
#endif