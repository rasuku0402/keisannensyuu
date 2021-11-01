#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double Combi(double n, int r){

  double i;
  double c = 1.0;  

  if(r==0){
    return c;
  }else{
    for(i=n;i>=n-r+1;i--){
      c=c*i;
  }
  for(i=1;i<=r;i++){
      c=c/i;
  }
  return c;
  }
}

double legendre(double n, int k){
    int i;
    int a = 1;
    for(i=n;i>0;i--){
        a = 2*a;
    }
    double c = Combi(n,k);
    double d = Combi(((n+k-1)/2), n);
    return a*c*d;
}

/* ルジャンドル多項式　P_n(x) */
double P(double x, double n){
    int k;
    double c[12];
    double d;
    for(k=0;k<=n;k++){
        c[k] = legendre(n,k);
        d = c[k]*pow(x,k);
    }
    return d;
}

/* 二分法 初期値 x1<x2 と 誤差限界 eps を入力 */
double bisec(double x1, double x2, double eps) {
    double c;
    while (x2 - x1 >= eps) {
        c = (x1+x2)/2.0; /* 中点計算 */
        if (P(x1,1)*P(c,1) > 0.0) {
             /* 同符号か判定 */
            x1 = c;
            printf("esp = %lf\n", x1);
        } else {

            x2 = c;
            printf("%.10lf", x2);
            break;
        }
    }
    return c;
}

int main(void){
    double eps=1e-9;
    printf("solution:%lf\n", bisec(-0.2,0.2,eps));
    return 0;
}