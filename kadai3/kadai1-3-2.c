#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*二項係数を計算するCombi関数*/

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

/*ルジャンドル多項式の係数a_n_kを計算する関数*/


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

/* ルジャンドル多項式P_n(x)を返す関数 */
double f(double x, double n){
    int k;
    double c[12]= {0};
    double d = 0;
    for(k=0;k<=n;k++){
        c[k] = legendre(n,k);
        d = d + c[k]*pow(x,k);
    }
    return d;
}

/*dP/dxを計算を返す関数*/

double df(double x, double n){
    double c[12] = {0};
    double d[12] = {0};
    int k;
    double e = 0;
    for(k=0;k<=n;k++){
        c[k] = legendre(n,k);
        d[k] = k*c[k];   /*微分した係数*/
        e = e + d[k]*pow(x,k-1);  /*xのk-1乗を加算していく*/
    }

    return e;
}

int main(){
    double x_st5 = 0.85;
    double x_st4 = 0.5;
    double x_st3 = 0.1;
    double x_st2 = -0.60;
    double x_st1 = -0.85;
    double esp = 100;

    double x;

    while(esp > 1e-9){
        x = x_st1;
        x_st1 = x_st1 - f(x_st1,5)/df(x_st1,5);
        esp = fabs(x_st1 - x);
        printf("eps = %.10lf\n", esp);
    }
    printf("solution1:%.9lf\n", x_st1);
    esp = 100;
    while(esp > 1e-9){
        x = x_st2;
        x_st2 = x_st2 - f(x_st2,5)/df(x_st2,5);
        esp = fabs(x_st2 - x);
        printf("eps = %.10lf\n", esp);
    }
    printf("solution2:%.9lf\n", x_st2);
    esp = 100;
    while(esp > 1e-9){
        x = x_st3;
        x_st3 = x_st3 - f(x_st3,5)/df(x_st3,5);
        esp = fabs(x_st3 - x);
        printf("eps = %.10lf\n", esp);
    }
    printf("solution3:%.9lf\n", x_st3);
    esp = 100;
    while(esp > 1e-9){
        x = x_st4;
        x_st4 = x_st4 - f(x_st4,5)/df(x_st4,5);
        esp = fabs(x_st4 - x);
        printf("eps = %.10lf\n", esp);
    }
    printf("solution4:%.9lf\n", x_st4);
    esp = 100;
    while(esp > 1e-9){
        x = x_st5;
        x_st5 = x_st5 - f(x_st5,5)/df(x_st5,5);
        esp = fabs(x_st5 - x);
        printf("eps = %.10lf\n", esp);
    }
    printf("solution5:%.9lf\n", x_st5);

    return 0;
}