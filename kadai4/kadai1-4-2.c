#include <stdio.h>
#include <math.h>
#include <string.h>

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

/* ルジャンドル多項式　P_n(x)を返す関数 */
double P(double x, double n){
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

double dPdx(double x, double n){
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
    int i,j;
    double x_0 = 0; /*ニュートン法の初期値をx_0とする*/
    double x;
    for(i=5;i<=10;i++){
        printf("n=%d\n", i);
            for(j=0;j<i;j++){
                double esp = 100;
                x_0 = cos((M_PI*(j+0.75))/(i+0.5)); /*ニュートン法の初期値を与式により更新*/

                while(esp > 1e-9){
                    x = x_0;
                    x_0 = x_0 - P(x_0,i)/dPdx(x_0,i);
                    esp = fabs(x_0 - x);
                }
                printf("%.9lf\n", x_0);
            }
    }
    
    return 0;
}