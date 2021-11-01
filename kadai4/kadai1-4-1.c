#include <stdio.h>
#include <math.h>
#include <string.h>

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
        d = d + c[k]*pow(x,k);
    }
    return d;
}

int main(){
    double x = 2;
    int n = 3;
    double c = P(x,n);
    printf("%.10lf\n", c);
    return 0;



}
