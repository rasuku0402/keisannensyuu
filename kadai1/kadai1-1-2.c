#include <stdio.h>
#include <math.h>

double kaijou(double k){
  int i;
  double c = 1;
  for(i=k;i>=0;i--){
    c = c*i;
  }
  return c;
}

/*二項係数を求める関数*/

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

/*ルジャンドル多項式の係数を求める関数*/
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

int main(){
    double a;
    int b;
    double c;
    printf("Input variable:");
    a=5;   /*後でscanfに修正する*/
    for(b=0;b<=a;b++){
        c = legendre(a,b);   /*forループでa_(n,0)からa_(n,n)まで求めて出力します*/
        printf("%lf\n", c);
    }
    return 0;
}

