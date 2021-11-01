#include <stdio.h>
#include <math.h>
#include <string.h>

void swap(double *a, double *b){
    double ter=*a;
    *a = *b;
    *b = ter;
}
void bouble(int n,double *o){
    int i;
    int j;
    for (i = 1; i < n;i++){
        for (j = 0; j < n - 1; j++)
        {
            if (o[j] > o[j + 1]){
                swap(&o[j], &o[j + 1]);
            }
        }
    }
}

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
    double c[12]= {0};
    double d = 0;
    for(k=0;k<=n;k++){
        c[k] = legendre(n,k);
        d = d + c[k]*pow(x,k);
    }
    return d;
}

/* 二分法 初期値 x1<x2 と 誤差限界 eps を入力 */
double bisec(double x1, double x2, double eps,int n) {
    double c;
    while ((x2 - x1) >= eps) {
        c = (x1+x2)/2.0; /* 中点計算 */
        if ((P(x1, n)*P(c,n))> 0) {
             /* 同符号か判定 */
            x1 = c;
        } else {

            x2 = c;
        }
        
    }
    return c;
}

/*配列を右にひとつずらす関数*/

void shiftright(double *a,int size){
    int i;
    double t = a[size -1];
    for(i=size;i>0;i--){
        a[i] = a[i-1];
    }
    a[0] = t;
}

void init(int n,double x,double*o){
    int i;
    for (i = 0; i <= n;i++){
        o[i] = x;
    }
}


int main(){
    double c;
    double eps=1e-9;
    int i,j,k;
    double x_str[12] = {-1,1,2,2,2,2,2,2,2,2,2,2}; /*二分法の区間指定する配列*/
    for(i=1;i<=10;i++){
        double x_s[10] = {0}; /*二分法によって得られる解を格納する配列*/
        printf("n=%d\n", i);

        /*ルジャンドル多項式の解を計算する*/

        for(j=0;j<i;j++){
          x_s[j] = bisec(x_str[j],x_str[j+1],eps,i);
          printf("%.9lf\n", x_s[j]);

        }
        /*x_strを更新する*/
        /*もう一度初期化*/
        x_str[0] = -1;
        x_str[1] = 1;
        for(j=2;j<12;j++){
            x_str[j] = 2;
        }
        for(j =0;j<i;j++){
          shiftright(x_str,12);/*配列を右シフト*/
          x_str[0] = x_s[j];
        }       
        bouble(12, x_str);
    }

    return 0;


}