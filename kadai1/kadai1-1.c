#include <stdio.h>
#include <math.h>

double Combi(double n, int r){

  int i;
  double c = 1.0;  

  if((1 <= r) && (r < n)){
    for(i=r;i>=1;i--){
      c *= ((n + 1.0 - (double)i) / (double)i);
    }
  }
  
  return c;
}

int main()
{
    double a;
    int b;
    double c;
    printf("Input variable:");
    scanf("%lf%d", &a,&b);
    c = Combi(a,b);
    printf("%f\n", c);
    return 0;
}