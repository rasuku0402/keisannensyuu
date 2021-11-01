#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double x){
    return (63/8)*x*x*x*x*x - (35/4)*x*x*x + (15/8)*x;
}

double df(double x){
    return (315/8)*x*x*x*x - (105/4)*x*x + 15/8 ;
}

int main(){
    double x_st1 = 0.85;
    double x_st2 = 0.5;
    double x_st3 = 0.1;
    double x_st4 = -0.60;
    double x_st5 = -0.85;
    double esp = 100;

    double x;

    while(esp > 1e-9){
        x = x_st1;
        x_st1 = x_st1 - f(x_st1)/df(x_st1);
        esp = fabs(x_st1 - x);
    }
    printf("%lf\n", x_st1);
    esp = 100;
    while(esp > 1e-9){
        x = x_st2;
        x_st2 = x_st2 - f(x_st2)/df(x_st2);
        esp = fabs(x_st2 - x);
    }
    printf("%lf\n", x_st2);
    esp = 100;
    while(esp > 1e-9){
        x = x_st3;
        x_st3 = x_st3 - f(x_st3)/df(x_st3);
        esp = fabs(x_st3 - x);
        printf("esp = %lf\n", esp);
    }
    printf("%lf\n", x_st3);
    esp = 100;
    while(esp > 1e-9){
        x = x_st4;
        x_st4 = x_st4 - f(x_st4)/df(x_st4);
        esp = fabs(x_st4 - x);
    }
    printf("%lf\n", x_st4);
    esp = 100;
    while(esp > 1e-9){
        x = x_st5;
        x_st5 = x_st5 - f(x_st5)/df(x_st5);
        esp = fabs(x_st5 - x);
    }
    printf("%lf\n", x_st5);

    return 0;
}