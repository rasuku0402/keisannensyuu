#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

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
    int i;
    for(i=1;i<10;i++){
        int size = sizeof(double);
        double *x_str = malloc((i+2)*size);
        init(i+2, 0, x_str);
        x_str[0] = -1;
        x_str[1] = 1;
        double x_s = 0;
        shiftright(x_str, i+2);
        x_str[0] = x_s;

    }
    return 0;
}
