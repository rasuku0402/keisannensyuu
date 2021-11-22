#include <stdlib.h>
#include <stdio.h>

/*a[i]~a[j]の領域を確保*/
double *dvector(int i, int j){
    double *a;
    if ((a = malloc(((j-i+1)*sizeof(double)))) == NULL){
        printf("メモリが確保できません(from dvector) \n");
        exit(1);
    }
    return(a-i);
}

void free_dvector(double *a, int i){
    free((void*)(a+i));  
}
int main(){
    double *a, *b;
    int i;
    int N =5;

    printf("aの領域確保\n"); /*配列aの添え字は1~N*/
    a = dvector(1,N);
    for (i=1;i<=N;i++){
        a[i] = (double)i/20.0;
        printf("a[%d] = %f \n", i, a[i]);
    }
    free_dvector(a,1);
    printf("aを解放\n");

    return 0;

}