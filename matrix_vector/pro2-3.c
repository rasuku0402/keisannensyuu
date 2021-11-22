#include <stdlib.h>
#include <stdio.h>

#define N  9

/*行列の領域確保*/
double **dmatrix(int nr1, int nr2, int nl1, int nl2){
    int i;
    double **a;
    int nrow,ncol;

    nrow = nr2-nr1+1; /*行の数*/
    ncol = nl2-nl1+1; /*列の数*/

    /*行の確保*/

    if((a = (double **)malloc(nrow*sizeof(double *)))==NULL){
        printf("メモリが確保できません(行列a) \n");
        exit(1);
    }
    a = a -nr1; /*行をずらす*/

    /*列の確保*/
    for(i=nr1;i<=nr2; i++){
       a[i] =  (double *)malloc(ncol*sizeof(double));
    }

    for(i=nr1;i<=nr2;i++) {
        a[i] = a[i] -nl1;  /*列をずらす*/
    }

    return(a);
}

/*行列の領域解放*/

void free_dmatrix(double **a,int nr1, int nr2, int nl1, int nl2){
    int i;

    /*メモリの確保*/
    for(i=nr1;i<=nr2;i++) free((void*)(a[i]+nl1));
    free((void *)(a+nr1));
}

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

/*a[1...N][1...N]の入力*/

void input_matrix(double **a, char c, FILE *fin, FILE *fout){
    int i,j;

    fprintf(fout, "行列%cは次の通りです\n", c);
    for(i=1;i<=N;i++){
        for(j=1;j<=N;j++){
            fscanf(fin, "%lf", &a[i][j]);
            fprintf(fout, "%5.2f\t", a[i][j]);
        }
        fprintf(fout, "\n");
    }
}

/*b[1...N]の入力*/
void input_vector(double *b, char c, FILE *fin, FILE *fout){
    int i;

    fprintf(fout, "ベクトル%cは次の通りです\n", c);
    for(i=1;i<=N;i++){
        fscanf(fin, "%lf", &b[i]);
        fprintf(fout, "%5.2f\t", b[i]);
        fprintf(fout, "\n");
    }
}

int main(){
    FILE *fin, *fout;
    double **a, *b;
    
    /*行列a及びベクトルbのための領域確保*/
    a = dmatrix(1,N,1,N); /*行列a[1...N][1...N]*/
    b = dvector(1,N); /*ベクトル b[1...N]*/

    /*ファイルのオープン*/
    if((fin = fopen("input.dat", "r"))==NULL){
        printf("ファイルが見つかりません：input.dat \n");
        exit(1);
    }
    if((fout = fopen("output.dat", "w")) == NULL){
        printf("ファイルが作成できません：output.dat \n");
        exit(1);
    }
    input_matrix(a,'A',fin,fout);  /*行列Aの入出力*/
    input_vector(b,'b',fin,fout); /*ベクトルbの入出力*/

    fclose(fin);
    fclose(fout);

    /*領域解放*/
    free_dmatrix(a,1,N,1,N);
    free_dvector(b,1);

    return 0;

}