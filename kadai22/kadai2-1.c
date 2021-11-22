#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define N  9 /*正方行列のサイズ*/
#define KMAX  10000 /*最大反復回数*/
#define EPS  1E-10 /*最大誤差*/

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



/*比較関数(配列を昇順に並べ替えるqsort関数に必要)*/
int double_comp(const void *s1, const void *s2){
    const double a1 = *((double *)s1); 
    const double a2 = *((double *)s2);

    if(a1 < a2){
        return -1;
    }else if(a1 == a2){
        return 0;
    }else{
        return 1;
    }
}

/*最大値ノルムの計算(課題2-1-3で使う)*/
double vector_norm_max(double *a, int m, int n){
    int i,tmp;
    tmp = n-m+1 ;/*全要素数の計算*/

    for(i=m;i<=n;i++) a[i] = fabs(a[i]);
    qsort(a+m, tmp, sizeof(a[0]), double_comp);

    return a[n];

}

/*ヤコビ法* (課題2-1-3)*/
double *jacobi(double **a, double *b, double *x){
    double eps, *xn;
    int i,j,k=0;

    xn = dvector(1,N); 

    do{
        for(i=1;i<=N;i++){
            xn[i] = b[i]; 
            for(j=1;j<=N;j++){
                xn[i] -= a[i][j]*x[j]; 
            }
            xn[i] += a[i][i]*x[i];
            xn[i] /=a[i][i];
        }

        for(i=1;i<=N;i++) x[i] = xn[i] -x[i];
        eps = vector_norm_max(x,1,N); /*最大値ノルムの計算*/
        printf("eps%d = %f", k, eps);
        for(i=1;i<=N;i++) x[i] = xn[i];
        k++;

    }while(eps > EPS && k < KMAX);
    free_dvector(xn, 1);
    if(k==KMAX){
        printf("答えが見つかりませんでした\n");
        exit(1);
    }else{
        printf("反復回数は%d回です\n", k);
        return x;
    }
}


int main(){
    FILE *fin, *fout;
    double **a, *b, *x;
    int i;

    a = dmatrix(1,N,1,N);
    b = dvector(1,N);
    x = dvector(1,N);

    /*ファイルのオープン*/
    if((fin = fopen("input_sp.dat", "r"))==NULL){
        printf("ファイルが見つかりません：input_sp.dat \n");
        exit(1);
    }
    if((fout = fopen("output_sp.dat", "w")) == NULL){
        printf("ファイルが作成できません：output_sp.dat \n");
        exit(1);
    }

    input_matrix(a,'A',fin,fout);  /*行列Aの入出力*/
    input_vector(b,'b',fin,fout); /*ベクトルbの入出力*/
    input_vector(x,'x',fin,fout); /*初期ベクトルx_0の入力*/

    x = jacobi(a,b,x);

    fprintf(fout, "Ax = bの解は次の通りです\n");
    for(i=1;i<=N;i++){
        fprintf(fout, "%f\n", x[i]);
    }
    fclose(fin); fclose(fout);
    free_dmatrix(a,1,N,1,N);
    free_dvector(b,1);
    free_dvector(x,1);

    return 0;
}