#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

/*行列とベクトルの積の計算を実行するサブルーチン (課題2-1-1)*/
void matrix_vector_product(double **a, double *b, double *c){
    int i,j;
    for(i=1;i<=N;i++){ 
        c[i] = 0;
        for(j=1;j<=N;j++){
            c[i] = a[i][j]*b[j] + c[i];
        }
    }
}

/*LU分解を行うサブルーチン(課題2-1-2)*/

void lu_decom(double **a, int *p){
    int i,j,k,ip;
    double alpha, tmp;
    double amax, eps = pow(2.0, -50.0); /*eps,amaxは入力された行列が正則であるかを判断するための変数*/

    for(k=1;k<=N-1;k++){
        amax = fabs(a[k][k]); ip = k;
        for(i=k+1;i<=N;i++){
            if(fabs(a[i][k])> amax){
                amax = fabs(a[i][k]); ip = i;
            }
        }
    /*正則性の判定*/
    if(amax < eps) printf("入力された行列は正則ではありません\n");
    /*ip を配列pに保存*/
    p[k] = ip;
    /*行交換*/
    if(ip != k){
        for(j=k;j<=N;j++){
            tmp = a[k][j]; a[k][j] = a[ip][j]; a[ip][j] = tmp;
        }
    }
    /*前進消去*/
    for(i=k+1;i<=N;i++){
        alpha = -a[i][k]/a[k][k];
        a[i][k] = alpha;
        for(j=k+1;j<=N;j++){
            a[i][j] = a[i][j] + alpha*a[k][j];
        }
    }
    }

}

/*LU分解により、連立1次方程式を解くサブルーチン(課題2-1-2)*/

double *lu_solve(double **a, double *b, int *p){
    int i,j,k;
    double tmp;

    /*右辺の行交換*/
    for(k=1;k<=N-1;k++){
        tmp=b[k]; b[k] = b[p[k]]; b[p[k]]= tmp;
        /*前進代入*/
        for(i=k+1;i<=N;i++){
            b[i] = b[i] +a[i][k]*b[k];
        }
    }

    /*後退代入*/
    b[N] = b[N]/a[N][N];
    for(k=N-1;k>=1;k--){
        tmp = b[k];
        for(j=k+1;j<=N;j++){
            tmp = tmp -a[k][j]*b[j];
        }
        b[k] = tmp/a[k][k];
    }
    /*連立方程式の解を出力*/
    printf("LU分解によって連立方程式Ax=bの解を求めた結果が以下の通り.\n");
    for(i=1;i<=N;i++){  
        printf("x[%d] = %.9f\n", i,b[i]);
    }
    return b;
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
    double **a, *α, *b;
    int i,j;
    FILE *fin, *fout;
    int p[N];

    a = dmatrix(1,N,1,N);
    α = dvector(1,N);
    b = dvector(1,N);

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
    input_vector(α,'a',fin,fout); /*ベクトルbの入出力*/

    fclose(fin);
    fclose(fout);


    /*行列の積の計算*/
    matrix_vector_product(a,α,b); /*aα=bを計算*/
    /*結果の表示*/
    printf("Aα=bで求めたベクトルbは以下の通りです \n");
    for(j=1;j<=N;j++){
        printf("b[%d] = %f\n",j,b[j]);
    }
    printf("\n");

    lu_decom(a, p); /*行列AをLU分解*/
    b = lu_solve(a,b,p); /*Ax = bを出力 (ここではbが解として出力されている)*/


    
    /*行列,ベクトル領域の開放*/
    free_dmatrix(a,1,N,1,N);
    free_dvector(α,1);
    free_dvector(b,1);

    return 0;

}
                    