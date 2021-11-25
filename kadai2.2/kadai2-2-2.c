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
            fprintf(fout, "%5.11f\t", a[i][j]);
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
        fprintf(fout, "%5.11f\t", b[i]);
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

/*ベクトルの差を求める関数(a-b)*/
void vector_sum(double *a, double *b, double *c, int n){
    int i;
    for(i=1;i<=n;i++){
        c[i] = a[i] - b[i];
    }
}

/*ベクトルのl1ノルムを求める関数*/
double vector_norm1(double *a, int n){
    int i;
    double norm = 0.0;
    for(i = 1; i<=n;i++){
        norm += fabs(a[i]);
    }
    return norm;
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

/*ヤコビ法* (課題2-1-3)*/
double *jacobi(double **a, double *b, double *x, FILE *fin, FILE *fout){
    double eps, eps2,*xn, *c, *r, *ans,*e; /*cはAx_k= cであり、残差ノルムの計算に使う,rは残差ベクトル*/
    int i,j,k=0;
    xn = dvector(1,N); 
    c = dvector(1,N);
    r = dvector(1,N);
    e = dvector(1,N);
    ans = dvector(1,N); /*真の解を格納*/
    for(i=1;i<=N;i++){
        ans[i] = 1.0;
    }

    do{
        for(i=1;i<=N;i++){
            xn[i] = b[i]; 
            for(j=1;j<=N;j++){
                xn[i] -= a[i][j]*x[j]; 
            }
            xn[i] += a[i][i]*x[i];
            xn[i] /=a[i][i];
        }

        
        matrix_vector_product(a, xn, c);
        vector_sum(b, c, r,N); /*残差ベクトルrを求めた*/
        vector_sum(ans,xn,e,N);
        eps = vector_norm1(r,N);
        eps2 = vector_norm1(e,N); /*eps2は真の解と近似解の差のノルム*/

        /*eps = vector_norm_max(x,1,N);  最大値ノルムの計算(多分間違ってる)*/
        fprintf(fout, "r%d = %.11f  ", k, eps);
        fprintf(fout, "e%d = %.11f\n", k, eps2);
        for(i=1;i<=N;i++) x[i] = xn[i]; /* x_kをx_k+1に更新*/
        k++;

    }while(eps > EPS && k < KMAX);
    free_dvector(xn, 1);
    free_dvector(c,1);
    free_dvector(r,1);
    free_dvector(e,1);
    free_dvector(ans,1);
    if(k==KMAX){
        printf("答えが見つかりませんでした\n");
        exit(1);
    }else{
        printf("反復回数は%d回です\n", k);
        printf("残差ノルムは%.11f\n", eps);
        printf("真の解と近似解の差のノルムは%.11f\n",eps2);
        return x;
    }
}

/*Gauss-Seidel法*(課題2-1-3)*/

double *gauss_seidel(double **a, double *b, double *x, FILE *fin, FILE *fout){
    double eps,eps2,*xo, *c, *r,*ans,*e, s, t;
    int i,j,k=0;
    xo = dvector(1,N);
    c = dvector(1,N);
    r = dvector(1,N);
    ans = dvector(1,N); /*真の解を格納*/
    e = dvector(1,N);
    for(i=1;i<=N;i++){
        ans[i] = 1.0;
    }

    do{
        for(i=1;i<=N;i++){
            xo[i] = x[i]; /*xoがx_k. x_kにx_k+1を代入*/
        } 
        /*i=1の処理だけべつで行う*/

        t = 0.0;
        for(j=2;j<=N;j++)  t +=a[1][j]*xo[j];
        x[1] = (b[1] - t)/a[1][1];

        for(i=2;i<=N;i++){
            s = 0.0; t = 0.0;
            for(j=1;j<i;j++) s += a[i][j]*x[j];
            for(j = i+1; j<=N;j++) t += a[i][j]*xo[j];
            x[i] = (b[i]-s-t)/a[i][i];
        }
        matrix_vector_product(a, x, c);
        vector_sum(b, c, r, N); /*残差ベクトルrを求めた*/
        vector_sum(ans,x,e,N);
        eps = vector_norm1(r,N);
        eps2 = vector_norm1(e,N);
        fprintf(fout, "r%d = %.11f ", k, eps);
        fprintf(fout, "e%d = %.11f\n", k, eps2);

        k++;
    }while(eps> EPS && k < KMAX);

    free_dvector(xo, 1);
    free_dvector(c, 1);
    free_dvector(r,1);
    free_dvector(ans,1);
    free_dvector(e,1);

    if(k == KMAX){
        printf("答えが見つかりませんでした\n");
        exit(1);
    }else{
        printf("反復回数は%d回です\n", k);
        printf("残差ノルムは%.11f\n", eps);
        printf("真の解と近似解の差のノルムは%.11f\n",eps2);
        return x;
    }
} 

/*SOR法 (課題2-1-3)*/
double *sor(double **a, double *b, double *x, double w, FILE *fin, FILE *fout){
    double eps,eps2,*xo,s,t,*c,*r,*ans,*e;
    int i,j,k=0;
    xo = dvector(1,N);
    c = dvector(1,N);
    r = dvector(1,N);
    ans = dvector(1,N); /*真の解を格納*/
    e = dvector(1,N);
    for(i=1;i<=N;i++){
        ans[i] = 1.0;
    }

    do{
        for(i=1;i<=N;i++) xo[i] = x[i];
        t = 0.0;
        for(j=2;j<=N;j++) t += a[1][j]*xo[j];
        x[1] = (b[1] -t)/a[1][1];
        for(i=2;i<=N;i++){
            s=0.0; t=0.0;
            for(j=1;j<i;j++) s += a[i][j]*x[j];
            for(j=i+1;j<=N;j++) t += a[i][j]*xo[j];
            x[i] = (b[i] -s-t)/a[i][i];
        }
        /*ここまではgauss法と同じ*/

        /*SOR法*/
        for(i=1;i<=N;i++){
            x[i] = xo[i] + w*(x[i] - xo[i]);
        }
        matrix_vector_product(a, x, c);
        vector_sum(b, c, r, N); /*残差ベクトルrを求めた*/
        vector_sum(ans, x, e,N);
        eps = vector_norm1(r,N);
        eps2 = vector_norm1(e,N);
        fprintf(fout, "r%d = %.11f  ", k, eps);
        fprintf(fout, "e%d = %.11f\n", k,eps2);
        k++;  
    }while(eps> EPS && k < KMAX);

    free_dvector(xo, 1);
    free_dvector(c, 1);
    free_dvector(r,1);
    free_dvector(ans,1);
    free_dvector(e,1);
    if(k == KMAX){
        printf("答えが見つかりませんでした\n");
        exit(1);
    }else{
        printf("反復回数は%d回です\n", k);
        printf("残差ノルムは%.11f\n", eps);
        printf("真の解と近似解の差のノルムは%.11f\n",eps2);

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
    if((fin = fopen("input_kadai2-2.dat", "r"))==NULL){
        printf("ファイルが見つかりません：input_sp.dat \n");
        exit(1);
    }
    if((fout = fopen("output_kadai2-2-2_jacobi.dat", "w")) == NULL){
        printf("ファイルが作成できません：output_sp.dat \n");
        exit(1);
    }

    input_matrix(a,'A',fin,fout);  /*行列Aの入出力*/
    input_vector(b,'b',fin,fout); /*ベクトルbの入出力*/
    input_vector(x,'x',fin,fout); /*初期ベクトルx_0の入力*/

    x = jacobi(a,b,x,fin,fout);

    fprintf(fout, "ヤコビによるAx = bの解は次の通りです\n");
    for(i=1;i<=N;i++){
        fprintf(fout, "%.11f\n", x[i]);
    }

    
    fclose(fin); fclose(fout); /*ヤコビ法終わり*/

    /*ここからSOR.ファイルのオープン*/
    if((fin = fopen("input_kadai2-2.dat", "r"))==NULL){
        printf("ファイルが見つかりません：input_sp.dat \n");
        exit(1);
    }
    if((fout = fopen("output_kadai2-2-2_SOR.dat", "w")) == NULL){
        printf("ファイルが作成できません：output_gauss.dat \n");
        exit(1);
    }

    input_matrix(a,'A',fin,fout);  /*行列Aの入出力*/
    input_vector(b,'b',fin,fout); /*ベクトルbの入出力*/
    input_vector(x,'x',fin,fout); /*初期ベクトルx_0の入力*/

    x = sor(a,b,x,1.79,fin,fout);
    fprintf(fout, "SOR法によるAx = bの解は次の通りです\n");
    for(i=1;i<=N;i++){
        fprintf(fout, "%.11f\n", x[i]);
    }
    fclose(fin); fclose(fout); /*SOR法終わり*/



    /*ここからGauss-Seidel.ファイルのオープン*/

    if((fin = fopen("input_kadai2-2.dat", "r"))==NULL){
        printf("ファイルが見つかりません：input_sp.dat \n");
        exit(1);
    }
    if((fout = fopen("output_kadai2-2-2_gauss.dat", "w")) == NULL){
        printf("ファイルが作成できません：output_gauss.dat \n");
        exit(1);
    }

    input_matrix(a,'A',fin,fout);  /*行列Aの入出力*/
    input_vector(b,'b',fin,fout); /*ベクトルbの入出力*/
    input_vector(x,'x',fin,fout); /*初期ベクトルx_0の入力*/

    x = gauss_seidel(a,b,x,fin,fout);
    fprintf(fout, "Gauss-seidel法によるAx = bの解は次の通りです\n");
    for(i=1;i<=N;i++){
        fprintf(fout, "%.11f\n", x[i]);
    }

    fclose(fin); fclose(fout); /*gauss法終わり*/

    

    free_dmatrix(a,1,N,1,N);
    free_dvector(b,1);
    free_dvector(x,1);

    return 0;
}