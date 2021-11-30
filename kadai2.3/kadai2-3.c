#include <stdio.h>
#include <math.h>
#include <stdlib.h>


#define N  9 /*正方行列のサイズ*/
#define KMAX  10000 /*最大反復回数*/
#define EPS  1E-10/*最大誤差*/

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

/*ベクトルa[m,...n], b[m,...,n]の内積を計算する*/
double inner_product (int m,int n, double *a, double *b){
    int i;
    double s = 0.0;
    for(i = m;i<=n;i++) s += a[i]*b[i];
    return s;
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
    
    return b;
}



/*べき乗法により絶対値最大の固有値を求めるプログラム*/

double power_method(double **a, double *x, FILE *fout){
    int i,k=0; /*kは反復回数*/
    double eps = pow(10.0, -10.0);/*eps = 10^{-8}*/
    double v2, v2s, *v, lambda;
    v = dvector(1,N);

    do{
        matrix_vector_product(a,x,v);
        lambda = inner_product(1,N,v,x);
        v2 = inner_product(1,N,v,v);
        v2s = sqrt(v2);
        for(i=1;i<=N;i++) x[i] = v[i]/v2s;
        ++k;
    }while(fabs(v2 - lambda*lambda) >= EPS);

    fprintf(fout, "反復回数は%d\n", k);
    fprintf(fout, "絶対値最大固有値lambdaは%.11f\n", lambda);
    fprintf(fout, "これに対応する固有ベクトルは次の通りです\n");

    for(i=1;i<=N;i++){
        fprintf(fout,"v[%d] = %.11f\n", i, x[i]);
    }
    free_dvector(v,1);
    return 1.0/lambda;
}


double power_method2(double **a, double *x, FILE *fout){
    int i,k=0; /*kは反復回数*/
    double eps = pow(10.0, -10.0);/*eps = 10^{-8}*/
    double v2, v2s, *v, lambda;
    v = dvector(1,N);

    for(i=1;i<=N;i++) x[i] = 0.5;

    do{
        matrix_vector_product(a,x,v);
        lambda = inner_product(1,N,v,x);
        v2 = inner_product(1,N,v,v);
        v2s = sqrt(v2);
        for(i=1;i<=N;i++) x[i] = v[i]/v2s;
        ++k;
    }while(fabs(v2 - lambda*lambda) >= EPS);

    fprintf(fout, "反復回数は%d\n", k);
    fprintf(fout, "絶対値最小固有値lambdaは%.11f\n", 1.0/lambda);
    fprintf(fout, "これに対応する固有ベクトルは次の通りです\n");

    for(i=1;i<=N;i++){
        fprintf(fout,"v[%d] = %.11f\n", i, x[i]);
    }
    free_dvector(v,1);

    return 1.0/lambda;
}

/*行列の逆行列を求めるサブルーチン*/
void inverse_matrix(double **a, double **a_inv){
    int i,j,k;
    double buf;
    //単位行列を作る
    for(i=1;i<=N;i++){
       for(j=1;j<=N;j++){
       a_inv[i][j]=(i==j)?1.0:0.0;
       }
    }
    //掃き出し法
    for(i=1;i<=N;i++){
        buf=1/a[i][i];
        for(j=1;j<=N;j++){
            a[i][j]*=buf;
            a_inv[i][j]*=buf;
        }
    
        for(j=1;j<=N;j++){
            if(i!=j){
                buf=a[j][i];
                for(k=1;k<=N;k++){
                    a[j][k]-=a[i][k]*buf;
                    a_inv[j][k]-=a_inv[i][k]*buf;
                }
            }
        }
    }



}




/*逆べき乗法により絶対値最小の固有値を求める*/

// void inverse(double **a, double *x,int *p, FILE *fout){
//     int i,k=0;
//     double lambda, v2,v2s,*v,x2,x2s;
//     v = dvector(1,N);


//     do{
//         lu_decom(a, p);  /*行列aをLU分解*/
//         v = lu_solve(a, x, p);  /*lu分解法により、Ax_(k+1) = x_(k) を解く*/
//         lambda = inner_product(1,N,v,x);
//         v2 = inner_product(1,N,v,v);
//         v2s = sqrt(v2);
//         for(i=1;i<=N;i++) x[i] = v[i]/v2s;
//         ++k;
        
//     }while(fabs(v2 - lambda*lambda) >= EPS && k < KMAX);

//     fprintf(fout, "反復回数は%d\n", k);
//     fprintf(fout, "絶対値最小固有値lambdaは%f\n", lambda);
//     fprintf(fout, "これに対応する固有ベクトルは次の通りです\n");

//     for(i=1;i<=N;i++){
//         fprintf(fout,"v[%d] = %f\n", i, x[i]);
//     }
//     free_dvector(v,1);
// }


// void inverse(double **ao, double **a, double eps, FILE *fout){
//     int i, j, k, p[N];
//     double v2, v2s, *v, *y, **lu, lambda, mu = 0.0, muo;

//     v=dvector(1,N); y=dvector(1,N);
//     lu = dmatrix(1,N,1,N);

//     for(i=1;i<=N;i++){
//         lambda = a[i][i]; /*近似固有値の代入*/
//         for(j=1;j<=N;j++) y[j] = 0.0; y[i] = 1.0; /*初期値設定*/

//         /*行列の作成及びLU分解*/
//         for(k=1;k<=N;k++){
//             for(j=1;j<=N;j++) lu[k][j] = ao[k][j];
//             lu[k][k] = ao[k][k] - lambda;
//         }
//         lu_decom(lu, p); /*lu分解*/

//         /*逆反復法*/
//         do{
//             muo = mu;
//             for(j=1;j<=N;j++) v[j] = y[j];
//             v = lu_solve(lu, v, p); /*固有ベクトルの計算*/
//             mu = inner_product(1,N,v,y); /*補正*/
//             v2 = inner_product(1,N,v,v);
//             v2s = sqrt(v2);
//             for(j = 1; j<=N; j++) y[j] = v[j]/v2s ;/*正規化*/
//         }while(fabs((mu-muo)/mu) >= eps);

//         /*結果の代入(固有ベクトルはaのi列に)*/
//         for(j=1;j<=N;j++) a[j][i] = y[j];

//     }
//     fprintf(fout, "絶対値最小固有値は%f\n", lambda);
//     fprintf(fout, "それに対応する固有ベクトルは以下\n");
//     for(i=1;i<=N;i++){
//         fprintf(fout,"v[%d] = %f\n", i, y[i]);
//     }

//     free_dmatrix(lu, 1,N,1,N); free_dvector(v,1); free_dvector(y, 1);

// }

int main(){
    FILE *fin, *fout;
    int i,j,p[N];
    double **a, *x,**ao,**a_inv; /*aoは行列保存用*/
    double lambda1,lambda2,lambda3,lambda4;
    double joukensuu1, joukennsuu2;
    a = dmatrix(1,N,1,N); ao = dmatrix(1,N,1,N);
    x = dvector(1,N); a_inv= dmatrix(1,N,1,N);

    /*ファイルのオープン(ここから1つ目の行列Aの最大固有値とその固有ベクトルを求める)*/
    if((fin = fopen("matrixA.dat", "r"))==NULL){
        printf("ファイルが見つかりません：input_sp.dat \n");
        exit(1);
    }
    if((fout = fopen("output_kadai2-3.dat", "w")) == NULL){
        printf("ファイルが作成できません：result_eigen_A.dat \n");
        exit(1);
    }
    input_matrix(a,'A',fin,fout);  /*行列Aの入出力*/
    input_vector(x,'x',fin,fout); /*ベクトルxの入力*/

    for(i=1; i<=N;i++) {
        for (j=1;j<=N; j++) ao[i][j] = a[i][j];
    }
    lambda1 = power_method(a,x,fout); /*べき乗法*/
    inverse_matrix(ao, a_inv); /*逆べき乗法*/
    lambda2 = power_method2(a_inv,x,fout);
    joukensuu1 = lambda1/lambda2;
    printf("行列Aの条件数は%.11f\n", joukensuu1);


    fclose(fin); fclose(fout);

    /*ファイルのオープン(ここから2つ目の行列Aの最大固有値を求める)*/
    if((fin = fopen("matrixA_hiruberut.dat", "r"))==NULL){
        printf("ファイルが見つかりません：input_sp.dat \n");
        exit(1);
    }
    if((fout = fopen("output_kadai2-3_hiruberut.dat", "w")) == NULL){
        printf("ファイルが作成できません：result_eigen_A.dat \n");
        exit(1);
    }

    a = dmatrix(1,N,1,N); ao = dmatrix(1,N,1,N);
    x = dvector(1,N);


    input_matrix(a,'A',fin,fout);  /*行列Aの入出力*/
    input_vector(x,'x',fin,fout); /*ベクトルxの入力*/

    for(i=1; i<=N;i++) {
        for (j=1;j<=N; j++) ao[i][j] = a[i][j];
    }

    lambda3 = power_method(a,x,fout);
    inverse_matrix(ao,a_inv);
    lambda4 = power_method2(a_inv,x,fout);

    joukennsuu2 = lambda3/lambda4;
    printf("ヒルベルト行列Aの条件数は%.11f\n", joukennsuu2);



    free_dvector(x,1);
    free_dmatrix(a,1,N,1,N);
    free_dmatrix(ao,1,N,1,N);

    fclose(fin); fclose(fout);


    

    return 0;
}