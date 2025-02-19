#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

#define MAX_ROWS 50  // (csvのデータ数が50まで対応今回はN=30)最大行数の仮設定(ゆとりをもつ)
// 3次元のパラメータ空間で(μ1,μ2,k)=(at,bt,ct)にそって最尤推定を行う
//a,b,x1,x2,x3はデータから事前に計算した(尤度関数のPfaffianに必要なため)
double a = 1, b, c, t0=0.1;//b = 4.266792032852375
double x1 = 0.5697013370597034, x2 = 0.43029866294029656, x3 = 0.16508059376762182;

// 尤度関数Lを計算(2次元の場合) c*tでκを表す, colsは変数に次元rowsはデータサイズ
double calculate_function(double **x, int rows, int cols, double a, double b, double c, double t) {
    double mu[cols];
    double sum = 0.0;
    mu[0] = a * t;  // μ1 = a * t
    mu[1] = b * t;  // μ2 = b * t

    for (int h = 0; h < rows; h++) {
        double inner_product = 0.0;
        for (int i = 0; i < cols; i++) {
            inner_product += mu[i] * x[h][i];
        }
        sum += inner_product * inner_product;
    }

    return exp(- (c * t) / rows * sum);
}

//階乗計算
unsigned long long factorial(int n) {
    if (n == 1) {
        return 1;
    }
    return n * factorial(n - 1);
}

// ポッホハマー記号 (a)_n の計算
double pochhammer(double a, int n) {
    double result = 1.0;
    for (int i = 0; i < n; i++) {
        result *= (a + i);
    }
    return result;
}

// クンマーの合流超幾何級数 1F1(a,b,x) の計算（第4引数で項数指定）
double kummer_M(double a, double b, double x, int terms) {
    double sum = 1.0; // 第0項
    for (int n = 1; n < terms; n++) {
        double term = (pochhammer(a, n) / pochhammer(b, n)) * (pow(x, n) / factorial(n));
        sum += term;
    }
    return sum;
}

// Pfaffian 方程式の定義 データから事前に計算したx1,x2,x3を必要とする.
int pfaffian_eqs(double t, const double y[], double f[], void *params) {
    (void)params; // 未使用パラメータを明示的に無視

    double A = -3 * c * t*t * (a * a * x1 + b * b * x2 + 2 * a * b * x3);
    double B = 3 * c;
    double C = (3 * a * a * t + 3 * b * b * t)/2 +A;
    double D = (3 * a * a * c * t*t * t + 3 * b * b * c * t*t * t - 1) / t;

    f[0] = A * y[0] + B * y[1];  // f[0]は CL
    f[1] = C * y[0] + D * y[1];  // f[1]は d_k CL

    return GSL_SUCCESS;
}

int main(void) {
    double min_y0=100000;// CLの初期値を大きな値で初期化
    double best_t =0.0;
    double best_k =0.0;
    double best_c=0.0;
    double best_b=0.0;
    char fname[]="half_circle_data.csv" ;//単位円周上のデータ

    FILE *file = fopen(fname, "r");
    if (!file) {
        perror("ファイルを開けません");
        return 1;
    }

    // CSVの最大列数を仮設定
    double *x[MAX_ROWS]={NULL};
    int rows = 0, cols = 0;
    char line[1024];

    // ファイルを1行ずつ読み取る(half_Circle_data)
    while (fgets(line, sizeof(line), file) && rows < MAX_ROWS) {
        x[rows] = malloc(4 * sizeof(double)); // 最大4列まで確保(S^3球面まで対応)
        char *ptr = line;
        if (!x[rows]) {
            fprintf(stderr, "メモリ確保に失敗\n");
            return 1;
        }

        int col = 0;
        while (sscanf(ptr, "%lf", &x[rows][col]) == 1) {
            ptr = strchr(ptr, ',');//sscanfのためにカンマのポインタを1つ進める
            if (!ptr) break;
            ptr++;
            col++;
        }
        cols = col + 1;
        rows++;
    }
    fclose(file);

    //
    for (double B=0.1; B<10.0; B+=0.1){
    for (int C=1; C<=1000;C++){
    // ユーザー入力 パラメータ空間の探索曲線(at,bt,ct)=(μ1,μ2,k) kappaに関する傾きt\in(0.1,1.5)と仮定
    c=(double)C;
    b=(double)B;

    // 尤度関数計算
    double result = calculate_function(x, rows, cols, a, b, c, t0);
    printf("関数の結果: %lf\n", result);

    // 1F1(1/2,1,k) の計算
    double M1 = kummer_M(0.5, 1.0, c*t0, 20);
    
    // dk1F1(1/2,1,k)=1/2 * 1F1(3/2,2,k) の計算
    double M2 = (1.0 / 2.0) * kummer_M(1.5, 2.0, c*t0, 20);
    
    printf("M(1/2,1,%.3f) = %.15f\n", c*t0, M1);
    printf("(1/2) * M(3/2,2,%.3f) = %.15f\n", c*t0, M2);

    double CL=M1*result;
    double KCL=M2*result;

    //Pfaffianのベクトル値関数の初期値を計算
    printf("CL = %.15f\n",  CL);
    printf("dkCL = %.15f\n", KCL);

    gsl_odeiv2_system sys = {pfaffian_eqs, NULL, 2, NULL}; // ホロノミックランク2

    // ルンゲクッタ4次 (RK4) に変更
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rk4, 1e-2, 1e-6, 1e-6);

    double t = t0, t1 = 1.1;  // tの範囲
    double y[2] = {CL, KCL};  // ベクトル値関数の初期値をクンマーの合流超幾何級数から計算
    for (double ti = t ; ti <= t1; ti += 0.01) {  // ステップサイズ0.01固定
        if (y[0] < 0.1 && y[1] < 0.1) {
            printf("f[0] and f[1] are both less than 0.1 at t = %g. Ending simulation.\n", t);
            break;
        }
        //パラメータμに関して単位円周上付近でMLEを計算
        if(ti*ti*(1+b*b)<0.1 || 1.1<ti*ti*(1+b*b)){
            break;

        }
        int status = gsl_odeiv2_driver_apply(d, &t, ti, y);
        if (status != GSL_SUCCESS) {
            printf("Error: solver failed at t = %g\n", ti);
            break;
        }

        printf("t = %g, CL(t) = %g, dkCL(t) = %g, k=%g\n", t, y[0], y[1],c*t);
    }
    if (y[0]<min_y0){
        min_y0=y[0]; //尤度関数CLの最小値
        best_t=t;   //パラメータtの最小値
        best_k=c*t; //最尤推定量κ
        best_c=C;  //κ=ct
        best_b=b;  //μ2=bt, μ1=t
    }

    gsl_odeiv2_driver_free(d);
}
    }
    for (int i = 0; i < rows; i++) {
        free(x[i]);  //　最後にメモリ解法
    }

    printf("尤度関数の最小値 CL=%g, at t=%g, μ1=%g, μ2=%g, k=%g, a=%d, b=%g, c=%g", min_y0, best_t, best_t, best_t*best_b, best_k, 1, best_b, best_c);
    return 0;
}
