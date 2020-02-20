/*
最終更新：2020/2/14
更新者：武者野

[内容]
・７章のコードを再現
・s(x,y)に対するpossion方程式の解f(x,y)を算出
・境界条件はf(x,y)=0
・s_func()に好きな関数形(境界条件は満たす)を入れることで任意の関数について poisson eq が解ける
・ここでは解析解の分かっているs(x,y)に大して計算を行い求めることで動作を確認
*/

#define _USE_MATH_DEFINES
#include <cstdio>
#include <cmath>
#include <time.h>
#include <vector>

using vd = std::vector<double>;
using vvd = std::vector<vd>;

/******************************計算条件******************************/
const int nx = 256+1;
const int ny = 256+1;
const double Lx = 1.0;
const double Ly = 1.0;
const double dx = Lx/double(nx-1);
const double dy = Ly/double(ny-1);
// SOR 法を何回反復したら残差を check するか
const int CHECK_INTERVAL = 100;
const double omega = 1.8;
const double eps = 1.0e-8;
/*******************************************************************/

// SOR法で1ステップだけf更新。誤差|f - fn|の最大値を返す
double sor(vvd &f, vvd &s);
// sor をたくさん呼び出して poisson 方程式を解く。返り値は反復回数
int poisson2d(vvd &f, vvd &s);
// 残差の最大値を計算。反復中に呼び出すことで精度の上昇を見れる。
double residual(vvd &f, vvd &s);
// 計算で get したfを解析解と比較
double error_f(vvd &f);
//解くべきs(x,y)の関数形
double s_func(double x, double y);
//上のs(x,y)の元での解析式f(x,y)
double f_analitic(double x, double y);

int main(int argc, char *argv[]){
  int it;
  double x, y;
  vvd s(ny,vd(nx, 0.0)), f(ny, vd(nx, 0.0));

  // 配列の初期化
  for(int jy = 0; jy < ny; jy++){
    for(int jx = 0; jx < nx; jx++){
      // mesh -> 座標 に変換
      x = dx*(double)jx; y = dy*(double)jy;
      f[jy][jx] = 0.0;
      s[jy][jx] = s_func(x,y);
    }
  }
  // check 用
  // for(int jy = 0; jy < ny; jy++) {
  //   for(int jx = 0; jx < nx; jx++) {
  //     if(jx < nx-1) printf("%f ", s[jy][jx]);
  //     else printf("%f\n", s[jy][jx]);
  //   }
  // }

  it = poisson2d(f,s);
  printf("\n%3dx%3d error=%9.3e iteration=%d\n",nx-1,ny-1,error_f(f),it);
  
  //check用
  // for(int jy = 0; jy < ny; jy++) {
  //   for(int jx = 0; jx < nx; jx++) {
  //     if(jx < nx-1) printf("%f ", f[jy][jx]);
  //     else printf("%f\n", f[jy][jx]);
  //   }
  // }
  // printf("%d\n");
  // for(int jy = 0; jy < ny; jy++) {
  //   for(int jx = 0; jx < nx; jx++) {
  //     x = dx*(double)jx; y = dy*(double)jy;
  //     if(jx < nx-1) printf("%f ", f_analitic(y,x));
  //     else printf("%f\n", f_analitic(y,x));
  //   }
  // }

  return 0;
}

double s_func(double x, double y){
  double kx = 2.0*M_PI, ky = 2.0*M_PI;
  return -(kx*kx + ky*ky)*f_analitic(x, y);
}

double f_analitic(double x, double y){
  double kx = 2.0*M_PI, ky = 2.0*M_PI;
  return sin(kx*x)*sin(ky*y);
}

double sor(vvd &f, vvd &s){
  double fn, err = 0.0;
  for (int jy = 1; jy < ny - 1; jy++){
    for (int jx = 1; jx < nx - 1; jx++){
      fn = ((f[jy][jx + 1] + f[jy][jx - 1]) / (dx * dx) 
         + (f[jy + 1][jx] + f[jy - 1][jx]) / (dy * dy) - s[jy][jx]) 
         * 0.5 * dx * dx * dy * dy / (dx * dx + dy * dy);
      err = fmax(fabs(fn - f[jy][jx]), err);
      f[jy][jx] = (1.0 - omega) * f[jy][jx] + omega * fn;
    }
  }

  return err;
}

int poisson2d(vvd &f, vvd &s){
  int icnt = 0, imax = 99999;
  double err;

  // 反復
  while(icnt++ < imax){
    //fの更新と更新前後の差の最大値確認
    err = sor(f, s);
    //更新してもほとんど変わらないようなら終了
    if(eps > err) return icnt;
    if(icnt%CHECK_INTERVAL == 0){
      double resi = residual(f, s);
      printf("reidual:%f\n", resi);
    }
    // printf("err:%f\n", err);
  }

  // imax 回反復しても収束しないなら-1を返して終了
  return -1;
}

double residual(vvd &f, vvd &s){
  double res, rmax = 0.0;
  //各格子点の(d^2f/dx^2 + d^2f/dy^2 - s)を計算
  for (int jy = 1; jy < ny - 1; jy++){
    for (int jx = 1; jx < nx - 1; jx++){
      res = (f[jy][jx + 1] - 2.0 * f[jy][jx] + f[jy][jx - 1]) / (dx * dx) 
          + (f[jy + 1][jx] - 2.0 * f[jy][jx] + f[jy - 1][jx]) / (dy * dy)
          - s[jy][jx];
      rmax = fmax(res, rmax);
    }
  }

  return rmax;
}

double error_f(vvd &f){
  double x, y, err = 0.0;

  for (int jy = 1; jy < ny - 1; jy++){
    for (int jx = 1; jx < nx - 1; jx++){
      //座標に変換
      x = dx * (double)jx;
      y = dy * (double)jy;

      err += fabs(f[jy][jx] - f_analitic(x,y));
    }
  }

  //err を格子点の数で平均
  return err / (double)((nx - 2) * (ny - 2));
}