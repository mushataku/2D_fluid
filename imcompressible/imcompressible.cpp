#define _USE_MATH_DEFINES
#include <cstdio>
#include <cmath>
#include <time.h>
#include <vector>
#include <random>

using vd = std::vector<double>;
using vvd = std::vector<vd>;


/******************************CONFIG******************************/
// 0:square 1:random 2:point 3:sin
#define INITIAL 0
// 0:fixed 1:periodic
#define BOUNDARY 1
// 0:upwind 1:central difference
#define METHOD 0
// 0:terminal に出力しない 1:出力する
#define TERMINAL 0
/*******************************************************************/


/******************************計算条件******************************/
const int nx = 100+4;
const int ny = 100+4;
const double Lx = 1.0;
const double Ly = 1.0;
const double dx = Lx/double(nx-4);
const double dy = Ly/double(ny-4);

const double Re = 2000.0;
//計算の安定性を決めるファクターμ, mu > 0.25 だと計算が爆発する
const double mu = 0.20;
/*******************************************************************/

/******************************出力枚数の調整用******************************/
/*
時刻 0 ~ endtime の間のプロファイルを時間 DT ごとに出力させるための変数達
TIME に(0, DT, 2DT, ..., endtime) をsetする
*/
const double T_EPS = 1.0e-10;
const double DT = 0.02;
const double ENDTIME = 5.0;
vd TIME;
//出力時刻をset
void TIME_set();
/*************************************************************************/

/*************************poisson eq 関連の関数・定数*************************/
// 仮の流速 u*, v* からpoisson 方程式の右辺 s を求める
void get_s(vvd &u, vvd &v, vvd &s, double dt);
// SOR法で1ステップだけf更新。誤差|f - fn|の最大値を返す
double sor(vvd &f, vvd &s);
// sor をたくさん呼び出して poisson 方程式を解く。返り値は反復回数
int poisson2d(vvd &f, vvd &s);
// 残差の最大値を計算。反復中に呼び出すことで精度の上昇を見れる。
double residual(vvd &f, vvd &s);

// 圧力から速度を補正
void correction(vvd &u, vvd &v, vvd &p, double dt);

const double omega = 1.8;
const double P_EPS = 1.0e-5;
/****************************************************************************/

/******************************burgersを解く関数******************************/
// 初期状態を決定
void initial(vvd &u, vvd &v);
// その時刻における f[jy][jx] の値をファイルとターミナルにアウトプット
void output(vvd &f, double t, FILE *data_fp);
// fn に境界条件を課す
void boundary(vvd &fn);

void diffusion(vvd &f, vvd &fn, double dt);
void x_advection(vvd &f, vvd &fn, vvd &u, double dt);
void y_advection(vvd &f, vvd &fn, vvd &v, double dt);
void rotation(vvd &u, vvd &v, vvd &rot);
void divergence(vvd &u, vvd &v, vvd &div);
/****************************************************************************/


/**************************ファイル**************************/
FILE *condition_fp = fopen("data/condition.txt","w");
FILE *u_fp = fopen("data/u.txt", "w");
FILE *v_fp = fopen("data/v.txt", "w");
FILE *div_fp = fopen("data/div.txt", "w");
FILE *rot_fp = fopen("data/rot.txt", "w");
FILE *p_fp = fopen("data/p.txt", "w");

/***********************************************************/

int main(){
  double dt, t = 0.0;
  int ti = 0, icnt = 0; //TIMEのindex
  int iteration;

  vvd u(ny,vd(nx, 0.0)), un(ny, vd(nx, 0.0));
  vvd v(ny,vd(nx, 0.0)), vn(ny, vd(nx, 0.0));
  vvd rot(ny,vd(nx, 0.0)), div(ny, vd(nx, 0.0));
  vvd p(ny,vd(nx, 0.0)), s(ny, vd(nx, 0.0));

  clock_t start_t, end_t;
  start_t = time(NULL);

  // ランダム変数のシードは時刻から取る、つまり毎回違うシード
  srand((unsigned)time(NULL));


  initial(u, v);
  boundary(u); boundary(v);
  TIME_set();

  dt = fmin(0.2*fmin(dx,dy)/1.0, mu*fmin(dx*dx, dy*dy)*Re);
  
  printf("NX:%d NY:%d\nRe:%f mu:%f\n", nx, ny, Re, mu);
  printf("dt:%.10f\n", dt);
  //ti = 3;
  do{
    if(ti < TIME.size() && t > TIME[ti] - T_EPS){
      if(ti == 0){
        get_s(u, v, s, dt);
        poisson2d(p, s);
      }
      divergence(u, v, div);
      rotation(u, v, rot);
      output(u, t, u_fp);
      output(v, t, v_fp);
      output(div, t, div_fp);
      output(rot, t, rot_fp);
      output(p, t, p_fp);
      ti++; icnt++;
    }

    x_advection(u, un, u, dt);
    x_advection(v, vn, u, dt);
    boundary(un); boundary(vn);
    u = un; v = vn;

    y_advection(u, un, v, dt);
    y_advection(v, vn, v, dt);
    boundary(un); boundary(vn);
    u = un; v = vn;

    diffusion(u, un, dt);
    diffusion(v, vn, dt);
    boundary(un); boundary(vn);
    u = un; v = vn;

    get_s(u, v, s, dt);
    iteration = poisson2d(p, s);
    if(iteration == -1){
      printf("poisson does not converged in time = %f\n", t);
      return 0;
    }
    correction(u, v, p, dt);
    boundary(u); boundary(v);

    printf("time:%f\n", t);

    t += dt;
  } while (t < ENDTIME + DT);

  printf("number of pictures:%d\n", icnt);
  fprintf(condition_fp, "%d %d %d\n%f %f %f %f\n", nx-4, ny-4, icnt, Lx, Ly, Re, mu);

  fclose(condition_fp);
  fclose(u_fp);
  fclose(v_fp);
  fclose(div_fp);
  fclose(rot_fp);
  fclose(p_fp);

  end_t = time(NULL);
  printf("This calculatioin took %ld second \n", end_t - start_t);

  return 0;
}

void initial(vvd &u, vvd &v){
  if(INITIAL == 0){
    for(int jy = 0; jy < ny; jy++) {
      for(int jx = 0; jx < nx; jx++) {
        if(0.3*nx < jx && jx < 0.7*nx && 0.3*ny < jy && jy < 0.7*ny){
          u[jy][jx] = 1.0;
          v[jy][jx] = 1.0;
        }/*
        else{
          u[jy][jx] = -1.0;
          v[jy][jx] = 0.0;
        }*/
      }
    }
  }
  if(INITIAL == 1){
    //全部の点の50%くらいをランダムに選んで大きい値を持たせる
    for(int i = 0; i < nx*ny*0.5; i++) {
      u[rand()%ny][rand()%nx] = 1.0;
      v[rand()%ny][rand()%nx] = 1.0;
    }
  }

  if(INITIAL == 2){
    u[ny/2][nx/2] = 1.0;
  }
  if(INITIAL == 3){
    double x, y, kx = 2.0*M_PI, ky = 2.0*M_PI;
    for(int jy = 0; jy < ny; jy++) {
      for(int jx = 0; jx < nx; jx++) {
        x = dx*(double)(jx-2);
        y = dy*(double)(jy-2);

        u[jy][jx] = -cos(kx*x)*sin(ky*y)/kx;
        v[jy][jx] = sin(kx*x)*cos(ky*y)/ky;
        
        x += 0.3; y+= 0.7;
        u[jy][jx] = -0.6*cos(2.0*kx*x)*sin(2.0*ky*y)/kx;
        v[jy][jx] = 0.6*sin(2.0*kx*x)*cos(2.0*ky*y)/ky;
        
      }
    }
  }
  
  //check用
  // for(int jy = 0; jy < ny; jy++) {
  //   for(int jx = 0; jx < nx; jx++) {
  //     if(jx < nx-1){
  //       printf("%.1f ", f[jy][jx]);
  //     }
  //     else{
  //       printf("%.1f\n", f[jy][jx]);
  //     }
  //   }
  // }

  return;
}

void output(vvd &f, double t, FILE *data_fp){

  if(TERMINAL == 1) printf("t:%f\n", t);
  fprintf(data_fp,"%f\n", t);

  //端っこの境界条件のためのダミーの格子点は出力しない
  for(int jy = 2; jy < ny-2; jy++) {
    for(int jx = 2; jx < nx-2; jx++) {
      if(jx < nx-3){
        if(TERMINAL == 1) printf("%.2f ", f[jy][jx]);
        fprintf(data_fp, "%f ", f[jy][jx]);
      }
      else{
        if(TERMINAL == 1) printf("%.2f\n", f[jy][jx]);
        fprintf(data_fp, "%f\n", f[jy][jx]);
      }
    }
  }
}

void diffusion(vvd &f, vvd &fn, double dt){
  //境界は更新しない
  for(int jy = 2; jy < ny-2; jy++) {
    for(int jx = 2; jx < nx-2; jx++) {
      fn[jy][jx] = f[jy][jx] + dt * 
      ( (f[jy][jx+1] - 2.0*f[jy][jx] + f[jy][jx-1])/dx/dx
         + (f[jy+1][jx] - 2.0*f[jy][jx] + f[jy-1][jx])/dy/dy
      )/Re;
    }
  }
  return;
}

void boundary(vvd &fn){
  if(BOUNDARY == 0){
    for(int jy=0 ; jy < ny; jy++) fn[jy][0] = 0.0;
    for(int jy=0 ; jy < ny; jy++) fn[jy][nx-1] = 0.0;
    for(int jx=0 ; jx < nx; jx++) fn[0][jx] = 0.0;
    for(int jx=0 ; jx < nx; jx++) fn[ny-1][jx] = 0.0;
  }
  if(BOUNDARY == 1){
    for(int jy=0 ; jy < ny; jy++) fn[jy][nx-1] = fn[jy][3];
    for(int jy=0 ; jy < ny; jy++) fn[jy][nx-2] = fn[jy][2];
    for(int jy=0 ; jy < ny; jy++) fn[jy][1] = fn[jy][nx-3];
    for(int jy=0 ; jy < ny; jy++) fn[jy][0] = fn[jy][nx-4];

    for(int jx=0 ; jx < nx; jx++) fn[ny-1][jx] = fn[3][jx];
    for(int jx=0 ; jx < nx; jx++) fn[ny-2][jx] = fn[2][jx];
    for(int jx=0 ; jx < nx; jx++) fn[1][jx] = fn[ny-3][jx];
    for(int jx=0 ; jx < nx; jx++) fn[0][jx] = fn[ny-4][jx];
  }
  //check 用
  // printf("boundary chech\n");
  // for(int jy = 0; jy < ny; jy++) {
  //   for(int jx = 0; jx < nx; jx++) {
  //     if(jx < nx-1){
  //       printf("%.1f ", f[jy][jx]);
  //     }
  //     else{
  //       printf("%.1f\n", f[jy][jx]);
  //     }
  //   }
  // }
  return;
}

void x_advection(vvd &f, vvd &fn, vvd &u, double dt){

  double a,b,c,z;

  for (int jy = 2; jy < ny - 2; jy++){
    for (int jx = 2; jx < nx - 2; jx++){
      if (u[jy][jx] > 0.0){
        a = (f[jy][jx+1] - 3.0*f[jy][jx] + 3.0*f[jy][jx-1] - f[jy][jx-2]) / (6.0*dx*dx*dx);
        b = (f[jy][jx+1] - 2.0*f[jy][jx] + f[jy][jx-1]) / (2.0*dx*dx);
        c = (2.0*f[jy][jx+1] + 3.0*f[jy][jx] - 6.0*f[jy][jx-1] + f[jy][jx-2]) / (6.0*dx);
      }
      else{
        a = (f[jy][jx+2] - 3.0*f[jy][jx+1] + 3.0*f[jy][jx] - f[jy][jx-1]) / (6.0*dx*dx*dx);
        b = (f[jy][jx+1] - 2.0*f[jy][jx] + f[jy][jx-1]) / (2.0*dx*dx);
        c = (-f[jy][jx+2] + 6.0*f[jy][jx+1] - 3.0*f[jy][jx] - 2.0*f[jy][jx-1]) / (6.0*dx);
      }
      z = -u[jy][jx]*dt;
      fn[jy][jx] = a*z*z*z + b*z*z + c*z + f[jy][jx];
    }
  }
}

void y_advection(vvd &f, vvd &fn, vvd &v, double dt){

  double a,b,c,z;

  for (int jy = 2; jy < ny - 2; jy++){
    for (int jx = 2; jx < nx - 2; jx++){
      if (v[jy][jx] > 0.0){
        a = (f[jy+1][jx] - 3.0*f[jy][jx] + 3.0*f[jy-1][jx] - f[jy-2][jx]) / (6.0*dy*dy*dy);
        b = (f[jy+1][jx] - 2.0*f[jy][jx] + f[jy-1][jx]) / (2.0*dy*dy);
        c = (2.0*f[jy+1][jx]+ 3.0*f[jy][jx] - 6.0*f[jy-1][jx] + f[jy-2][jx]) / (6.0*dy);
      }
      else{
        a = (f[jy+2][jx] - 3.0*f[jy+1][jx] + 3.0*f[jy][jx] - f[jy-1][jx]) / (6.0*dy*dy*dy);
        b = (f[jy+1][jx] - 2.0*f[jy][jx] + f[jy-1][jx]) / (2.0*dy*dy);
        c = (-f[jy+2][jx] + 6.0*f[jy+1][jx] - 3.0*f[jy][jx] - 2.0*f[jy-1][jx]) / (6.0*dy);
      }
      z = -v[jy][jx]*dt;
      fn[jy][jx] = a*z*z*z + b*z*z + c*z + f[jy][jx];
    }
  }
}

void TIME_set(){
  double tmp = 0.0;
  while(tmp < ENDTIME + T_EPS){
    TIME.push_back(tmp);
    tmp += DT;
  }
}

void rotation(vvd &u, vvd &v, vvd &rot){
  for(int jy = 2; jy < ny-2; jy++) {
    for(int jx = 2; jx < nx-2; jx++) {
      rot[jy][jx] = 0.5*(v[jy+1][jx+1] - v[jy+1][jx]+ v[jy][jx+1] - v[jy][jx])/dx
                  - 0.5*(u[jy+1][jx+1] - u[jy][jx+1] + u[jy+1][jx] - u[jy][jx])/dy;
    }
  }
}

void divergence(vvd &u, vvd &v, vvd &div){
  for (int jy = 2; jy < ny-2; jy++){
    for (int jx = 2; jx < nx-2; jx++){
      div[jy][jx] = 0.5*(u[jy+1][jx+1] - u[jy+1][jx] + u[jy][jx+1] - u[jy][jx])/dx
                  + 0.5*(v[jy+1][jx+1] - v[jy][jx+1] + v[jy+1][jx] - v[jy][jx])/dy;
    }
  }
}

void get_s(vvd &u, vvd &v, vvd &s, double dt){
  for (int jy = 2; jy < ny-2; jy++){
    for (int jx = 2; jx < nx-2; jx++){
      s[jy][jx] = 0.5*(u[jy+1][jx+1] - u[jy+1][jx] + u[jy][jx+1] - u[jy][jx])/dx
                  + 0.5*(v[jy+1][jx+1] - v[jy][jx+1] + v[jy+1][jx] - v[jy][jx])/dy;
      s[jy][jx] /= dt;
    }
  }
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
    if(P_EPS > err) return icnt;
  }

  // imax 回反復しても収束しないなら-1を返して終了
  return -1;
}

void correction(vvd &u, vvd &v, vvd &p, double dt){
  for (int j = 2; j < ny-2; j++){
    for (int i = 2; i < nx-2; i++){
      u[j][i] += -0.5*(p[j][i] - p[j][i-1] + p[j-1][i] - p[j-1][i-1])/dx*dt;
      v[j][i] += -0.5*(p[j][i] - p[j-1][i] + p[j][i-1] - p[j-1][i-1])/dy*dt;
    }
  }
}