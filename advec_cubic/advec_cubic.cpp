#define _USE_MATH_DEFINES
#include <cstdio>
#include <cmath>
#include <time.h>
#include <vector>
#include <random>

// 0:fixed 1:random 2:point 3:sin
#define INITIAL 3
// 0:fixed 1:periodic
#define BOUNDARY 1
// 0:upwind 1:central difference
#define METHOD 0

using vd = std::vector<double>;
using vvd = std::vector<vd>;

/******************************計算条件******************************/
const int nx = 100+4;
const int ny = 100+4;
const double Lx = 1.0;
const double Ly = 1.0;
const double dx = Lx/double(nx-4);
const double dy = Ly/double(ny-4);

const double kappa = 1.0;
//計算の安定性を決めるファクターμ, mu > 0.25 だと計算が爆発する
const double mu = 0.20;

/*******************************************************************/

//初期状態を決定
void initial(vvd &f);
//その時刻における f[jy][jx] の値をファイルとターミナルにアウトプット
void output(vvd &f, double t, FILE *data_fp);
//f[jy][jx] をもとにワンタイムステップ後の状態 fn[jy][jx] を計算。
//fn は fn[0][jx] など（つまり境界）は更新されないことに注意
void diffusion(vvd &f, vvd &fn, double dt);
//fn に境界条件を課す
void boundary(vvd &fn);

void x_advection(vvd &f, vvd &fn, double u, double dt);
void y_advection(vvd &f, vvd &fn, double v, double dt);

int main(){
  double dt, t = 0.0, u = 1.0, v = 1.0;
  int icnt = 0, ocnt = 0; //写真の枚数を数える
  vvd f(ny,vd(nx, 0.0)), fn(ny, vd(nx, 0.0));
  FILE *data_fp, *picnum_fp;//二つ目は写真の枚数を出力するファイル
  clock_t start_t, end_t;
  start_t = time(NULL);

  // ランダム変数のシードは時刻から取る、つまり毎回違うシード
  srand((unsigned)time(NULL));

  data_fp = fopen("data/advection_cubic.txt", "w");
  picnum_fp = fopen("data/picture_number.txt", "w");
  //printf("NX:%d NY:%d\nk:%f mu:%f\n", nx, ny, kappa, mu);
  fprintf(data_fp,"%d %d\n%f %f %f %f\n", nx-4, ny-4, Lx, Ly, kappa, mu);

  initial(f);
  boundary(f);
  
  dt = 0.2 * fmin(dx/fabs(u),dy/fabs(v));
  printf("dt:%f\n", dt);

  do{
    if(icnt%2 == 0){
      output(f, t, data_fp);
      ocnt++;
    }
    x_advection(f, fn, u, dt);
    boundary(fn);
    f = fn;

    y_advection(f, fn, v, dt);
    boundary(fn);
    f = fn;

    t += dt;
    icnt++;
  } while (t < 2.0 + 1e-9);

  //写真の枚数を出力するで～
  printf("number of pictures:%d\n", ocnt);
  fprintf(picnum_fp,"%d", ocnt);

  fclose(picnum_fp);
  fclose(data_fp);

  end_t = time(NULL);
  printf("This calculatioin took %ld second \n", end_t - start_t);

  return 0;
}

void initial(vvd &f){
  if(INITIAL == 0){
    for(int jy = 0; jy < ny; jy++) {
      for(int jx = 0; jx < nx; jx++) {
        if(0.3*nx < jx && jx < 0.7*nx && 0.3*ny < jy && jy < 0.7*ny){
          f[jy][jx] = 1.0;
        }
        // else if(0.8*nx < jx && jx < 0.9*nx && 0.8*ny < jy && jy < 0.9*ny){
        //   f[jy][jx] = 5.0;
        // }
        else{
          f[jy][jx] = 0.0;
        }
      }
    }
  }
  if(INITIAL == 1){
    //全部の点の10%くらいをランダムに選んで大きい値を持たせる
    for(int i = 0; i < nx*ny*0.8; i++) {
      f[rand()%ny][rand()%nx] = 1.0;
    }
  }

  if(INITIAL == 2){
    for(int jy = 0; jy < ny; jy++) {
      for(int jx = 0; jx < nx; jx++) {
        if(0.3*nx < jx && jx < 0.7*nx && 0.3*ny < jy && jy < 0.7*ny){
          f[jy][jx] = 1e5;
        }
        // else if(0.8*nx < jx && jx < 0.9*nx && 0.8*ny < jy && jy < 0.9*ny){
        //   f[jy][jx] = 5.0;
        // }
        else{
          f[jy][jx] = 1.0;
        }
      }
    }
  }
  if(INITIAL == 3){
    double x, y, kx = 2.0*M_PI, ky = 2.0*M_PI;
    for(int jy = 0; jy < ny; jy++) {
      for(int jx = 0; jx < nx; jx++) {
        x = dx*(double)(jx-2);
        y = dy*(double)(jy-2);
        f[jy][jx] = sin(kx*x)*sin(ky*y);
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
  //printf("t:%f\n", t);
  fprintf(data_fp,"%f\n", t);

  //端っこの境界条件のためのダミーの格子点は出力しない
  for(int jy = 2; jy < ny-2; jy++) {
    for(int jx = 2; jx < nx-2; jx++) {
      if(jx < nx-3){
        //printf("%.2f ", f[jy][jx]);
        fprintf(data_fp, "%f ", f[jy][jx]);
      }
      else{
        //printf("%.2f\n", f[jy][jx]);
        fprintf(data_fp, "%f\n", f[jy][jx]);
      }
    }
  }
}

void diffusion(vvd &f, vvd &fn, double dt){
  //jy は通常 0 から ny-1 まで走るが、境界条件であるところのjy = 0, jy = ny-1 は更新しない。jxも同様
  for(int jy = 1; jy < ny-1; jy++) {
    for(int jx = 1; jx < nx-1; jx++) {
      fn[jy][jx] = f[jy][jx] + kappa * dt * 
      ( (f[jy][jx+1] - 2.0*f[jy][jx] + f[jy][jx-1])/dx/dx
         + (f[jy+1][jx] - 2.0*f[jy][jx] + f[jy-1][jx])/dy/dy
      );
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

void x_advection(vvd &f, vvd &fn, double u, double dt){

  double a,b,c,z;

  for (int jy = 2; jy < ny - 2; jy++){
    for (int jx = 2; jx < nx - 2; jx++){
      if (u > 0.0){
        a = (f[jy][jx+1] - 3.0*f[jy][jx] + 3.0*f[jy][jx-1] - f[jy][jx-2]) / (6.0*dx*dx*dx);
        b = (f[jy][jx+1] - 2.0*f[jy][jx] + f[jy][jx-1]) / (2.0*dx*dx);
        c = (2.0*f[jy][jx+1] + 3.0*f[jy][jx] - 6.0*f[jy][jx-1] + f[jy][jx-2]) / (6.0*dx);
      }
      else{
        a = (f[jy][jx+2] - 3.0*f[jy][jx+1] + 3.0*f[jy][jx] - f[jy][jx-1]) / (6.0*dx*dx*dx);
        b = (f[jy][jx+1] - 2.0*f[jy][jx] + f[jy][jx-1]) / (2.0*dx*dx);
        c = (-f[jy][jx+2] + 6.0*f[jy][jx+1] - 3.0*f[jy][jx] - 2.0*f[jy][jx-1]) / (6.0*dx);
      }
      z = -u*dt;
      fn[jy][jx] = a*z*z*z + b*z*z + c*z + f[jy][jx];
    }
  }
}

void y_advection(vvd &f, vvd &fn, double v, double dt){

  double a,b,c,z;

  for (int jy = 2; jy < ny - 2; jy++){
    for (int jx = 2; jx < nx - 2; jx++){
      if (v > 0.0){
        a = (f[jy+1][jx] - 3.0*f[jy][jx] + 3.0*f[jy-1][jx] - f[jy-2][jx]) / (6.0*dy*dy*dy);
        b = (f[jy+1][jx] - 2.0*f[jy][jx] + f[jy-1][jx]) / (2.0*dy*dy);
        c = (2.0*f[jy+1][jx]+ 3.0*f[jy][jx] - 6.0*f[jy-1][jx] + f[jy-2][jx]) / (6.0*dy);
      }
      else{
        a = (f[jy+2][jx] - 3.0*f[jy+1][jx] + 3.0*f[jy][jx] - f[jy-1][jx]) / (6.0*dy*dy*dy);
        b = (f[jy+1][jx] - 2.0*f[jy][jx] + f[jy-1][jx]) / (2.0*dy*dy);
        c = (-f[jy+2][jx] + 6.0*f[jy+1][jx] - 3.0*f[jy][jx] - 2.0*f[jy-1][jx]) / (6.0*dy);
      }
      z = -v*dt;
      fn[jy][jx] = a*z*z*z + b*z*z + c*z + f[jy][jx];
    }
  }
}
