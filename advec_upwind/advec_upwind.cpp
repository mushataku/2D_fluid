#include <cstdio>
#include <cmath>
#include <time.h>
#include <vector>
#include <random>

// 0:fixed 1:random
#define INITIAL 0
// 0:fixed 1:periodic
#define BOUNDARY 1
// 0:upwind 1:central difference
#define METHOD 1

using vd = std::vector<double>;
using vvd = std::vector<vd>;

/******************************計算条件******************************/
const int nx = 50+3;
const int ny = 50+3;
const double Lx = 1.0;
const double Ly = 1.0;
const double dx = Lx/double(nx-3);
const double dy = Ly/double(ny-3);

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

void advection_upwind(vvd &f, vvd &fn, double u, double v, double dt);
void advection_central(vvd &f, vvd &fn, double u, double v, double dt);

int main(){
  double dt, t = 0.0, u = 1.0, v = 1.0;
  int icnt = 0, ocnt = 0; //写真の枚数を数える
  vvd f(ny,vd(nx, 0.0)), fn(ny, vd(nx, 0.0));
  FILE *data_fp, *picnum_fp;//二つ目は写真の枚数を出力するファイル
  clock_t start_t, end_t;
  start_t = time(NULL);

  // ランダム変数のシードは時刻から取る、つまり毎回違うシード
  srand((unsigned)time(NULL));

  data_fp = fopen("data/advection_upwind.txt", "w");
  picnum_fp = fopen("data/picture_number.txt", "w");
  //printf("NX:%d NY:%d\nk:%f mu:%f\n", nx, ny, kappa, mu);
  fprintf(data_fp,"%d %d\n%f %f\n", nx-2, ny-2, kappa, mu);

  initial(f);
  boundary(f);
  
  dt = 0.2 * fmin(dx/fabs(u),dy/fabs(u));
  printf("dt:%f\n", dt);

  do{
    if(icnt%1 == 0){
      output(f, t, data_fp);
      ocnt++;
    }
    if(METHOD == 0)
      advection_upwind(f, fn, u, v, dt);
    if(METHOD == 1)
      advection_central(f, fn, u, v, dt);
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
  for(int jy = 1; jy < ny-1; jy++) {
    for(int jx = 1; jx < nx-1; jx++) {
      if(jx < nx-2){
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
    for(int jy=0 ; jy < ny; jy++) fn[jy][0] = fn[jy][nx-2];
    for(int jy=0 ; jy < ny; jy++) fn[jy][nx-1] = fn[jy][1];
    for(int jx=0 ; jx < nx; jx++) fn[0][jx] = fn[ny-2][jx];
    for(int jx=0 ; jx < nx; jx++) fn[ny-1][jx] = fn[1][jx];    
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

void advection_upwind(vvd &f, vvd &fn, double u, double v, double dt){

  double cflx = u*dt/dx, cfly = v*dt/dy;

  for(int jy = 1; jy < ny-1; jy++) {
    for(int jx = 1; jx < nx-1; jx++) {
      fn[jy][jx] = f[jy][jx];

      if(u > 0.0)
        fn[jy][jx] += -cflx*(f[jy][jx] - f[jy][jx-1]);
      else
        fn[jy][jx] += -cflx*(f[jy][jx+1] - f[jy][jx]);

      if (u > 0.0)
        fn[jy][jx] += -cflx*(f[jy][jx] - f[jy-1][jx]);
      else
        fn[jy][jx] += -cflx*(f[jy+1][jx] - f[jy][jx]);
    }
  }
}

void advection_central(vvd &f, vvd &fn, double u, double v, double dt){

  double cflx = u*dt/dx, cfly = v*dt/dy;

  for(int jy = 1; jy < ny-1; jy++) {
    for(int jx = 1; jx < nx-1; jx++) {
      fn[jy][jx] = f[jy][jx];

      fn[jy][jx] += -0.5*cflx*(f[jy][jx+1] - f[jy][jx-1]);
      fn[jy][jx] += -0.5*cfly*(f[jy+1][jx] - f[jy-1][jx]);
    }
  }
}
