#include <cstdio>
#include <cmath>
#include <time.h>
#include <vector>
#include <random>

/*
更新者：武者野
最終更新：1/24 3:30
・fの値を新たなディレクトリ"data"に出力するようにした
・写真の枚数を別ファイル"picture_number"に出力
*/

// 0:固定 1:ランダム
#define INITIAL 1
using vd = std::vector<double>;
using vvd = std::vector<vd>;

/******************************計算条件******************************/
const int nx = 50+1;
const int ny = 50+1;
const double Lx = 1.0;
const double Ly = 1.0;
const double dx = Lx/double(nx-1);
const double dy = Ly/double(ny-1);

const double kappa = 1.0;
//計算の安定性を決めるファクターμ, mu > 0.25 だと計算が爆発する
const double mu = 0.20;

/*******************************************************************/
/**********************一様乱数の定義*************************/
//非決定的な乱数生成機、これを使ってシードをランダムに選ぶ
std::random_device rnddev;
//シード値の指定に使う乱数
std::mt19937 mt(rnddev());
//一様乱数
std::uniform_int_distribution<int> rndnx(0, nx-1);
std::uniform_int_distribution<int> rndny(0, ny-1);
std::uniform_real_distribution<double> rndd(0.0, 1.0);
/************************************************************/

//初期状態を決定
void initial(vvd &f);
//その時刻における f[jy][jx] の値をファイルとターミナルにアウトプット
void output(vvd &f, double t, FILE *data_fp);
//f[jy][jx] をもとにワンタイムステップ後の状態 fn[jy][jx] を計算。
//fn は fn[0][jx] など（つまり境界）は更新されないことに注意
void diffusion(vvd &f, vvd &fn, double dt);
//fn に境界条件を課す
void boundary(vvd &fn);



int main(){

  int icnt = 0;
  double dt, t = 0.0;
  vvd f(ny,vd(nx, 0.0)), fn(ny, vd(nx, 0.0));
  FILE *data_fp, *picnum_fp;//二つ目は写真の枚数を出力するファイル
  clock_t start_t, end_t;
  start_t = time(NULL);

  // ランダム変数のシードは時刻から取る、つまり毎回違うシード
  srand((unsigned)time(NULL));

  data_fp = fopen("data/diffusion.txt", "w");
  picnum_fp = fopen("data/picture_number.txt", "w");
  //printf("NX:%d NY:%d\nk:%f mu:%f\n", nx, ny, kappa, mu);
  fprintf(data_fp,"%d %d\n%f %f\n", nx, ny, kappa, mu);

  initial(f);

  dt = mu * fmin(dx*dx,dy*dy)/kappa;

  do{
    output(f, t, data_fp);
    diffusion(f, fn, dt);
    boundary(fn);

    f = fn;

    t += dt;
    icnt++;
  } while (t < 0.02 + 1e-9);

  //写真の枚数を出力するで～
  printf("number of pictures:%d\n", icnt);
  fprintf(picnum_fp,"%d", icnt);

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
      f[rndny(mt)][rndnx(mt)] = 1.0;
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
  for(int jy = 0; jy < ny; jy++) {
    for(int jx = 0; jx < nx; jx++) {
      if(jx < nx-1){
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
  for(int jy=0 ; jy < ny; jy++) fn[jy][0] = 0.0;
  for(int jy=0 ; jy < ny; jy++) fn[jy][nx-1] = 0.0;
  for(int jx=0 ; jx < nx; jx++) fn[0][jx] = 0.0;
  for(int jx=0 ; jx < nx; jx++) fn[ny-1][jx] = 0.0;
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

