#include <vector>
using vd = std::vector<double>;
using vvd = std::vector<vd>;

/******************************計算条件******************************/
const int nx = 200+5;
const int ny = 200+5;
const double Lx = 1.0;
const double Ly = 1.0;
const double dx = Lx/double(nx-5);
const double dy = Ly/double(ny-5);
// 計算終了時間
const double ENDTIME = 2.0;
/*********************************************************************/

/******************************関数******************************/
// 初期状態を決定
void initial(vvd &u, vvd &v);
// un, vn に境界条件を課す
void boundary(vvd &un, vvd &vn);

// Burgers eq を解く
void x_advection(vvd &f, vvd &fn, vvd &u, double dt);
void y_advection(vvd &f, vvd &fn, vvd &v, double dt);
void diffusion(vvd &f, vvd &fn, double dt);

// 仮の流速 u*, v* からpoisson 方程式の右辺 s を求める
void get_s(vvd &u, vvd &v, vvd &s, double dt);
// poisson eq を解く
int poisson2d(vvd &f, vvd &s);
// 圧力から速度を補正
void correction(vvd &u, vvd &v, vvd &p, double dt);

void output_all(/*u, v, rot など*/);
/****************************************************************/

int main(){
  vvd u(ny,vd(nx, 0.0)), un(ny, vd(nx, 0.0));
  vvd v(ny,vd(nx, 0.0)), vn(ny, vd(nx, 0.0));
  vvd p(ny,vd(nx, 0.0)), s(ny, vd(nx, 0.0));

  double t, dt; //dt は適宜決定
  initial(u, v);

  do{

    {// x方向に移流
      x_advection(u, un, u, dt);
      x_advection(v, vn, u, dt);
      boundary(un, vn);
      u = un; v = vn;
    }

    {// y方向に移流
      y_advection(u, un, v, dt);
      y_advection(v, vn, v, dt);
      boundary(un, vn);
      u = un; v = vn;
    }

    {// 拡散
      diffusion(u, un, dt);
      diffusion(v, vn, dt);
      boundary(un, vn);
      u = un; v = vn;
    }

    {// poisson 方程式と結合して速度の修正
      get_s(u, v, s, dt);
      poisson2d(p, s);
      correction(u, v, p, dt);
      boundary(u, v);
    }

    t += dt;
    output_all(/*u, v, rot など*/);
  } while (t < ENDTIME);

  return 0;
}