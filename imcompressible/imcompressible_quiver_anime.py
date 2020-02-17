#最終更新： 2020/2/15 16:00

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize  # Normalizeをimport
from IPython.display import HTML

###################### CONFIG ##################
# 動画保存方法の config, 使いたい保存方法を１にする
MP4       = 0
GIF       = 1
HTML_SHOW = 0
PLT       = 0

# 圧力と発散を表示させるかどうか、1なら表示
P_DIV = 0
###################### CONFIG ##################

################### PARAMETER ##################
TITLE = '2D Navier-Stokes equation'

# quiver で描画するときに領域を何区間に分割するか（あんまり多いと見辛いので

# NX-1. Ny-1 の約数じゃないと壊れる
QNX = 20
QNY = 20

cfp = open('data/condition.txt')

NX, NY, FRAMES = map(int, cfp.readline().split())
Lx, Ly, Re, mu = map(float, cfp.readline().split())
cfp.close
################### PARAMETER ##################

####################################描画のための関数####################################
def get_data(f):
  time = float(f.readline())
  data=[]
  for _ in range(NY):
    tmp = list(map(float, f.readline().split()))
    data.append(tmp)
  return data, time

def init_im(ax, f, title):
  data, _ = get_data(f)
  ax.tick_params(labelsize=14)
  ax.set_title(title, fontsize=20)
  if(title == "divergence"):
    MAX = 0.1
  elif(title == "rotation"):
    MAX = 20
  elif(title == "pressure"):
    MAX = 0.5
  im = ax.imshow(data, extent=(0,Lx,0,Ly), origin="lower", animated=True, cmap='jet', norm=Normalize(vmin=-MAX, vmax=MAX))
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cbar = fig.colorbar(im, cax=cax)
  cbar.ax.tick_params(labelsize=12)

  return im

def im_reset(im, f):
  data, time = get_data(f)
  im.set_data(data)
  return time

def get_UV():
  _ = float(f_u.readline())
  U = []
  for _ in range(NY):
    tmp = list(map(float, f_u.readline().split()))
    tmp = tmp[::NX//QNX]
    U.append(tmp)
  
  U = U[::NY//QNY]

  _ = float(f_v.readline())
  V = []
  for _ in range(NY):
    tmp = list(map(float, f_v.readline().split()))
    tmp = tmp[::NX//QNX]
    V.append(tmp)
  V = V[::NY//QNY]


  return U, V


def init_quiver(ax):
  X, Y = np.meshgrid(np.linspace(0, Lx, QNX+1), np.linspace(0, Ly, QNY+1))
  U, V = get_UV()
  C = np.hypot(U, V)
  ax.set_title("flow vector", fontsize=20)
  ax.set_aspect("equal")
  ax.tick_params(labelsize=14)
  im = ax.quiver(X, Y, U, V, C, scale_units="width", cmap="jet")

  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cbar = fig.colorbar(im, cax=cax)
  cbar.ax.tick_params(labelsize=12)

  return im

def quiver_reset(im):
  U, V = get_UV()
  C = np.hypot(U, V)
  im.set_UVC(U, V, C)

######################################################################################


f_u = open('data/u.txt')
f_v = open('data/v.txt')
f_div = open('data/div.txt')
f_rot = open('data/rot.txt')
f_p = open('data/p.txt')

if(P_DIV == 0):
  fig = plt.figure(figsize=(10, 6))
else:
  fig = plt.figure(figsize=(10, 10))

fig.suptitle(TITLE, fontsize=20)
#fig.subplots_adjust(bottom=0.10)

if(P_DIV == 1):
  ax_q = fig.add_subplot(221)
  ax_rot = fig.add_subplot(222)
  ax_p = fig.add_subplot(223)
  ax_div = fig.add_subplot(224)
else:
  ax_q = fig.add_subplot(121)
  ax_rot = fig.add_subplot(122)

fig.text(0, 0.01, "Re="+str(Re),
          backgroundcolor="black",color="white", size=20)


time_text = fig.text(0.02, 0.98, '', size=20, color="white", horizontalalignment='left',
            verticalalignment='top', backgroundcolor='black')

################### アニメの初期画像 ###################
# im_uの初期化
im_q = init_quiver(ax_q)
im_rot = init_im(ax_rot, f_rot, "rotation")

if(P_DIV == 1):
  im_p = init_im(ax_p, f_p, "pressure")
  im_div = init_im(ax_div, f_div, "divergence")

time_text.set_text("time = 0.00")

######################################################
plt.tight_layout()


# 更新用関数
def update(cnt):
  # 最初は更新せずに初期の画像を出力させる
  if cnt == 0:
    return

  print("cnt:",cnt)
  ######################################################

  # 変更箇所だけreset
  quiver_reset(im_q)
  time = im_reset(im_rot, f_rot)
  if(P_DIV == 1):
    im_reset(im_p, f_p)
    im_reset(im_div, f_div)

  time_text.set_text("time = %.2f"%time)
ani = FuncAnimation(fig, update, repeat=True, interval=20, frames=FRAMES)

  # オプション frames 何枚の画像か，つまり関数を呼び出す回数

# mp4 ファイルとして保存
if(MP4 == 1):
  ani.save("animes/imcompressible_anime.mp4", writer="ffmpeg", fps=20)
# gif ファイルとして保存
if(GIF == 1):
  ani.save("animes/imcompressible_anime.gif", writer="pillow", fps=20)
# HTML上で表示
if(HTML_SHOW == 1):
  HTML(ani.to_jshtml())
if(PLT == 1):
  plt.show()

f_u.close
f_v.close
f_div.close
f_rot.close
f_p.close
