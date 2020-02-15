#最終更新： 2020/2/15 13:00

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import FuncAnimation
from IPython.display import HTML

###################### CONFIG ##################
# 動画保存方法の config, 使いたい保存方法を１にする
MP4       = 0
GIF       = 1
HTML_SHOW = 0
PLT       = 0
###################### CONFIG ##################

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
  im = ax.imshow(data, extent=(0,Lx,0,Ly), origin="lower", animated=True, cmap='jet')
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)

  cbar = fig.colorbar(im, cax=cax)
  cbar.ax.tick_params(labelsize=12)

  return im

def im_reset(im, f):
  data, time = get_data(f)
  im.set_data(data)
  return time
######################################################################################


################### PARAMETER ##################
TITLE = '2D burgers'

cfp = open('data/condition.txt')

NX, NY ,FRAMES= map(int, cfp.readline().split())
Lx, Ly, Re, mu = map(float, cfp.readline().split())
cfp.close
################### PARAMETER ##################

f_u = open('data/u.txt')
f_v = open('data/v.txt')
f_div = open('data/div.txt')
f_rot = open('data/rot.txt')

fig = plt.figure(figsize=(10, 10))
fig.suptitle(TITLE, fontsize=20)
#fig.subplots_adjust(bottom=0.10)

ax_u = fig.add_subplot(221)
ax_v = fig.add_subplot(222)
ax_div = fig.add_subplot(223)
ax_rot = fig.add_subplot(224)

fig.text(0, 0.01, r"$\mu$="+ str(mu) +", Re="+str(Re),
          backgroundcolor="black",color="white", size=20)


time_text = fig.text(0.02, 0.98, '', size=20, color="white", horizontalalignment='left',
            verticalalignment='top', backgroundcolor='black')

################### アニメの初期画像 ###################
# im_uの初期化
im_u = init_im(ax_u, f_u, "u")
im_v = init_im(ax_v, f_v, "v")
im_div = init_im(ax_div, f_div, "divergence")
im_rot = init_im(ax_rot, f_rot, "rotation")

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
  time = im_reset(im_u, f_u)
  im_reset(im_v, f_v)
  im_reset(im_div, f_div)
  im_reset(im_rot, f_rot)
  time_text.set_text("time = %.2f"%time)
ani = FuncAnimation(fig, update, repeat=True, interval=20, frames=FRAMES)

  # オプション frames 何枚の画像か，つまり関数を呼び出す回数

# mp4 ファイルとして保存
if(MP4 == 1):
  ani.save("animes/burgers_anime.mp4", writer="ffmpeg", fps=20)
# gif ファイルとして保存
if(GIF == 1):
  ani.save("animes/burgers_anime.gif", writer="pillow", fps=20)
# HTML上で表示
if(HTML_SHOW == 1):
  HTML(ani.to_jshtml())
if(PLT == 1):
  plt.show()

f_u.close
