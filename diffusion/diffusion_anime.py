#最終更新： 2020/2/13 16:50

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython.display import HTML

###################### CONFIG ##################
# 動画保存方法の config。 使いたい保存方法を１にする
MP4       = 0
GIF       = 1
HTML_SHOW = 0
PLT       = 0
###################### CONFIG ##################

################### PARAMETER ##################
TITLE = '2D diffusion'

picnum_file = open('data/picture_number.txt')
FRAMES = int(picnum_file.readline())

f = open('data/diffusion.txt')
NX, NY = map(int, f.readline().split())
kappa, mu = map(float, f.readline().split())
################### PARAMETER ##################


fig = plt.figure(figsize=(6, 6))
fig.subplots_adjust(left=0.20)
ax = fig.add_subplot(111)
ax.set_xlim(0, NX-1)
ax.set_ylim(0, NY-1)
ax.set_xlabel('X', fontsize=16)
ax.set_ylabel('y', fontsize=16)
ax.tick_params(labelsize=14)
ax.set_title(TITLE, fontsize=20)
ax.text(0.02, 0.035, r"$\mu$="+str(mu)+"\n"+r"$\kappa$="+str(kappa),
          size=20, color="white", backgroundcolor='black', transform = ax.transAxes)


time_text = ax.text(0.02, 0.91, '', size=20, color="white",
          backgroundcolor='black', transform=ax.transAxes)


################### アニメの初期画像 ###################
time = float(f.readline())
time_text.set_text("time = " + str(time))

data=[]
for j in range(NY):
  tmp = list(map(float, f.readline().split()))
  data.append(tmp)

im = ax.imshow(data, animated=True, cmap='jet')
cbar = fig.colorbar(im, shrink=0.8)
cbar.ax.tick_params(labelsize=16)
######################################################


# 更新用関数
def update(cnt):
  # 最初は更新せずに初期の画像を出力させる
  if cnt == 0:
    return

  print("cnt:",cnt)  
  time = float(f.readline())
  
  # 二次元配列受け取り
  data=[]
  for j in range(NY):
    tmp = list(map(float, f.readline().split()))
    data.append(tmp)
  
  # 変更箇所だけset
  time_text.set_text("time = " + str(time))
  im.set_data(data)

ani = FuncAnimation(fig, update, repeat=True, interval=20, frames=FRAMES)

  # オプション frames 何枚の画像か，つまり関数を呼び出す回数

# mp4 ファイルとして保存
if(MP4 == 1):
  ani.save("animes/dif_anime.mp4", writer="ffmpeg")
# gif ファイルとして保存
if(GIF == 1):
  ani.save("animes/dif_anime.gif", writer="pillow")
# HTML上で表示
if(HTML_SHOW == 1):
  HTML(ani.to_jshtml())
if(PLT == 1):
  plt.show()

f.close
picnum_file.close
