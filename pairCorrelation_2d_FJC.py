# Pair Correlation Function (2D)
# 2Dの粒子配置から対相関関数g(r)を計算し、プロットするコード
# これはRadial Distribution Function (RDF)とも呼ばれる

# オリジナルはAIによって生成された3D版を参考にしている。
# 元のコード: code/PairCorrelation/pairCorrelation_3d_by_AI.py

# 2D格子上の自由連結鎖の場合
# 格子上のランダムコイルだと格子上に規則的に配置された場合とあまり見分けがつかない

import random as rd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def calculate_rdf(positions, L, dr=0.1):
    """
    2D系の対相関関数 g(r) を計算する簡単な関数
    positions: (N, 2)の配列 (粒子の座標)
    L: ボックスサイズ (正方形)
    dr: ビンの幅
    """
    N = len(positions)
    max_r = L / 2
    bins = np.arange(0, max_r, dr)
    hist = np.zeros(len(bins)-1)
    
    count=0                                         # ちゃんとペア数を数えてるかのチェック用

    # 周期境界条件を考慮した距離の計算（簡略化のため全ペア計算）
    for i in range(N):
        for j in range(i + 1, N):
            # 最小イメージ法（周期的境界条件）
            diff = positions[i] - positions[j]
            diff = diff - L * np.round(diff / L)    # ここが周期的境界条件の適用部分だろうが、とりあえず認める
            r = np.linalg.norm(diff)                # ノルムの計算
            
#            count+=1                               # ちゃんとペア数を数えてるかのチェック用

            if r < max_r:                           # rを残さず、binに入れてしまう
                bin_idx = int(r / dr)
                if bin_idx < len(hist):
                    hist[bin_idx] += 2 # ペアは2つ(i->j, j->i)
    
#    print("Total pairs considered:", count)        # ちゃんとペア数を数えてるかのチェック用

    # 正規化
    rho = N / (L**2)                                # 全体の密度
    r = (bins[1:] + bins[:-1]) / 2
    shell_volume = 2 * np.pi * r * dr               # 2Dのシェル面積    
    g_r = hist / (N * rho * shell_volume)
    """
    g(r)は、ある粒子から距離rにある粒子の数密度を、全体の平均密度で割ったもの
    2Dでは、シェルの面積が2πr drとなるため、これを用いて正規化を行う
    この正規化はFJCではまだうまくできていない（将来考える）
    """

    return r, g_r

# N個のセグメントからなる自由連結鎖をボックスサイズLの箱の中に配置

try:
    N = int(input('Number of particles (default=1000): '))
except ValueError:
    N = 1000       # 粒子数

# 自由連結鎖の生成
def coordinates_FJC(N):
    x, y = 0, 0
    coordinate_list = [[x,y]]
    for i in range(N):
        theta = 2*np.pi*rd.random()
        x = x + np.cos(theta)
        y = y + np.sin(theta)
        coordinate = [x, y]
        coordinate_list.append(coordinate)
    return coordinate_list

positions = np.array(coordinates_FJC(N))
L = int(np.sqrt(N))       # ボックスサイズ

r, g_r = calculate_rdf(positions, L)

# フィッティング（ここはまだ上手くいってない）
def func(x, a, b):
    return a / x + b

param, cov = curve_fit(func, r, g_r)
a_fit, b_fit = param
fit_curve = func(r, a_fit, b_fit)

# プロット
fig = plt.figure(figsize=(16.0, 8.0))
ax1 = fig.add_subplot(121)
ax1.plot(positions[:,0], positions[:,1], marker='o', markersize=3, mfc='red', mec='none', linestyle='-', c='blue', linewidth=1)
ax1.set_xlim(-L, L)
ax1.set_ylim(-L, L)
ax1.set_title('Distribution of FJC segments')

result_text1 = "Number of segments: {}".format(N)
result_text2 = "Box size: {}".format(L)

fig.text(0.75, 0.8, result_text1)
fig.text(0.75, 0.75, result_text2)

ax2 = fig.add_subplot(122)
ax2.plot(r, g_r)
# ax2.plot(r, fit_curve, linestyle='--', color='orange')
ax2.axhline(1, color='red', linestyle='--')         # 無相関だとg(r)=1
ax2.set_title('Pair Correlation Function of FJC')
ax2.set_xlabel('Distance r')
ax2.set_ylabel('Pair Correlation Function g(r)') 

savefile = "./png/pairCorrelation_FJC_N{0}xL{1}.png".format(N, L)
plt.savefig(savefile)

plt.tight_layout()
plt.show()