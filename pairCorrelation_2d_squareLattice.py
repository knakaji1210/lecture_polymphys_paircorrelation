# Pair Correlation Function (2D)
# 2Dの粒子配置から対相関関数g(r)を計算し、プロットするコード
# これはRadial Distribution Function (RDF)とも呼ばれる

# オリジナルはAIによって生成された3D版を参考にしている。
# 元のコード: code/PairCorrelation/pairCorrelation_3d_by_AI.py

# 粒子が格子点上に規則正しく配置されている場合

import numpy as np
import matplotlib.pyplot as plt

def calculate_rdf(positions, L, dr=0.1):
    """
    2D系の対相関関数 g(r) を計算する簡単な関数
    positions: ((N+1)^2, 2)の配列 (粒子の座標)
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
    """

    return r, g_r

# ボックスサイズLの箱の中の全ての格子点に粒子を配置

try:
    L = float(input('Box size (default=10.0): '))
except ValueError:
    L = 10.0       # ボックスサイズ

# 粒子が格子点上に規則正しく配置されている場合
def regular_lattice_points_in_box(N, M):
    positions = []
    for x in range(N + 1):
        for y in range(M + 1):
            positions.append((x - L/2, y - L/2))  # 中心を原点に揃える
    
    return positions

positions = np.array(regular_lattice_points_in_box(int(L), int(L)))

r, g_r = calculate_rdf(positions, L)

# プロット
fig = plt.figure(figsize=(16.0, 8.0))
ax1 = fig.add_subplot(121)
ax1.scatter(positions[:,0], positions[:,1], s=1)
ax1.set_title('Particle Distribution')

result_text1 = "Number of particles: {}".format(len(positions))
result_text2 = "Box size: {}".format(L)

fig.text(0.75, 0.8, result_text1)
fig.text(0.75, 0.75, result_text2)

ax2 = fig.add_subplot(122)
ax2.plot(r, g_r)
# ax2.axhline(1, color='red', linestyle='--')         # 無相関だとg(r)=1
ax2.set_title('Pair Correlation Function of Particles')
ax2.set_xlabel('Distance r')
ax2.set_ylabel('Pair Correlation Function g(r)') 

savefile = "./png/pairCorrelation_squareLattice_N{0}xL{1}.png".format(len(positions), L)
plt.savefig(savefile)

plt.tight_layout()
plt.show()