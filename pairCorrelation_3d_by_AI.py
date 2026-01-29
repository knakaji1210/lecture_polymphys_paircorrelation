import numpy as np
import matplotlib.pyplot as plt

def calculate_rdf(positions, L, dr=0.1):
    """
    3D系の対相関関数 g(r) を計算する簡単な関数
    positions: (N, 3)の配列 (原子の座標)
    L: ボックスサイズ (立方体)
    dr: ビンの幅
    """
    N = len(positions)
    max_r = L / 2
    bins = np.arange(0, max_r, dr)
    hist = np.zeros(len(bins)-1)
    
    # 周期境界条件を考慮した距離の計算（簡略化のため全ペア計算）
    for i in range(N):
        for j in range(i + 1, N):
            # 最小イメージ法（周期的境界条件）
            diff = positions[i] - positions[j]
            diff = diff - L * np.round(diff / L)
            r = np.linalg.norm(diff)
            
            if r < max_r:
                bin_idx = int(r / dr)
                if bin_idx < len(hist):
                    hist[bin_idx] += 2 # ペアは2つ(i->j, j->i)
    
    # 正規化
    rho = N / (L**3)
    r = (bins[1:] + bins[:-1]) / 2
    shell_volume = 4 * np.pi * r**2 * dr
    g_r = hist / (N * rho * shell_volume)
    
    return r, g_r

# --- 実行例 ---
# 1000個のランダムな粒子 (L=10)
N = 1000
L = 10.0
np.random.seed(42)
positions = np.random.rand(N, 3) * L

r, g_r = calculate_rdf(positions, L)

# プロット
plt.plot(r, g_r)
plt.axhline(1, color='red', linestyle='--') # 理想気体はg(r)=1
plt.xlabel('Distance r')
plt.ylabel('Pair Correlation Function g(r)')
plt.show()
