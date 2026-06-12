"""matplotlib 可视化：测量直方图、态振幅图、电路图、Shor 分布图。

matplotlib 为惰性导入：未安装时其余功能不受影响。
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from matplotlib.figure import Figure


def _require_matplotlib():
    try:
        import matplotlib.pyplot as plt
        return plt
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "visualization requires matplotlib: pip install matplotlib"
        ) from exc


def plot_histogram(counts: dict, title: str = "Measurement counts",
                   ax=None, color: str = "#5B8FF9") -> "Figure":
    """测量结果柱状图（类 Qiskit）。counts: {'0101': 312, ...}"""
    plt = _require_matplotlib()
    if ax is None:
        fig, ax = plt.subplots(figsize=(max(6, 0.6 * len(counts)), 4))
    else:
        fig = ax.figure
    keys = sorted(counts)
    values = [counts[k] for k in keys]
    total = sum(values)
    bars = ax.bar(keys, values, color=color)
    for bar, v in zip(bars, values):
        ax.annotate(f"{v / total:.3f}", (bar.get_x() + bar.get_width() / 2, v),
                    ha="center", va="bottom", fontsize=8)
    ax.set_xlabel("outcome")
    ax.set_ylabel("counts")
    ax.set_title(title)
    if len(keys) > 8:
        ax.tick_params(axis="x", rotation=70)
    fig.tight_layout()
    return fig


def plot_state(statevector, title: str = "Statevector",
               max_states: int = 64) -> "Figure":
    """态振幅图：柱高 = |amp|²，柱色 = 相位（HSV 色环）。"""
    plt = _require_matplotlib()
    amps = np.asarray(statevector.to_numpy())
    n = int(math.log2(len(amps)))
    probs = np.abs(amps) ** 2
    phases = np.angle(amps)

    if len(amps) > max_states:  # 只画概率最大的 max_states 个基态
        idx = np.sort(np.argsort(probs)[-max_states:])
    else:
        idx = np.arange(len(amps))

    fig, ax = plt.subplots(figsize=(max(6, 0.35 * len(idx)), 4))
    cmap = plt.get_cmap("hsv")
    colors = cmap((phases[idx] + np.pi) / (2 * np.pi))
    ax.bar(range(len(idx)), probs[idx], color=colors)
    ax.set_xticks(range(len(idx)))
    ax.set_xticklabels([f"|{i:0{n}b}>" for i in idx], rotation=70, fontsize=7)
    ax.set_ylabel("probability")
    ax.set_title(f"{title} (color = phase)")

    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(vmin=-np.pi, vmax=np.pi))
    fig.colorbar(sm, ax=ax, label="phase [rad]")
    fig.tight_layout()
    return fig


def plot_circuit(circuit, title: str = "") -> "Figure":
    """电路图：渲染 circuit.draw() 的 ASCII 版图（等宽字体）。"""
    plt = _require_matplotlib()
    art = circuit.draw()
    lines = art.rstrip("\n").split("\n")
    width = max(len(line) for line in lines)
    fig, ax = plt.subplots(
        figsize=(max(4, 0.11 * width), max(2, 0.3 * len(lines))))
    ax.axis("off")
    ax.text(0, 1, art, family="monospace", fontsize=9, va="top")
    if title:
        ax.set_title(title)
    fig.tight_layout()
    return fig


def plot_shor(shor_result, title: str | None = None) -> "Figure":
    """Shor 专用：计数寄存器测量分布 + 理论峰位置 s·2^t/r 标注。"""
    plt = _require_matplotlib()
    if not shor_result.counts:
        raise ValueError(
            "no quantum measurement data (N was factored classically)")

    counts = shor_result.counts
    t = len(next(iter(counts)))          # 计数比特数
    dim = 2 ** t
    xs = [int(k, 2) for k in counts]
    ys = [counts[k] for k in counts]

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(xs, ys, width=max(1.0, dim / 256), color="#5B8FF9",
           label="measured")

    r = shor_result.period
    if r:
        peaks = [s * dim / r for s in range(r)]
        for i, peak in enumerate(peaks):
            ax.axvline(peak, color="#E8684A", linestyle="--", alpha=0.7,
                       label="theory peaks s·2^t/r" if i == 0 else None)
    ax.set_xlabel(f"counting register value (t = {t} bits)")
    ax.set_ylabel("counts")
    ax.set_title(title or
                 f"Shor: N={shor_result.N}={shor_result.factor1}x"
                 f"{shor_result.factor2}, a={shor_result.a}, r={r}")
    ax.legend()
    fig.tight_layout()
    return fig
