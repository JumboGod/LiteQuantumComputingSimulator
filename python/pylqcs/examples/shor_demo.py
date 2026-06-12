"""Shor 算法端到端演示 + 全套可视化。

用法::

    PYTHONPATH=python python3 python/pylqcs/examples/shor_demo.py [输出目录]

交互环境（Jupyter 等）直接显示图像；无显示环境保存 PNG。
"""
import sys

import matplotlib

if not sys.flags.interactive:
    matplotlib.use("Agg")

import pylqcs as lq
from pylqcs.visualization import plot_histogram, plot_shor, plot_state

out_dir = sys.argv[1] if len(sys.argv) > 1 else "."

# —— 1. Shor 分解 N=15（扫种子保证走量子求周期路径，便于画分布图）——
res = None
for seed in range(1, 50):
    r = lq.algorithms.shor(15, seed=seed)
    if r.success and r.period != 0:
        res = r
        break
assert res is not None
print(f"Shor: {res.N} = {res.factor1} x {res.factor2}  "
      f"(a={res.a}, r={res.period}, attempts={res.attempts})")
fig = plot_shor(res)
fig.savefig(f"{out_dir}/shor_15.png", dpi=120)
print(f"  -> {out_dir}/shor_15.png")

# —— 2. Bell 态：测量直方图 + 态振幅图 ——
qc = lq.QuantumCircuit(2, 2)
qc.h(0).cx(0, 1).measure_all()
print("\nBell circuit:")
print(qc.draw())

result = lq.StatevectorSimulator(seed=7).run(qc, shots=2048)
print(f"counts: {result.counts()}")
plot_histogram(result.counts(), title="Bell state, 2048 shots") \
    .savefig(f"{out_dir}/bell_histogram.png", dpi=120)
print(f"  -> {out_dir}/bell_histogram.png")

sv = lq.StatevectorSimulator().run_statevector(
    lq.QuantumCircuit(3).h(0).cx(0, 1).cx(1, 2).t(2))
plot_state(sv, title="GHZ + T") \
    .savefig(f"{out_dir}/ghz_state.png", dpi=120)
print(f"  -> {out_dir}/ghz_state.png")

if sys.flags.interactive:
    import matplotlib.pyplot as plt
    plt.show()
