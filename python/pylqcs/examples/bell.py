"""最小示例：Bell 态制备、采样与统计。"""
import pylqcs as lq

qc = lq.QuantumCircuit(2, 2)
qc.h(0).cx(0, 1).measure_all()
print(qc.draw())

result = lq.StatevectorSimulator(seed=42).run(qc, shots=1024)
print("counts:", result.counts())
print("probabilities:", {k: round(v, 3) for k, v in result.probabilities().items()})
