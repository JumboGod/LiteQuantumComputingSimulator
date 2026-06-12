"""pylqcs — LiteQuantumComputingSimulator 的 Python 接口。

C++ 负责全部数值计算，Python 侧提供交互与 matplotlib 可视化。

快速上手::

    import pylqcs as lq

    qc = lq.QuantumCircuit(2, 2)
    qc.h(0).cx(0, 1).measure_all()
    result = lq.StatevectorSimulator().run(qc, shots=1024)
    print(result.counts())          # {'00': ~512, '11': ~512}

    res = lq.algorithms.shor(15, seed=42)
    print(res.factor1, res.factor2)  # 3 5
"""

from ._pylqcs import (  # noqa: F401
    DensityMatrix,
    DensityMatrixSimulator,
    Gate,
    GateType,
    KrausChannel,
    NoiseModel,
    QuantumCircuit,
    Result,
    Statevector,
    StatevectorSimulator,
    algorithms,
    io,
)

__version__ = "0.1.0"

__all__ = [
    "DensityMatrix",
    "DensityMatrixSimulator",
    "Gate",
    "GateType",
    "KrausChannel",
    "NoiseModel",
    "QuantumCircuit",
    "Result",
    "Statevector",
    "StatevectorSimulator",
    "algorithms",
    "io",
]
