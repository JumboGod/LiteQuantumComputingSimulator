"""pylqcs 绑定层测试：链式 API、零拷贝、算法调用"""
import math

import numpy as np
import pytest

import pylqcs as lq


def test_fluent_chain_and_counts():
    qc = lq.QuantumCircuit(2, 2)
    qc.h(0).cx(0, 1).measure_all()
    assert qc.num_qubits() == 2
    assert qc.depth() == 3

    sim = lq.StatevectorSimulator(seed=42)
    result = sim.run(qc, shots=1024)
    counts = result.counts()
    assert isinstance(counts, dict)
    assert set(counts) <= {"00", "11"}
    assert sum(counts.values()) == 1024


def test_statevector_numpy_zero_copy():
    qc = lq.QuantumCircuit(2)
    qc.h(0).cx(0, 1)
    sv = lq.StatevectorSimulator().run_statevector(qc)
    arr = sv.to_numpy()
    assert arr.dtype == np.complex128
    assert arr.shape == (4,)
    s = 1 / math.sqrt(2)
    np.testing.assert_allclose(arr, [s, 0, 0, s], atol=1e-12)
    # 零拷贝验证：通过 numpy 写入可从 C++ 侧读出
    arr[1] = 0.5 + 0.25j
    assert sv[1] == pytest.approx(0.5 + 0.25j)
    # 生命周期：numpy 数组持有 Statevector
    assert arr.base is sv


def test_statevector_rejects_bad_index():
    sv = lq.Statevector(2)
    with pytest.raises(IndexError):
        _ = sv[4]


def test_build_time_validation_raises():
    qc = lq.QuantumCircuit(2)
    with pytest.raises(Exception):
        qc.h(5)
    with pytest.raises(Exception):
        qc.unitary(np.array([[1, 1], [0, 1]], dtype=complex), [0])


def test_unitary_from_numpy():
    h = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)
    qc = lq.QuantumCircuit(1)
    qc.unitary(h, [0])
    sv = lq.StatevectorSimulator().run_statevector(qc)
    np.testing.assert_allclose(sv.to_numpy(),
                               [1 / math.sqrt(2), 1 / math.sqrt(2)], atol=1e-12)


def test_circuit_draw_and_repr():
    qc = lq.QuantumCircuit(2, 2)
    qc.h(0).cx(0, 1).measure_all()
    art = qc.draw()
    assert "q_0:" in art and "[h]" in art
    assert "instructions=4" in repr(qc)


def test_shor_factors_15():
    res = lq.algorithms.shor(15, seed=42)
    assert res.success
    assert {res.factor1, res.factor2} == {3, 5}
    assert res.factor1 * res.factor2 == 15


def test_shor_result_counts_dict():
    # 扫种子找一个量子路径成功的（period != 0），验证 counts/circuit 暴露
    for seed in range(1, 30):
        res = lq.algorithms.shor(15, seed=seed)
        assert res.success
        if res.period != 0:
            assert isinstance(res.counts, dict) and res.counts
            assert res.circuit is not None
            assert res.circuit.num_qubits() == 12
            return
    pytest.fail("no seed exercised the quantum path")


def test_deutsch_jozsa_with_python_lambda():
    const = lq.algorithms.deutsch_jozsa(3, lambda x: True)
    assert const.is_constant
    balanced = lq.algorithms.deutsch_jozsa(3, lambda x: bool(x & 1))
    assert not balanced.is_constant


def test_bernstein_vazirani():
    secret = 0b1011
    res = lq.algorithms.bernstein_vazirani(
        4, lambda x: bin(x & secret).count("1") % 2 == 1)
    assert res.secret == secret


def test_simon():
    s = 0b110
    res = lq.algorithms.simon(3, lambda x: min(x, x ^ s), seed=42)
    assert res.secret == s


def test_qpe_t_gate():
    prep = lq.QuantumCircuit(1)
    prep.x(0)
    res = lq.algorithms.phase_estimation(lq.Gate(lq.GateType.T), prep, 3)
    assert res.phase == pytest.approx(0.125)


def test_grover():
    res = lq.algorithms.grover(4, [9])
    assert res.best_state == 9
    assert res.success_probability > 0.9


def test_number_theory():
    nt = lq.algorithms.number_theory
    assert nt.gcd(48, 18) == 6
    assert nt.pow_mod(7, 4, 15) == 1
    assert nt.is_prime(17)
    assert 4 in nt.continued_fraction_denominators(192, 256, 15)


def test_expectation_value():
    qc = lq.QuantumCircuit(2)
    qc.h(0).cx(0, 1)
    sv = lq.StatevectorSimulator().run_statevector(qc)
    assert sv.expectation_value("ZZ") == pytest.approx(1.0)
    assert sv.expectation_value("XX") == pytest.approx(1.0)
    assert sv.expectation_value("YY") == pytest.approx(-1.0)


H2 = [
    (-1.052373245772859, "II"),
    (+0.39793742484318045, "IZ"),
    (-0.39793742484318045, "ZI"),
    (-0.01128010425623538, "ZZ"),
    (+0.18093119978423156, "XX"),
]


def test_vqe_h2_with_python_ansatz():
    def ansatz(theta):
        qc = lq.QuantumCircuit(2)
        qc.x(0).ry(theta[0], 1).cx(1, 0)
        return qc

    res = lq.algorithms.vqe(H2, ansatz, 1, max_iterations=50)
    assert res.energy == pytest.approx(-1.8572750302023824, abs=1e-9)
    assert len(res.parameters) == 1


def test_qasm_round_trip():
    qc = lq.QuantumCircuit(2, 2)
    qc.h(0).cx(0, 1).rz(0.5, 1).measure_all()
    qasm = lq.io.to_qasm(qc)
    assert "OPENQASM 2.0;" in qasm
    parsed = lq.io.from_qasm(qasm)
    assert parsed.num_qubits() == 2
    a = lq.StatevectorSimulator(seed=3).run(qc, shots=100).counts()
    b = lq.StatevectorSimulator(seed=3).run(parsed, shots=100).counts()
    assert a == b


def test_noisy_density_matrix():
    bell = lq.QuantumCircuit(2)
    bell.h(0).cx(0, 1)
    ideal = lq.StatevectorSimulator().run_statevector(bell)

    rho = lq.DensityMatrixSimulator().run_density_matrix(bell)
    assert rho.fidelity(ideal) == pytest.approx(1.0)

    p = 0.1
    rho.apply_channel(lq.KrausChannel.depolarizing(p), 0)
    rho.apply_channel(lq.KrausChannel.depolarizing(p), 1)
    expected = (1 - p) ** 2 + p * (2 - p) / 4
    assert rho.fidelity(ideal) == pytest.approx(expected, abs=1e-12)

    noise = lq.NoiseModel()
    noise.add_all_qubit_channel(lq.KrausChannel.depolarizing(0.1))
    bell.measure_all()
    counts = lq.DensityMatrixSimulator(seed=7, noise=noise).run(
        bell, shots=2048).counts()
    assert set(counts) == {"00", "01", "10", "11"}


def test_trajectory_noise_on_statevector_backend():
    noise = lq.NoiseModel()
    noise.add_all_qubit_channel(lq.KrausChannel.depolarizing(0.1))

    bell = lq.QuantumCircuit(2, 2)
    bell.h(0).cx(0, 1).measure_all()
    shots = 20000
    counts = lq.StatevectorSimulator(seed=7, noise=noise).run(
        bell, shots=shots).counts()
    # 与密度矩阵精确解交叉验证
    exact = lq.DensityMatrixSimulator(noise=noise).run_density_matrix(
        lq.QuantumCircuit(2).h(0).cx(0, 1)).probabilities()
    for i, key in enumerate(["00", "01", "10", "11"]):
        assert counts.get(key, 0) / shots == pytest.approx(exact[i], abs=0.02)

    # 16 比特含噪（密度矩阵后端不可能）
    big = lq.QuantumCircuit(16, 16)
    big.h(0)
    for q in range(15):
        big.cx(q, q + 1)
    big.measure_all()
    r = lq.StatevectorSimulator(seed=1, noise=noise).run(big, shots=8)
    assert sum(r.counts().values()) == 8


def test_pauli_rotation_matches_rzz():
    a = lq.QuantumCircuit(2)
    a.h(0).rp(0.7, "ZZ", [0, 1])
    b = lq.QuantumCircuit(2)
    b.h(0).rzz(0.7, 0, 1)
    sa = lq.StatevectorSimulator(fuse_gates=False).run_statevector(a).to_numpy()
    sb = lq.StatevectorSimulator(fuse_gates=False).run_statevector(b).to_numpy()
    np.testing.assert_allclose(sa, sb, atol=1e-12)


def test_parametric_circuit_and_gradient_descent():
    H2 = [
        (-1.052373245772859, "II"),
        (+0.39793742484318045, "IZ"),
        (-0.39793742484318045, "ZI"),
        (-0.01128010425623538, "ZZ"),
        (+0.18093119978423156, "XX"),
    ]
    pc = lq.ParametricCircuit(2)
    pc.x(0)
    p = pc.ry(1)            # 参数索引 0
    assert p == 0
    pc.cx(1, 0)
    assert pc.num_parameters() == 1

    # parameter-shift 梯度 vs 有限差分
    grad = lq.algorithms.parameter_shift_gradient(H2, pc, [0.3])
    sim = lq.StatevectorSimulator()

    def energy(theta):
        return lq.algorithms.expectation(
            sim.run_statevector(pc.bind(theta)), H2)

    eps = 1e-6
    fd = (energy([0.3 + eps]) - energy([0.3 - eps])) / (2 * eps)
    assert grad[0] == pytest.approx(fd, abs=1e-5)

    # 梯度下降收敛到 H2 基态
    res = lq.algorithms.vqe_gradient_descent(
        H2, pc, max_iterations=500, learning_rate=0.3, tol=1e-7)
    assert res.energy == pytest.approx(-1.8572750302023824, abs=1e-5)


def test_visualization_figures():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from pylqcs.visualization import (plot_circuit, plot_histogram,
                                      plot_shor, plot_state)

    qc = lq.QuantumCircuit(2, 2)
    qc.h(0).cx(0, 1).measure_all()
    result = lq.StatevectorSimulator(seed=1).run(qc, shots=256)
    assert plot_histogram(result.counts()) is not None

    sv = lq.StatevectorSimulator().run_statevector(
        lq.QuantumCircuit(2).h(0).cx(0, 1))
    assert plot_state(sv) is not None
    assert plot_circuit(qc) is not None

    for seed in range(1, 30):
        res = lq.algorithms.shor(15, seed=seed)
        if res.period != 0:
            assert plot_shor(res) is not None
            break
