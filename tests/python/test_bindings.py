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
