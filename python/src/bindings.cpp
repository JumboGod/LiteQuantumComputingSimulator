// pylqcs Python 绑定：薄封装，所有计算留在 C++ 核心
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <optional>
#include <vector>

#include <lqcs/lqcs.hpp>

namespace py = pybind11;
using namespace py::literals;
using namespace lqcs;

namespace {

using NpComplexArray =
    py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast>;

std::vector<complex_t> to_matrix(const NpComplexArray& arr) {
    const auto buf = arr.request();
    const auto* p = static_cast<const complex_t*>(buf.ptr);
    return {p, p + buf.size};
}

}  // namespace

PYBIND11_MODULE(_pylqcs, m) {
    m.doc() = "LiteQuantumComputingSimulator Python bindings";

    // —— GateType / Gate ——
    py::enum_<GateType>(m, "GateType")
        .value("I", GateType::I).value("X", GateType::X)
        .value("Y", GateType::Y).value("Z", GateType::Z)
        .value("H", GateType::H).value("S", GateType::S)
        .value("Sdg", GateType::Sdg).value("T", GateType::T)
        .value("Tdg", GateType::Tdg).value("SX", GateType::SX)
        .value("RX", GateType::RX).value("RY", GateType::RY)
        .value("RZ", GateType::RZ).value("P", GateType::P)
        .value("U", GateType::U).value("SWAP", GateType::SWAP)
        .value("iSWAP", GateType::iSWAP).value("RXX", GateType::RXX)
        .value("RYY", GateType::RYY).value("RZZ", GateType::RZZ)
        .value("Unitary", GateType::Unitary)
        .value("Permutation", GateType::Permutation)
        .value("Measure", GateType::Measure).value("Reset", GateType::Reset)
        .value("Barrier", GateType::Barrier);

    py::class_<Gate>(m, "Gate")
        .def(py::init<GateType, std::vector<double>, std::size_t>(),
             "type"_a, "params"_a = std::vector<double>{}, "n_controls"_a = 0)
        .def_readonly("type", &Gate::type)
        .def_readonly("params", &Gate::params)
        .def_readonly("n_controls", &Gate::n_controls)
        .def("num_qubits", &Gate::num_qubits)
        .def("name", &Gate::name)
        .def("inverse", &Gate::inverse)
        .def("matrix",
             [](const Gate& g) {
                 const auto v = g.matrix();
                 const py::ssize_t dim = py::ssize_t{1} << g.num_qubits();
                 py::array_t<std::complex<double>> arr({dim, dim});
                 std::copy(v.begin(), v.end(),
                           static_cast<std::complex<double>*>(arr.request().ptr));
                 return arr;
             })
        .def("__repr__", [](const Gate& g) { return "<Gate " + g.name() + ">"; });

    // —— QuantumCircuit（流式 API：方法返回自身引用，支持链式调用）——
    auto rp = py::return_value_policy::reference_internal;
    py::class_<QuantumCircuit>(m, "QuantumCircuit")
        .def(py::init<std::size_t, std::size_t>(), "num_qubits"_a,
             "num_clbits"_a = 0)
        .def("i", &QuantumCircuit::i, rp).def("x", &QuantumCircuit::x, rp)
        .def("y", &QuantumCircuit::y, rp).def("z", &QuantumCircuit::z, rp)
        .def("h", &QuantumCircuit::h, rp).def("s", &QuantumCircuit::s, rp)
        .def("sdg", &QuantumCircuit::sdg, rp).def("t", &QuantumCircuit::t, rp)
        .def("tdg", &QuantumCircuit::tdg, rp).def("sx", &QuantumCircuit::sx, rp)
        .def("rx", &QuantumCircuit::rx, rp).def("ry", &QuantumCircuit::ry, rp)
        .def("rz", &QuantumCircuit::rz, rp).def("p", &QuantumCircuit::p, rp)
        .def("u", &QuantumCircuit::u, rp)
        .def("cx", &QuantumCircuit::cx, rp).def("cy", &QuantumCircuit::cy, rp)
        .def("cz", &QuantumCircuit::cz, rp).def("ch", &QuantumCircuit::ch, rp)
        .def("cp", &QuantumCircuit::cp, rp).def("crx", &QuantumCircuit::crx, rp)
        .def("cry", &QuantumCircuit::cry, rp).def("crz", &QuantumCircuit::crz, rp)
        .def("ccx", &QuantumCircuit::ccx, rp)
        .def("cswap", &QuantumCircuit::cswap, rp)
        .def("mcx",
             [](QuantumCircuit& qc, std::vector<qubit_t> controls,
                qubit_t target) -> QuantumCircuit& {
                 return qc.mcx(controls, target);
             },
             rp, "controls"_a, "target"_a)
        .def("mcp",
             [](QuantumCircuit& qc, double lam, std::vector<qubit_t> controls,
                qubit_t target) -> QuantumCircuit& {
                 return qc.mcp(lam, controls, target);
             },
             rp, "lam"_a, "controls"_a, "target"_a)
        .def("swap", &QuantumCircuit::swap, rp)
        .def("iswap", &QuantumCircuit::iswap, rp)
        .def("rxx", &QuantumCircuit::rxx, rp)
        .def("ryy", &QuantumCircuit::ryy, rp)
        .def("rzz", &QuantumCircuit::rzz, rp)
        .def("unitary",
             [](QuantumCircuit& qc, const NpComplexArray& mat,
                std::vector<qubit_t> qubits) -> QuantumCircuit& {
                 return qc.unitary(to_matrix(mat), qubits);
             },
             rp, "matrix"_a, "qubits"_a)
        .def("permutation",
             [](QuantumCircuit& qc, std::vector<std::size_t> table,
                std::vector<qubit_t> qubits) -> QuantumCircuit& {
                 return qc.permutation(table, qubits);
             },
             rp, "table"_a, "qubits"_a)
        .def("measure", &QuantumCircuit::measure, rp)
        .def("measure_all", &QuantumCircuit::measure_all, rp)
        .def("reset", &QuantumCircuit::reset, rp)
        .def("barrier", &QuantumCircuit::barrier, rp)
        .def("compose",
             [](QuantumCircuit& qc, const QuantumCircuit& other,
                std::vector<qubit_t> qubit_map) -> QuantumCircuit& {
                 return qc.compose(other, qubit_map);
             },
             rp, "other"_a, "qubit_map"_a = std::vector<qubit_t>{})
        .def("inverse", &QuantumCircuit::inverse)
        .def("power", &QuantumCircuit::power)
        .def("num_qubits", &QuantumCircuit::num_qubits)
        .def("num_clbits", &QuantumCircuit::num_clbits)
        .def("num_instructions", &QuantumCircuit::num_instructions)
        .def("depth", &QuantumCircuit::depth)
        .def("draw", &QuantumCircuit::draw)
        .def("__str__", &QuantumCircuit::draw)
        .def("__repr__", [](const QuantumCircuit& qc) {
            return "<QuantumCircuit qubits=" + std::to_string(qc.num_qubits()) +
                   " clbits=" + std::to_string(qc.num_clbits()) +
                   " instructions=" + std::to_string(qc.num_instructions()) + ">";
        });

    // —— Statevector ——
    py::class_<Statevector>(m, "Statevector")
        .def(py::init<std::size_t>(), "num_qubits"_a)
        .def("num_qubits", &Statevector::num_qubits)
        .def("norm", &Statevector::norm)
        .def("__len__", &Statevector::size)
        .def("__getitem__",
             [](const Statevector& sv, std::size_t i) {
                 if (i >= sv.size()) throw py::index_error();
                 return static_cast<std::complex<double>>(sv[i]);
             })
        .def("probabilities",
             [](const Statevector& sv) {
                 const auto p = sv.probabilities();
                 py::array_t<double> arr(static_cast<py::ssize_t>(p.size()));
                 std::copy(p.begin(), p.end(),
                           static_cast<double*>(arr.request().ptr));
                 return arr;
             })
        .def("expectation_value", &Statevector::expectation_value, "pauli"_a)
        // 零拷贝 numpy 视图：base 持有 Statevector 保证生命周期
        .def("to_numpy", [](py::object self) {
            auto& sv = self.cast<Statevector&>();
            return py::array_t<std::complex<double>>(
                {static_cast<py::ssize_t>(sv.size())},
                {static_cast<py::ssize_t>(sizeof(complex_t))},
                reinterpret_cast<std::complex<double>*>(sv.data()), self);
        });

    // —— Result ——
    py::class_<Result>(m, "Result")
        .def("counts", &Result::counts)
        .def("probabilities", &Result::probabilities)
        .def("shots", &Result::shots)
        .def("__repr__", [](const Result& r) {
            std::string s = "<Result shots=" + std::to_string(r.shots()) + " {";
            bool first = true;
            for (const auto& [k, v] : r.counts()) {
                if (!first) s += ", ";
                s += "'" + k + "': " + std::to_string(v);
                first = false;
            }
            return s + "}>";
        });

    // —— StatevectorSimulator ——
    py::class_<StatevectorSimulator>(m, "StatevectorSimulator")
        .def(py::init([](std::optional<std::uint64_t> seed) {
                 return StatevectorSimulator({.seed = seed});
             }),
             "seed"_a = py::none())
        .def("run", &StatevectorSimulator::run, "circuit"_a, "shots"_a = 1024,
             py::call_guard<py::gil_scoped_release>())
        .def("run_statevector", &StatevectorSimulator::run_statevector,
             "circuit"_a, py::call_guard<py::gil_scoped_release>())
        .def("name", &StatevectorSimulator::name);

    // —— 算法库 ——
    auto alg = m.def_submodule("algorithms", "quantum algorithms library");
    using namespace lqcs::algorithms;

    alg.def("qft", &qft, "n"_a, "do_swaps"_a = true);
    alg.def("qft_inverse", &qft_inverse, "n"_a, "do_swaps"_a = true);
    alg.def("xor_oracle", &xor_oracle, "n_in"_a, "n_out"_a, "f"_a);
    alg.def("phase_oracle",
            [](std::size_t n, std::vector<std::uint64_t> marked) {
                return phase_oracle(n, marked);
            },
            "n"_a, "marked"_a);
    alg.def("modular_exponentiation", &modular_exponentiation, "a"_a, "N"_a,
            "n_counting"_a, "n_work"_a);

    py::class_<DeutschJozsaResult>(alg, "DeutschJozsaResult")
        .def_readonly("is_constant", &DeutschJozsaResult::is_constant)
        .def_readonly("p_all_zero", &DeutschJozsaResult::p_all_zero)
        .def_readonly("circuit", &DeutschJozsaResult::circuit);
    alg.def("deutsch_jozsa", &deutsch_jozsa, "n"_a, "f"_a);

    py::class_<BernsteinVaziraniResult>(alg, "BernsteinVaziraniResult")
        .def_readonly("secret", &BernsteinVaziraniResult::secret)
        .def_readonly("circuit", &BernsteinVaziraniResult::circuit);
    alg.def("bernstein_vazirani", &bernstein_vazirani, "n"_a, "f"_a);

    py::class_<SimonResult>(alg, "SimonResult")
        .def_readonly("secret", &SimonResult::secret)
        .def_readonly("quantum_samples", &SimonResult::quantum_samples)
        .def_readonly("circuit", &SimonResult::circuit);
    alg.def("simon", &simon, "n"_a, "f"_a, "seed"_a = py::none());

    py::class_<QPEResult>(alg, "QPEResult")
        .def_readonly("phase", &QPEResult::phase)
        .def_readonly("measured", &QPEResult::measured)
        .def_readonly("n_counting", &QPEResult::n_counting)
        .def_readonly("circuit", &QPEResult::circuit);
    alg.def("phase_estimation", &phase_estimation, "u"_a, "eigenstate_prep"_a,
            "n_counting"_a);

    py::class_<GroverResult>(alg, "GroverResult")
        .def_readonly("best_state", &GroverResult::best_state)
        .def_readonly("success_probability", &GroverResult::success_probability)
        .def_readonly("iterations", &GroverResult::iterations)
        .def_readonly("circuit", &GroverResult::circuit);
    alg.def("grover",
            [](std::size_t n, std::vector<std::uint64_t> marked,
               std::optional<std::size_t> iterations) {
                return grover(n, marked, iterations);
            },
            "n"_a, "marked"_a, "iterations"_a = py::none());

    py::class_<ShorResult>(alg, "ShorResult")
        .def_readonly("success", &ShorResult::success)
        .def_readonly("N", &ShorResult::N)
        .def_readonly("factor1", &ShorResult::factor1)
        .def_readonly("factor2", &ShorResult::factor2)
        .def_readonly("a", &ShorResult::a)
        .def_readonly("period", &ShorResult::period)
        .def_readonly("attempts", &ShorResult::attempts)
        .def_readonly("counts", &ShorResult::counts)
        .def_readonly("circuit", &ShorResult::circuit)
        .def("__repr__", [](const ShorResult& r) {
            if (!r.success) {
                return "<ShorResult N=" + std::to_string(r.N) + " failed>";
            }
            return "<ShorResult " + std::to_string(r.N) + " = " +
                   std::to_string(r.factor1) + " x " + std::to_string(r.factor2) +
                   " (a=" + std::to_string(r.a) +
                   ", r=" + std::to_string(r.period) + ")>";
        });
    alg.def(
        "shor",
        [](std::uint64_t N, std::size_t shots, std::size_t max_attempts,
           std::optional<std::uint64_t> seed) {
            py::gil_scoped_release release;
            return shor(N, {shots, max_attempts, seed});
        },
        "N"_a, "shots"_a = 64, "max_attempts"_a = 10, "seed"_a = py::none());

    // —— VQE ——
    py::class_<VQEResult>(alg, "VQEResult")
        .def_readonly("energy", &VQEResult::energy)
        .def_readonly("parameters", &VQEResult::parameters)
        .def_readonly("iterations", &VQEResult::iterations)
        .def_readonly("history", &VQEResult::history);
    alg.def(
        "expectation",
        [](const Statevector& sv,
           const std::vector<std::pair<double, std::string>>& terms) {
            Hamiltonian h;
            for (const auto& [c, p] : terms) h.push_back({c, p});
            return expectation(sv, h);
        },
        "statevector"_a, "hamiltonian"_a,
        "<psi|H|psi>, hamiltonian = [(coeff, pauli_str), ...]");
    alg.def(
        "vqe",
        [](const std::vector<std::pair<double, std::string>>& terms,
           const std::function<QuantumCircuit(std::vector<double>)>& ansatz,
           std::size_t n_params, std::size_t max_iterations, double tol,
           std::vector<double> initial_parameters) {
            Hamiltonian h;
            for (const auto& [c, p] : terms) h.push_back({c, p});
            const Ansatz wrapped = [&ansatz](std::span<const double> theta) {
                return ansatz({theta.begin(), theta.end()});
            };
            return vqe(h, wrapped, n_params,
                       {max_iterations, tol, std::move(initial_parameters)});
        },
        "hamiltonian"_a, "ansatz"_a, "n_params"_a, "max_iterations"_a = 100,
        "tol"_a = 1e-9, "initial_parameters"_a = std::vector<double>{});

    auto nt = alg.def_submodule("number_theory", "classical number theory");
    nt.def("gcd", &number_theory::gcd);
    nt.def("pow_mod", &number_theory::pow_mod);
    nt.def("mul_mod", &number_theory::mul_mod);
    nt.def("is_prime", &number_theory::is_prime);
    nt.def("perfect_power_base", &number_theory::perfect_power_base);
    nt.def("continued_fraction_denominators",
           &number_theory::continued_fraction_denominators, "num"_a, "den"_a,
           "max_denominator"_a);

    // —— OpenQASM 2.0 ——
    auto io_mod = m.def_submodule("io", "OpenQASM import/export");
    io_mod.def("to_qasm", &io::to_qasm, "circuit"_a);
    io_mod.def("from_qasm", &io::from_qasm, "source"_a);

    // —— 噪声与密度矩阵后端 ——
    py::class_<KrausChannel>(m, "KrausChannel")
        .def_static("depolarizing", &KrausChannel::depolarizing, "p"_a)
        .def_static("amplitude_damping", &KrausChannel::amplitude_damping,
                    "gamma"_a)
        .def_static("phase_damping", &KrausChannel::phase_damping, "gamma"_a)
        .def_static("bit_flip", &KrausChannel::bit_flip, "p"_a)
        .def_static("phase_flip", &KrausChannel::phase_flip, "p"_a);

    py::class_<NoiseModel>(m, "NoiseModel")
        .def(py::init<>())
        .def("add_all_qubit_channel", &NoiseModel::add_all_qubit_channel,
             py::return_value_policy::reference_internal, "channel"_a)
        .def("empty", &NoiseModel::empty);

    py::class_<DensityMatrix>(m, "DensityMatrix")
        .def(py::init<std::size_t>(), "num_qubits"_a)
        .def("num_qubits", &DensityMatrix::num_qubits)
        .def("dim", &DensityMatrix::dim)
        .def("element", &DensityMatrix::element, "row"_a, "col"_a)
        .def("apply_channel", &DensityMatrix::apply_channel, "channel"_a,
             "qubit"_a)
        .def("trace", &DensityMatrix::trace)
        .def("purity", &DensityMatrix::purity)
        .def("fidelity", &DensityMatrix::fidelity, "psi"_a)
        .def("probabilities", [](const DensityMatrix& rho) {
            const auto p = rho.probabilities();
            py::array_t<double> arr(static_cast<py::ssize_t>(p.size()));
            std::copy(p.begin(), p.end(),
                      static_cast<double*>(arr.request().ptr));
            return arr;
        });

    py::class_<DensityMatrixSimulator>(m, "DensityMatrixSimulator")
        .def(py::init([](std::optional<std::uint64_t> seed, NoiseModel noise) {
                 return DensityMatrixSimulator({seed, std::move(noise)});
             }),
             "seed"_a = py::none(), "noise"_a = NoiseModel{})
        .def("run", &DensityMatrixSimulator::run, "circuit"_a, "shots"_a = 1024,
             py::call_guard<py::gil_scoped_release>())
        .def("run_density_matrix", &DensityMatrixSimulator::run_density_matrix,
             "circuit"_a, py::call_guard<py::gil_scoped_release>())
        .def("name", &DensityMatrixSimulator::name);
}
