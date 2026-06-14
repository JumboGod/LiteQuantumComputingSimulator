#include "lqcs/io/qasm.hpp"

#include <cctype>
#include <cstdio>
#include <map>
#include <memory>
#include <numbers>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace lqcs::io {

// ============================ 导出 ============================

namespace {

std::string fmt(double v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    return buf;
}

std::string params_str(const std::vector<double>& params) {
    std::string s = "(";
    for (std::size_t i = 0; i < params.size(); ++i) {
        if (i) s += ",";
        s += fmt(params[i]);
    }
    return s + ")";
}

// 基门 → qelib1.inc 门名（含控制位组合）；不可表示时返回空
std::string qasm_gate_name(const Gate& g) {
    const std::size_t nc = g.n_controls;
    switch (g.type) {
        case GateType::I: return nc == 0 ? "id" : "";
        case GateType::X:
            return nc == 0 ? "x" : nc == 1 ? "cx" : nc == 2 ? "ccx" : "";
        case GateType::Y: return nc == 0 ? "y" : nc == 1 ? "cy" : "";
        case GateType::Z: return nc == 0 ? "z" : nc == 1 ? "cz" : "";
        case GateType::H: return nc == 0 ? "h" : nc == 1 ? "ch" : "";
        case GateType::S: return nc == 0 ? "s" : "";
        case GateType::Sdg: return nc == 0 ? "sdg" : "";
        case GateType::T: return nc == 0 ? "t" : "";
        case GateType::Tdg: return nc == 0 ? "tdg" : "";
        case GateType::SX: return nc == 0 ? "sx" : "";
        case GateType::RX: return nc == 0 ? "rx" : nc == 1 ? "crx" : "";
        case GateType::RY: return nc == 0 ? "ry" : nc == 1 ? "cry" : "";
        case GateType::RZ: return nc == 0 ? "rz" : nc == 1 ? "crz" : "";
        case GateType::P: return nc == 0 ? "u1" : nc == 1 ? "cu1" : "";
        case GateType::U: return nc == 0 ? "u3" : "";
        case GateType::SWAP: return nc == 0 ? "swap" : nc == 1 ? "cswap" : "";
        default: return "";
    }
}

}  // namespace

std::string to_qasm(const QuantumCircuit& circuit) {
    std::ostringstream os;
    os << "OPENQASM 2.0;\n"
       << "include \"qelib1.inc\";\n"
       << "qreg q[" << circuit.num_qubits() << "];\n";
    if (circuit.num_clbits() > 0) {
        os << "creg c[" << circuit.num_clbits() << "];\n";
    }

    for (const auto& inst : circuit.instructions()) {
        const Gate& g = inst.gate;
        if (g.type == GateType::Barrier) {
            os << "barrier q;\n";
            continue;
        }
        if (g.type == GateType::Measure) {
            os << "measure q[" << inst.qubits[0] << "] -> c[" << inst.clbits[0]
               << "];\n";
            continue;
        }
        if (g.type == GateType::Reset) {
            os << "reset q[" << inst.qubits[0] << "];\n";
            continue;
        }
        const std::string name = qasm_gate_name(g);
        if (name.empty()) {
            throw std::invalid_argument(
                "to_qasm: gate '" + g.name() +
                "' is not representable in the OpenQASM 2.0 qelib1 gate set");
        }
        os << name;
        if (!g.params.empty()) os << params_str(g.params);
        for (std::size_t i = 0; i < inst.qubits.size(); ++i) {
            os << (i ? "," : " ") << "q[" << inst.qubits[i] << "]";
        }
        os << ";\n";
    }
    return os.str();
}

// ============================ 导入 ============================

namespace {

// 极简词法/语法分析器：覆盖 to_qasm 输出与 Qiskit 常见导出子集
class Parser {
   public:
    explicit Parser(std::string_view src) : src_(src) {}

    QuantumCircuit parse() {
        skip_ws();
        expect_keyword("OPENQASM");
        expect_token("2.0");
        expect_token(";");
        // 头部之后逐条语句处理
        while (skip_ws(), pos_ < src_.size()) {
            parse_statement();
        }
        if (!circuit_) {
            throw std::invalid_argument("from_qasm: missing qreg declaration");
        }
        return std::move(*circuit_);
    }

   private:
    // —— 词法 ——
    void skip_ws() {
        while (pos_ < src_.size()) {
            if (std::isspace(static_cast<unsigned char>(src_[pos_]))) {
                ++pos_;
            } else if (src_.substr(pos_, 2) == "//") {
                while (pos_ < src_.size() && src_[pos_] != '\n') ++pos_;
            } else {
                break;
            }
        }
    }

    char peek() { return pos_ < src_.size() ? src_[pos_] : '\0'; }

    bool consume(char c) {
        skip_ws();
        if (peek() == c) {
            ++pos_;
            return true;
        }
        return false;
    }

    void expect(char c) {
        if (!consume(c)) {
            fail(std::string("expected '") + c + "'");
        }
    }

    std::string identifier() {
        skip_ws();
        std::size_t start = pos_;
        while (pos_ < src_.size() &&
               (std::isalnum(static_cast<unsigned char>(src_[pos_])) ||
                src_[pos_] == '_')) {
            ++pos_;
        }
        if (start == pos_) fail("expected identifier");
        return std::string(src_.substr(start, pos_ - start));
    }

    void expect_keyword(const std::string& kw) {
        if (identifier() != kw) fail("expected '" + kw + "'");
    }

    void expect_token(std::string_view tok) {
        skip_ws();
        if (src_.substr(pos_, tok.size()) != tok) {
            fail("expected '" + std::string(tok) + "'");
        }
        pos_ += tok.size();
    }

    [[noreturn]] void fail(const std::string& msg) {
        throw std::invalid_argument("from_qasm: " + msg + " near position " +
                                    std::to_string(pos_));
    }

    // —— 参数表达式：number | pi | -x | (expr) 与 + - * / ——
    double parse_expr() {
        double v = parse_term();
        while (true) {
            skip_ws();
            if (consume('+'))
                v += parse_term();
            else if (consume('-'))
                v -= parse_term();
            else
                return v;
        }
    }
    double parse_term() {
        double v = parse_factor();
        while (true) {
            skip_ws();
            if (consume('*'))
                v *= parse_factor();
            else if (consume('/'))
                v /= parse_factor();
            else
                return v;
        }
    }
    double parse_factor() {
        skip_ws();
        if (consume('-')) return -parse_factor();
        if (consume('(')) {
            const double v = parse_expr();
            expect(')');
            return v;
        }
        if (std::isalpha(static_cast<unsigned char>(peek()))) {
            if (identifier() != "pi") fail("unknown symbol in expression");
            return std::numbers::pi;
        }
        std::size_t consumed = 0;
        const double v = std::stod(std::string(src_.substr(pos_)), &consumed);
        if (consumed == 0) fail("expected number");
        pos_ += consumed;
        return v;
    }

    // —— 寄存器引用 q[i]；whole=true 时允许不带下标的整寄存器 ——
    std::optional<std::size_t> parse_ref(const std::string& expected_reg) {
        const std::string reg = identifier();
        if (reg != expected_reg) fail("unknown register '" + reg + "'");
        skip_ws();
        if (peek() != '[') return std::nullopt;  // 整寄存器
        expect('[');
        const auto idx = static_cast<std::size_t>(parse_expr());
        expect(']');
        return idx;
    }

    // —— 语句 ——
    void parse_statement() {
        const std::string head = identifier();
        if (head == "include") {
            skip_ws();
            expect('"');
            while (pos_ < src_.size() && src_[pos_] != '"') ++pos_;
            expect('"');
            expect(';');
            return;
        }
        if (head == "qreg" || head == "creg") {
            const std::string name = identifier();
            expect('[');
            const auto size = static_cast<std::size_t>(parse_expr());
            expect(']');
            expect(';');
            if (head == "qreg") {
                if (circuit_) fail("only one qreg is supported");
                qreg_name_ = name;
                qreg_size_ = size;
                make_circuit();
            } else {
                if (creg_size_ > 0) fail("only one creg is supported");
                if (circuit_ && circuit_->num_instructions() > 0) {
                    fail("creg must be declared before any gate");
                }
                creg_name_ = name;
                creg_size_ = size;
                circuit_.reset();  // 重建以带上经典寄存器
                make_circuit();
            }
            return;
        }
        if (!circuit_) fail("statement before qreg declaration");

        if (head == "measure") {
            const auto q = parse_ref(qreg_name_);
            skip_ws();
            expect_token("->");
            const auto c = parse_ref(creg_name_);
            expect(';');
            if (q && c) {
                circuit_->measure(static_cast<qubit_t>(*q),
                                  static_cast<clbit_t>(*c));
            } else if (!q && !c) {  // measure q -> c; 整寄存器
                if (qreg_size_ != creg_size_) {
                    fail("whole-register measure needs equal register sizes");
                }
                for (std::size_t i = 0; i < qreg_size_; ++i) {
                    circuit_->measure(static_cast<qubit_t>(i),
                                      static_cast<clbit_t>(i));
                }
            } else {
                fail("mixed whole-register and indexed measure");
            }
            return;
        }
        if (head == "barrier") {
            // 忽略参数列表，作为全局屏障
            while (pos_ < src_.size() && src_[pos_] != ';') ++pos_;
            expect(';');
            circuit_->barrier();
            return;
        }
        if (head == "reset") {
            const auto q = parse_ref(qreg_name_);
            expect(';');
            if (!q) fail("reset requires an indexed qubit");
            circuit_->reset(static_cast<qubit_t>(*q));
            return;
        }
        parse_gate(head);
    }

    void parse_gate(const std::string& name) {
        std::vector<double> params;
        skip_ws();
        if (consume('(')) {
            params.push_back(parse_expr());
            while (consume(',')) params.push_back(parse_expr());
            expect(')');
        }
        std::vector<qubit_t> qubits;
        do {
            const auto q = parse_ref(qreg_name_);
            if (!q) fail("gate operands must be indexed qubits");
            qubits.push_back(static_cast<qubit_t>(*q));
        } while (consume(','));
        expect(';');
        apply_gate(name, params, qubits);
    }

    void apply_gate(const std::string& name, const std::vector<double>& p,
                    const std::vector<qubit_t>& q) {
        QuantumCircuit& qc = *circuit_;
        // (门名, 参数个数, 比特数) → 电路调用
        if (name == "id" && q.size() == 1) {
            qc.i(q[0]);
            return;
        }
        if (name == "x" && q.size() == 1) {
            qc.x(q[0]);
            return;
        }
        if (name == "y" && q.size() == 1) {
            qc.y(q[0]);
            return;
        }
        if (name == "z" && q.size() == 1) {
            qc.z(q[0]);
            return;
        }
        if (name == "h" && q.size() == 1) {
            qc.h(q[0]);
            return;
        }
        if (name == "s" && q.size() == 1) {
            qc.s(q[0]);
            return;
        }
        if (name == "sdg" && q.size() == 1) {
            qc.sdg(q[0]);
            return;
        }
        if (name == "t" && q.size() == 1) {
            qc.t(q[0]);
            return;
        }
        if (name == "tdg" && q.size() == 1) {
            qc.tdg(q[0]);
            return;
        }
        if (name == "sx" && q.size() == 1) {
            qc.sx(q[0]);
            return;
        }
        if (name == "rx" && p.size() == 1 && q.size() == 1) {
            qc.rx(p[0], q[0]);
            return;
        }
        if (name == "ry" && p.size() == 1 && q.size() == 1) {
            qc.ry(p[0], q[0]);
            return;
        }
        if (name == "rz" && p.size() == 1 && q.size() == 1) {
            qc.rz(p[0], q[0]);
            return;
        }
        if ((name == "p" || name == "u1") && p.size() == 1 && q.size() == 1) {
            qc.p(p[0], q[0]);
            return;
        }
        if (name == "u2" && p.size() == 2 && q.size() == 1) {
            qc.u(std::numbers::pi / 2, p[0], p[1], q[0]);
            return;
        }
        if ((name == "u" || name == "u3") && p.size() == 3 && q.size() == 1) {
            qc.u(p[0], p[1], p[2], q[0]);
            return;
        }
        if (name == "cx" && q.size() == 2) {
            qc.cx(q[0], q[1]);
            return;
        }
        if (name == "cy" && q.size() == 2) {
            qc.cy(q[0], q[1]);
            return;
        }
        if (name == "cz" && q.size() == 2) {
            qc.cz(q[0], q[1]);
            return;
        }
        if (name == "ch" && q.size() == 2) {
            qc.ch(q[0], q[1]);
            return;
        }
        if ((name == "cp" || name == "cu1") && p.size() == 1 && q.size() == 2) {
            qc.cp(p[0], q[0], q[1]);
            return;
        }
        if (name == "crx" && p.size() == 1 && q.size() == 2) {
            qc.crx(p[0], q[0], q[1]);
            return;
        }
        if (name == "cry" && p.size() == 1 && q.size() == 2) {
            qc.cry(p[0], q[0], q[1]);
            return;
        }
        if (name == "crz" && p.size() == 1 && q.size() == 2) {
            qc.crz(p[0], q[0], q[1]);
            return;
        }
        if (name == "swap" && q.size() == 2) {
            qc.swap(q[0], q[1]);
            return;
        }
        if (name == "ccx" && q.size() == 3) {
            qc.ccx(q[0], q[1], q[2]);
            return;
        }
        if (name == "cswap" && q.size() == 3) {
            qc.cswap(q[0], q[1], q[2]);
            return;
        }
        fail("unsupported gate '" + name + "' with " +
             std::to_string(p.size()) + " params / " +
             std::to_string(q.size()) + " qubits");
    }

    void make_circuit() {
        if (qreg_size_ == 0) return;
        circuit_ = std::make_unique<QuantumCircuit>(qreg_size_, creg_size_);
    }

    std::string_view src_;
    std::size_t pos_ = 0;
    std::unique_ptr<QuantumCircuit> circuit_;
    std::string qreg_name_, creg_name_;
    std::size_t qreg_size_ = 0, creg_size_ = 0;
};

}  // namespace

QuantumCircuit from_qasm(std::string_view source) {
    return Parser(source).parse();
}

}  // namespace lqcs::io
