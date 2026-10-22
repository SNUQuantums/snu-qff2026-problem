"""The judge's own checks. `pytest tests -q` must pass before a release."""

import math
from pathlib import Path

import numpy as np
import pytest
import xtim
from qiskit import QuantumCircuit

import tprep
from tprep.check import CHUNK, _sample, final_round
from tprep.qasm import load, to_stim

ROOT = Path(__file__).resolve().parent.parent
INJECTION = ROOT / "examples" / "injection.qasm"
INJECTION_DATA = [0, 6, 13, 2, 7, 8, 1]
HEADER = 'OPENQASM 2.0;\ninclude "qelib1.inc";\n'


@pytest.fixture(scope="module")
def injection():
    return INJECTION.read_text()


def qasm(body, n=36, creg=0):
    return HEADER + f"qreg q[{n}];\n" + (f"creg c[{creg}];\n" if creg else "") + body


def format_errors(text):
    return tprep.check(text, INJECTION_DATA).format_errors


def test_grid_is_square_and_connected():
    assert tprep.N_QUBITS == tprep.ROWS * tprep.COLS
    assert len(tprep.couplers()) == tprep.ROWS * (tprep.COLS - 1) + tprep.COLS * (tprep.ROWS - 1)
    assert tprep.adjacent(0, 1) and tprep.adjacent(0, tprep.COLS)
    assert not tprep.adjacent(0, tprep.COLS + 1)  # diagonal is not a coupler


def test_format_rejects_bad_circuits():
    assert "OpenQASM 2.0" in format_errors("H 0\nM 0")[0]                     # Stim is not accepted
    assert "not neighbours" in format_errors(qasm("cx q[0],q[35];"))[0]
    assert "not an allowed operation" in format_errors(qasm("ccx q[0],q[1],q[2];"))[0]
    assert "not an allowed operation" in format_errors(qasm("rx(0.3) q[0];"))[0]
    assert "not an allowed operation" in format_errors(qasm("gate g a { h a; }\ng q[0];"))[0]
    assert "not an allowed operation" in format_errors(
        qasm("measure q[0] -> c[0];\nif (c==1) x q[1];", creg=1))[0]
    assert "chip has 36" in format_errors(qasm("h q[0];", n=40))[0]
    assert "exactly one quantum register" in format_errors(HEADER + "qreg a[2];\nqreg b[2];\nh a[0];")[0]
    assert "not a measured classical bit" in format_errors(qasm("h q[0];\n// DETECTOR c[0]", creg=1))[0]
    assert "could not be parsed" in format_errors(qasm("h q[0]"))[0]
    assert "empty" in format_errors(qasm(""))[0]


def test_conversion_to_the_stim_the_judge_runs():
    text = qasm("h q[0];\ncs q[0],q[1];\ncsdg q[1],q[2];\nch q[2],q[3];\nsx q[0];\n"
                "measure q[0] -> c[0];\nmeasure q[1] -> c[1];\nreset q[0];\nmeasure q[2] -> c[0];\n"
                "// DETECTOR c[0] c[1]\n// DETECTOR c[1]", creg=2)
    assert "cs is not allowed" in to_stim(text)[1][0]      # the official gate set is CNOT only
    stim, errors = to_stim(text, noise=tprep.LEGACY)      # judge 2.x's gate set
    assert errors == []
    # c[0] means its last write; a detector goes right after its last measurement
    assert stim.splitlines() == ["H 0", "CS 0 1", "CS_DAG 1 2", "CH 2 3", "SQRT_X 0", "M 0", "M 1",
                                 "DETECTOR rec[-1]", "R 0", "M 2", "DETECTOR rec[-1] rec[-2]"]


def test_x_readout_runs_as_mx_only_when_the_qubit_is_not_reused():
    text = qasm("h q[0];\nmeasure q[0] -> c[0];\nh q[1];\nmeasure q[1] -> c[1];\nreset q[1];", creg=2)
    assert to_stim(text)[0].splitlines() == ["MX 0", "MX 1", "R 1"]    # native X readout: no h
    stim, _ = to_stim(text, noise=tprep.LEGACY)                         # h kept as I, with its noise
    assert stim.splitlines() == ["I 0", "MX 0", "I 1", "MX 1", "R 1"]
    # measured again (the counterexample to rewriting it always): kept as H + M
    stim, _ = to_stim(qasm("h q[0];\nmeasure q[0] -> c[0];\nmeasure q[0] -> c[1];\n// DETECTOR c[0] c[1]", creg=2))
    assert stim.splitlines() == ["H 0", "M 0", "M 0", "DETECTOR rec[-2] rec[-1]"]
    # h, measure, h leaves |+>/|->, as MX does: run as I, MX, I (the same channel, noise included)
    text = qasm("h q[0];\nt q[0];\nh q[0];\nmeasure q[0] -> c[0];\nh q[0];\nmeasure q[0] -> c[1];\n", creg=2)
    assert to_stim(text)[0].splitlines()[2:] == ["MX 0", "M 0"]
    stim, _ = to_stim(text, noise=tprep.LEGACY)
    assert stim.splitlines()[2:] == ["I 0", "MX 0", "I 0", "M 0"]
    other = text.replace("measure q[0] -> c[0];\nh q[0];", "measure q[0] -> c[0];\nx q[0];")
    assert to_stim(other)[0].splitlines()[2:5] == ["H 0", "M 0", "X 0"]   # next is not h: kept


def test_submission_text_round_trips_detectors():
    qc = QuantumCircuit(3, 2)
    qc.h(0)
    qc.measure([0, 1], [0, 1])
    text = tprep.submission_text(qc, [[qc.clbits[0], qc.clbits[1]], [("c", 0)]])
    assert text.splitlines()[-2:] == ["// DETECTOR c[0] c[1]", "// DETECTOR c[0]"]
    assert to_stim(text)[0].splitlines() == ["MX 0", "DETECTOR rec[-1]", "M 1", "DETECTOR rec[-2] rec[-1]"]


def test_multiqubit_gate_count():
    assert tprep.count_multiqubit_gates("CX 0 1 1 2\nH 0\nCZ 0 6") == 3


def test_noise_insertion():
    """The official model (problem.md 1.5): uniform depolarizing noise at p."""
    out = tprep.apply_noise("H 0\nCX 0 1\nM 2\nR 2\nMX 0\nRX 1", 1e-3).splitlines()
    assert out == [
        "X_ERROR(0.001) 0", "H 0", "DEPOLARIZE1(0.001) 0",          # first use: a noisy |0>
        "X_ERROR(0.001) 1", "CX 0 1", "DEPOLARIZE2(0.001) 0 1",
        "X_ERROR(0.001) 2", "X_ERROR(0.001) 2", "M 2", "DEPOLARIZE1(0.001) 2",   # |0>, flip, after
        "R 2", "X_ERROR(0.001) 2",               # M then R: no wait before a reset
        "DEPOLARIZE1(0.0009999999999999731) 2",  # q[2] was reset at step 1 and waits step 2
        "Z_ERROR(0.001) 0", "MX 0", "DEPOLARIZE1(0.001) 0",
        "RX 1", "Z_ERROR(0.001) 1"]
    legacy = tprep.apply_noise("H 0\nCX 0 1\nM 2\nR 2", 1e-3, tprep.LEGACY).splitlines()
    assert legacy == ["H 0", "DEPOLARIZE1(0.001) 0", "CX 0 1", "DEPOLARIZE2(0.001) 0 1",
                      "X_ERROR(0.001) 2", "M 2", "R 2", "X_ERROR(0.001) 2"]
    assert tprep.apply_noise("H 0", 1e-3, tprep.NoiseModel(gate1=0)) == "H 0"


def test_idle_noise_option():
    """Off in the official model; IDLE charges DEPOLARIZE1(p) per waiting time step."""
    text = "H 0\nCX 0 1\nCX 1 2\nCX 2 3\nM 3\nDETECTOR rec[-1]\nH 0"
    assert "DEPOLARIZE1(0.001) 1" not in tprep.apply_noise(text, 1e-3)
    out = tprep.apply_noise(text, 1e-3, tprep.IDLE).splitlines()
    idle = {(i, l.split()[1]): float(l[12:].split(")")[0]) for i, l in enumerate(out)
            if l.startswith("DEPOLARIZE1") and out[i - 1].split()[0] != "H"}
    two_steps = 0.75 * (1 - (1 - 4e-3 / 3) ** 2)
    # steps: H0 0 | CX01 1 | CX12 2, H0 2 | CX23 3 | M3 4. q1 waits 3,4; q2 waits 4; q0 waits 3,4.
    assert idle == {(6, "1"): pytest.approx(two_steps), (9, "2"): pytest.approx(1e-3),
                    (15, "0"): pytest.approx(two_steps)}
    assert [l for l in out if l.startswith(("H", "CX", "M", "DETECTOR"))] == text.splitlines()


def test_final_round_corrects_one_error_and_not_two(injection):
    """The judge's noiseless round fixes any weight-1 error on the output; weight 2 becomes logical."""
    stim = to_stim(injection)[0]
    d = INJECTION_DATA
    tail = final_round(d) + "".join(
        f"\nPAULI_EXPECTATION({k}) " + "*".join(f"{P}{q}" for q in d) for k, P in enumerate("XY"))
    run = lambda err: _sample(xtim.Circuit(stim + err + "\n" + tail), 100, 0)["score"]
    for err in (f"X_ERROR(1) {d[3]}", f"Z_ERROR(1) {d[4]}", f"Y_ERROR(1) {d[6]}"):
        assert run(err) == pytest.approx(0, abs=1e-12)
    assert run(f"X_ERROR(1) {d[3]} {d[4]}") == pytest.approx(0.5)   # X-bar: Y-bar flips, X-bar does not
    assert run(f"Z_ERROR(1) {d[3]} {d[4]}") == pytest.approx(1.0)   # Z-bar: both flip


def test_correctness_wants_exactly_t_in_the_code_space(injection):
    first = lambda text: tprep.check(text, INJECTION_DATA, shots=10**4).correctness_errors[0]
    assert "does not prepare" in first(injection.replace("t q[0];", "tdg q[0];", 1))
    assert "code space" in first(injection + "x q[6];\n")
    random_det = injection.replace("qreg q[36];", "qreg q[36];\ncreg r[1];") + \
        "h q[30];\nmeasure q[30] -> r[0];\n// DETECTOR r[0]\n"
    assert "not deterministic" in first(random_det)
    always_one = random_det.replace("h q[30];", "x q[30];")      # parity 1 every run: a fine detector
    assert tprep.check(always_one, INJECTION_DATA, shots=10**4).passed
    assert "not |T_L>" in first(qasm("h q[0];\nt q[0];"))


def test_pauli_frame_from_measurements_is_tracked(injection):
    """Z_L^m, X_L^m or Y_L^m times |T_L>, with m measured, is |T_L> up to a known Pauli: the judge
    finds m and corrects for it (problem.md section 2). With m reset unread, the output is mixed."""
    head, rest = injection.split("t q[0];\n", 1)
    cz = "h q[0];\ncx q[1],q[0];\nh q[0];\n"
    for gates in (cz, "cx q[1],q[0];\n", cz + "cx q[1],q[0];\n"):  # Z, X, Y
        byproduct = "h q[1];\n" + gates                    # q[1] is still unused here
        read = byproduct + "measure q[1] -> f[0];\nreset q[1];\n"
        text = head.replace("qreg q[36];", "qreg q[36];\ncreg f[1];") + "t q[0];\n" + read + rest
        r = tprep.check(text, INJECTION_DATA, p=0.0, shots=10**4)
        assert r.correct_ok and r.score == pytest.approx(0, abs=1e-12), (gates, str(r))
        unread = text.replace("measure q[1] -> f[0];\n", "")
        assert not tprep.check(unread, INJECTION_DATA, p=0.0, shots=10**4).correct_ok, gates


def test_injection_example_passes(injection):
    r = tprep.check(injection, INJECTION_DATA, shots=10**5)
    assert r.passed, str(r)
    assert r.n_2q == 28
    assert r.acceptance == pytest.approx(1.0)


def test_quantum_circuit_and_text_score_alike(injection):
    assert tprep.check(load(injection), INJECTION_DATA, shots=10**4).score == \
        tprep.check(injection, INJECTION_DATA, shots=10**4).score


def test_injection_is_exact_without_noise(injection):
    r = tprep.check(injection, INJECTION_DATA, p=0.0, shots=10**5)
    assert r.score == 0.0, "the example must prepare |T_L> exactly when noiseless"


def test_score_is_first_order_in_p(injection):
    """Injection cannot beat a c1*p floor: the error it injects is a logical one."""
    scores = [tprep.check(injection, INJECTION_DATA, p=p, shots=10**5).score for p in (1e-4, 1e-3)]
    assert scores[1] / scores[0] == pytest.approx(10, rel=0.15)


def test_score_is_deterministic(injection):
    a = tprep.check(injection, INJECTION_DATA, shots=10**5).score
    b = tprep.check(injection, INJECTION_DATA, shots=10**5).score
    assert a == b


def test_shots_run_in_seeded_chunks(injection):
    """Chunk k uses seed + k, so the first 1e6 shots of the official run are the 1e6-shot run."""
    a = tprep.check(injection, INJECTION_DATA, shots=2 * CHUNK)
    b = tprep.check(injection, INJECTION_DATA, shots=CHUNK, seed=0)
    c = tprep.check(injection, INJECTION_DATA, shots=CHUNK, seed=1)
    assert a.score == pytest.approx((b.score + c.score) / 2, rel=1e-9)


def test_wrong_placement_fails_correctness(injection):
    r = tprep.check(injection, [0, 1, 2, 3, 4, 5, 6], shots=10**5)
    assert not r.passed
    assert r.format_ok and not r.correct_ok


def test_data_qubits_must_be_seven_distinct(injection):
    assert not tprep.check(injection, [0, 0, 13, 2, 7, 8, 1]).format_ok
    assert not tprep.check(injection, [0, 6, 13]).format_ok
    assert not tprep.check(injection, INJECTION_DATA + [1]).format_ok        # 8 entries, one repeated
    assert not tprep.check(injection, [0, 6, 13, 2, 7, 8, 1.0]).format_ok    # not an integer
    assert not tprep.check(injection, [0, 6, 13, 2, 7, 8, True]).format_ok
    assert not tprep.check(injection, [0, 6, 13, 2, 7, 8, 36]).format_ok     # off the chip


def test_size_limits(injection):
    from tprep.qasm import MAX_BYTES, MAX_STATEMENTS
    assert "statements" in format_errors(qasm("id q[0];\n" * (MAX_STATEMENTS + 1)))[0]
    assert "bytes" in format_errors(qasm("// " + "x" * MAX_BYTES + "\nh q[0];"))[0]
    from tprep.qasm import MAX_MEASURE_RESET
    assert "measure and reset" in format_errors(qasm("reset q[0];\n" * (MAX_MEASURE_RESET + 1)))[0]
    assert tprep.check(qasm("id q[0];\n" * 10 + "// a; comment; with; semicolons\n"), INJECTION_DATA).format_errors \
        == []  # passes format; it fails later for not preparing |T_L>


def test_barrier_delays_a_qubit_for_idle_noise():
    """A barrier holds its qubits as in Qiskit, so a fresh qubit can start late and not wait."""
    only_idle = tprep.NoiseModel(gate1=0, gate2=0, meas=0, reset=0, idle=1)
    body = "h q[0];\ncx q[0],q[1];\ncx q[0],q[1];\ncx q[0],q[1];\n{}reset q[2];\ncx q[1],q[2];\n"

    def idle(barrier):
        stim, errors = to_stim(qasm(body.format(barrier)))
        assert not errors
        out = tprep.apply_noise(stim, 1e-3, only_idle).splitlines()
        return {l.split()[1]: float(l[12:].split(")")[0]) for l in out if l.startswith("DEPOLARIZE1")}

    steps = lambda n: pytest.approx(0.75 * (1 - (1 - 4e-3 / 3) ** n))
    # No barrier: reset q2 at step 0, cx q1,q2 at step 4, so q2 waits steps 1-3 (q0 waits step 4).
    assert idle("") == {"2": steps(3), "0": steps(1)}
    # With it: reset q2 at step 4 once q1 is free, cx at 5; q1 waits step 4 instead (q0 steps 4-5).
    assert idle("barrier q[1],q[2];\n") == {"1": steps(1), "0": steps(2)}
    legacy = lambda b: tprep.apply_noise(to_stim(qasm(body.format(b)))[0], 1e-3, tprep.LEGACY)
    assert legacy("barrier q[1],q[2];\n") == legacy("")   # no effect without idle noise


def test_uniform_noise_model():
    """CNOT only; h + measure and reset + h are native X-basis operations; measurements depolarize."""
    stim, errors = to_stim(qasm("reset q[0];\nh q[0];\ncx q[0],q[1];\nh q[0];\nmeasure q[0] -> c[0];\n",
                                creg=1), noise=tprep.UNIFORM)
    assert not errors and stim.split("\n")[:3] == ["RX 0", "CX 0 1", "MX 0"]
    assert "DEPOLARIZE1(0.001) 0" in tprep.apply_noise("M 0", 1e-3, tprep.UNIFORM)
    assert "swap is not allowed" in to_stim(qasm("swap q[0],q[1];\n"), noise=tprep.UNIFORM)[1][0]


# A |T> made by measurement: |+>, then T^dag, an X readout through q[1] (m[0]), then T leaves
# Z^m |T>, m = 1 with probability sin^2(pi/8). The injection example then encodes it, so the output
# is Z_L^m |T_L>: correct with FRAME_Z on m[0], or with POSTSELECT(0) m[0].
MEASURED_T = ("h q[0];\ntdg q[0];\nh q[1];\ncx q[1],q[0];\nh q[1];\nmeasure q[1] -> m[0];\n"
              "reset q[1];\nt q[0];\n")
COS2 = math.cos(math.pi / 8) ** 2


def measured_t(injection, tail="", creg=2):
    head = injection.replace("qreg q[36];", f"qreg q[36];\ncreg m[{creg}];")
    return head.replace("h q[0];\nt q[0];\n", MEASURED_T, 1) + tail


def test_declared_frame_corrects_a_magic_measurement(injection):
    noiseless = lambda text: tprep.check(text, INJECTION_DATA, p=0.0, shots=10**4)
    r = noiseless(measured_t(injection, "// FRAME_X(0)\n// FRAME_Z(0) m[0]\n"))
    assert r.correct_ok and r.score == pytest.approx(0, abs=1e-12), str(r)
    assert r.ideal_acceptance == pytest.approx(1) and r.frame == {"X": (0, []), "Z": (0, [0])}
    # The judge's own frame search misses this frame (the outcome is not uniformly random).
    errors = noiseless(measured_t(injection)).correctness_errors
    assert "does not prepare" in errors[0] and "declare it" in errors[1]
    assert errors[0].endswith(f"{1 - COS2:.6g}")
    wrong = noiseless(measured_t(injection, "// FRAME_X(0) m[0]\n// FRAME_Z(0)\n"))
    assert "does not prepare" in wrong.correctness_errors[0]


def test_postselection_may_be_random(injection):
    noiseless = lambda tail: tprep.check(measured_t(injection, tail), INJECTION_DATA, p=0.0, shots=10**4)
    r = noiseless("// POSTSELECT(0) m[0]\n")              # keep m = 0, where no frame is needed
    assert r.correct_ok and r.ideal_acceptance == pytest.approx(COS2, abs=1e-12), str(r)
    assert r.acceptance == pytest.approx(COS2, abs=0.02) and r.score == pytest.approx(0, abs=1e-12)
    # The constant term: keep m = 1 and correct it with FRAME_Z(1). Correct, but under the floor.
    r = noiseless("// POSTSELECT(1) m[0]\n// FRAME_X(0)\n// FRAME_Z(1)\n")
    assert r.correct_ok and r.ideal_acceptance == pytest.approx(1 - COS2, abs=1e-12)
    assert not r.acceptance_ok
    # Keeping m = 1 without the correction: the accepted runs are Z_L |T_L>.
    assert "does not prepare" in noiseless("// POSTSELECT(1) m[0]\n").correctness_errors[0]


def test_conditional_detector(injection):
    """A second readout of the same magic observable (m[1]) equals m[0]. Alone it is random, but
    among the runs POSTSELECT(0) m[0] keeps it is always 0: a detector there, read against 0."""
    again = "h q[1];\ncx q[1],q[0];\nh q[1];\nmeasure q[1] -> m[1];\nreset q[1];\n"
    text = measured_t(injection).replace("reset q[1];\nt q[0];", "reset q[1];\n" + again + "t q[0];", 1)
    run = lambda tail, p=0.0: tprep.check(text + tail, INJECTION_DATA, p=p, shots=10**5)
    r = run("// POSTSELECT(0) m[0]\n// DETECTOR m[1]\n")
    assert r.correct_ok and r.conditional_detectors == [0], str(r)
    assert r.ideal_acceptance == pytest.approx(COS2) and r.acceptance == pytest.approx(COS2, abs=0.01)
    # Without the POSTSELECT it is a random detector, as before.
    assert "not deterministic" in run("// DETECTOR m[1]\n// FRAME_X(0)\n// FRAME_Z(0) m[0]\n").correctness_errors[0]
    # It catches errors: with noise, the runs where m[1] = 1 are dropped.
    noisy = run("// POSTSELECT(0) m[0]\n// DETECTOR m[1]\n", p=1e-2)
    plain = run("// POSTSELECT(0) m[0]\n", p=1e-2)
    assert noisy.acceptance < plain.acceptance and noisy.score < plain.score
    # The same comparison written as an ordinary (unconditional) detector stays on xtim's path.
    r = run("// DETECTOR m[0] m[1]\n// FRAME_X(0)\n// FRAME_Z(0) m[0]\n")
    assert r.correct_ok and r.conditional_detectors == [] and r.ideal_acceptance == pytest.approx(1)


def test_only_the_accepted_runs_must_be_right(injection):
    """A random X on a data qubit, heralded by m[0]: wrong (outside the code space) unless the
    runs that carry it are postselected away."""
    head = injection.replace("qreg q[36];", "qreg q[36];\ncreg m[1];")
    text = head + "h q[3];\ncx q[3],q[2];\nmeasure q[3] -> m[0];\n"
    bad = tprep.check(text, INJECTION_DATA, p=0.0, shots=10**4)
    assert "code space" in bad.correctness_errors[0]
    good = tprep.check(text + "// POSTSELECT(0) m[0]\n", INJECTION_DATA, p=0.0, shots=10**4)
    assert good.correct_ok and good.ideal_acceptance == pytest.approx(0.5), str(good)


def test_declarations_read_the_last_measurement_of_a_bit(injection):
    """As for DETECTOR: a bit measured twice means its last measurement."""
    overwrite = "reset q[2];\nmeasure q[2] -> m[0];\n"     # q[2] is unused until later; reads 0
    text = measured_t(injection).replace("reset q[1];\n", "reset q[1];\n" + overwrite, 1)
    lost = tprep.check(text + "// FRAME_X(0)\n// FRAME_Z(0) m[0]\n", INJECTION_DATA, p=0.0, shots=10**4)
    assert "does not prepare" in lost.correctness_errors[0]
    kept = text.replace("measure q[1] -> m[0];", "measure q[1] -> m[1];")
    r = tprep.check(kept + "// FRAME_X(0)\n// FRAME_Z(0) m[1]\n", INJECTION_DATA, p=0.0, shots=10**4)
    assert r.correct_ok, str(r)


def test_declared_frame_matches_the_found_one(injection):
    """Existing submissions keep their score: declaring the frame the judge finds changes nothing."""
    head, rest = injection.split("t q[0];\n", 1)
    text = head.replace("qreg q[36];", "qreg q[36];\ncreg f[1];") + \
        "t q[0];\nh q[1];\ncx q[1],q[0];\nmeasure q[1] -> f[0];\nreset q[1];\n" + rest   # X_L^f
    found = tprep.check(text, INJECTION_DATA, shots=10**4)
    declared = tprep.check(text + "// FRAME_X(0) f[0]\n// FRAME_Z(0)\n", INJECTION_DATA, shots=10**4)
    assert found.frame == declared.frame == {"X": (0, [0]), "Z": (0, [])}
    assert found.score == declared.score and found.acceptance == declared.acceptance


def test_declaration_format(injection):
    errors = lambda tail: tprep.check(measured_t(injection, tail), INJECTION_DATA, shots=10**3).format_errors
    assert "exactly once each" in errors("// FRAME_Z(0) m[0]\n")[0]
    assert "exactly once each" in errors("// FRAME_X(0)\n// FRAME_X(0)\n// FRAME_Z(0)\n")[0]
    assert "names no bits" in errors("// POSTSELECT(0)\n")[0]
    assert "not a measured classical bit" in errors("// POSTSELECT(0) m[1]\n")[0]
    assert "FRAME_Z(0) or FRAME_Z(1)" in errors("// FRAME_X(0)\n// FRAME_Z(2) m[0]\n")[0]


def test_quantum_circuit_declarations_match_text(injection):
    text = measured_t(injection)
    qc = load(text)
    m = qc.cregs[0]
    from_qc = tprep.check(qc, INJECTION_DATA, frame={"X": (0, []), "Z": (0, [m[0]])},
                          postselect=[(0, [("m", 0)])], shots=10**4)
    from_text = tprep.check(text + "// POSTSELECT(0) m[0]\n// FRAME_X(0)\n// FRAME_Z(0) m[0]\n",
                            INJECTION_DATA, shots=10**4)
    assert from_qc.correct_ok and (from_qc.score, from_qc.acceptance) == (from_text.score, from_text.acceptance)
    with pytest.raises(TypeError):
        tprep.check(text, INJECTION_DATA, postselect=[(0, [("m", 0)])])


def test_exact_check_matches_a_statevector():
    """tprep.exact against a 2-qubit statevector, for every non-Clifford gate it rewrites."""
    import itertools

    import stim
    from tprep.exact import run

    I2, X, Y, Z = np.eye(2), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])
    H, S, T = (X + Z) / np.sqrt(2), np.diag([1, 1j]), np.diag([1, np.exp(1j * np.pi / 4)])
    one = {"H": H, "T": T, "T_DAG": T.conj().T, "S": S}
    ctrl = lambda U: np.kron(np.diag([1, 0]), I2) + np.kron(np.diag([0, 1]), U)
    two = {"CS": ctrl(S), "CS_DAG": ctrl(S.conj().T), "CH": ctrl(H)}
    for text in ("H 0\nT 0", "H 0\nH 1\nCS 0 1\nT_DAG 1", "H 0\nCS_DAG 0 1\nH 1\nCH 1 0", "H 0\nCH 0 1",
                 "H 0\nT 0\nH 0\nT 0\nS 0"):
        psi = np.array([1, 0, 0, 0], complex)
        for line in text.splitlines():
            name, *qs = line.split()
            if name in two:   # matrices above have control = qubit 0; swap for `CH 1 0`
                SW = np.eye(4)[[0, 2, 1, 3]]
                U = two[name] if qs == ["0", "1"] else SW @ two[name] @ SW
            else:
                U = np.kron(*[one[name] if q == int(qs[0]) else I2 for q in (0, 1)])
            psi = U @ psi
        st, inv, _, _, wire = run(text, 2)
        for paulis in itertools.product("IXYZ", repeat=2):
            M = np.kron(*[{"I": I2, "X": X, "Y": Y, "Z": Z}[c] for c in paulis])
            ps = stim.PauliString(st.n)
            for q, c in enumerate(paulis):
                if c != "I":
                    ps[wire[q]] = c
            assert st.expectation(ps, inv) == pytest.approx((psi.conj() @ M @ psi).real, abs=1e-12), (text, paulis)


def test_fault_scan_with_postselection(injection):
    again = "h q[1];\ncx q[1],q[0];\nh q[1];\nmeasure q[1] -> m[1];\nreset q[1];\n"
    text = measured_t(injection, "// POSTSELECT(0) m[0]\n// DETECTOR m[1]\n").replace(
        "reset q[1];\nt q[0];", "reset q[1];\n" + again + "t q[0];", 1)          # with a conditional detector
    scan = tprep.fault_scan(text, INJECTION_DATA, processes=1)
    slope = tprep.check(text, INJECTION_DATA, p=1e-5, shots=10**6).score / 1e-5
    assert scan.c1 == pytest.approx(slope, rel=0.1)


def test_fault_scan_matches_the_score(injection):
    """c1 from single faults is the slope of 1-F at small p."""
    scan = tprep.fault_scan(injection, INJECTION_DATA, processes=1)
    slope = tprep.check(injection, INJECTION_DATA, p=1e-5, shots=10**6).score / 1e-5
    assert scan.c1 == pytest.approx(slope, rel=0.1)
    assert "instruction" in str(scan)


def test_fault_scan_matches_the_score_with_idle_noise(injection):
    scan = tprep.fault_scan(injection, INJECTION_DATA, noise=tprep.IDLE, processes=1)
    slope = tprep.check(injection, INJECTION_DATA, p=1e-5, shots=10**6, noise=tprep.IDLE).score / 1e-5
    assert scan.c1 == pytest.approx(slope, rel=0.1)


def test_draw_returns_a_figure(injection):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    assert tprep.draw(injection, INJECTION_DATA).axes
