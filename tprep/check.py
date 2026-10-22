"""The four validity gates and the score (problem.md section 2)."""

import numbers
import warnings
from dataclasses import dataclass, field

import numpy as np

from .grid import ACCEPTANCE_FLOOR, N_CODE, N_QUBITS, SUPPORTS
from .noise import STANDARD, apply_noise
from .parse import check_format, count_measurements, count_multiqubit_gates, declarations, parse

CHUNK = 10**6  # shots per sampler; chunk k runs with seed + k (problem.md section 2.3)
# Target expectations of |T_L>. Y-bar = i X-bar Z-bar = -Y^7 (problem.md section 0.3).
TARGET = {"X": 2**-0.5, "Y": -(2**-0.5)}


def final_round(data_qubits):
    """The judge's noiseless syndrome round on the output: X checks, then Z checks (6 records)."""
    return "\n".join(f"MPP {'*'.join(f'{P}{data_qubits[i]}' for i in sorted(s))}"
                     for P in "XZ" for s in SUPPORTS)


@dataclass
class Report:
    p: float
    shots: int
    seed: int
    format_errors: list = field(default_factory=list)
    class_error: str | None = None
    correctness_errors: list = field(default_factory=list)
    score: float | None = None        # infidelity 1 - F after the ideal final correction
    sem: float | None = None          # its Monte Carlo standard error
    acceptance: float | None = None
    ideal_acceptance: float | None = None   # exact, noiseless; below 1 only with POSTSELECT
    frame: dict | None = None         # the logical frame used: declared, or found by the judge
    conditional_detectors: list = field(default_factory=list)  # 0-based; deterministic given POSTSELECT
    n_2q: int | None = None
    diagnosis: object = None

    @property
    def format_ok(self):
        return not self.format_errors

    @property
    def class_ok(self):
        return self.format_ok and self.class_error is None

    @property
    def correct_ok(self):
        return self.class_ok and not self.correctness_errors

    @property
    def acceptance_ok(self):
        return self.correct_ok and self.acceptance is not None and self.acceptance >= ACCEPTANCE_FLOOR

    @property
    def passed(self):
        return self.acceptance_ok

    def __str__(self):
        mark = lambda ok: "PASS" if ok else "FAIL"
        lines = [f"p={self.p:g}  shots={self.shots:,}  seed={self.seed}",
                 f"  format:      {mark(self.format_ok)}"]
        lines += [f"    - {e}" for e in self.format_errors]
        if self.format_ok:
            lines.append(f"  class:       {mark(self.class_ok)}")
            if self.class_error:
                lines.append(f"    - {self.class_error}")
        if self.class_ok:
            lines.append(f"  correctness: {mark(self.correct_ok)}")
            lines += [f"    - {e}" for e in self.correctness_errors]
        if self.correct_ok:
            lines.append(f"  acceptance:  {mark(self.acceptance_ok)}"
                         f"  ({self.acceptance:.4f}, floor {ACCEPTANCE_FLOOR})")
            if self.ideal_acceptance is not None and self.ideal_acceptance < 1 - 1e-9:
                lines.append(f"    (noiselessly {self.ideal_acceptance:.6f}: POSTSELECT keeps only part of the runs)")
            if self.conditional_detectors:
                lines.append("    (detector(s) " + ", ".join(f"#{k + 1}" for k in self.conditional_detectors) +
                             " are deterministic only among the POSTSELECT-ed runs, and are read against "
                             "their value there)")
            lines.append(f"  score 1-F = {self.score:.6g} ± {self.sem:.2g}   (3 significant figures: {self.score:.3g})")
            lines.append(f"  2q+ gates = {self.n_2q}")
        lines.append(f"  RESULT: {mark(self.passed)}")
        return "\n".join(lines)


def check(circuit, data_qubits, *, detectors=(), frame=None, postselect=(), p=1e-3, shots=10**6,
          seed=0, noise=STANDARD, n_qubits=N_QUBITS):
    """Grade one submission. This is the code the scoring server runs.

    `circuit` is the submission text (OpenQASM 2.0 with `// DETECTOR` lines, problem.md section 1.3,
    and optional `// POSTSELECT`, `// FRAME_X`, `// FRAME_Z` lines, see tprep.qasm) or a Qiskit
    QuantumCircuit, whose detectors, frame and postselect conditions are then given here (see
    `submission_text`). `shots` defaults to 1e6 for fast local iteration; the official score uses
    1e8 and seed 0.
    """
    from .qasm import submission_text, to_stim

    if not isinstance(circuit, str):
        circuit = submission_text(circuit, detectors, frame, postselect)
    elif frame is not None or postselect:
        raise TypeError("for submission text, write FRAME_X/FRAME_Z/POSTSELECT lines in the text")
    stim_text, errors = to_stim(circuit, n_qubits, noise)
    if errors:
        report = Report(p=p, shots=shots, seed=seed)
        report.format_errors = errors
        return report
    return check_stim(stim_text, data_qubits, p=p, shots=shots, seed=seed, noise=noise, n_qubits=n_qubits)


def check_stim(circuit_text, data_qubits, *, p=1e-3, shots=10**6, seed=0,
               noise=STANDARD, n_qubits=N_QUBITS):
    """`check` for a circuit already in the Stim text xtim reads. Submissions go through `check`.
    FRAME_X/FRAME_Z/POSTSELECT arrive as `#!` lines with absolute record indices (parse.declarations).
    """
    import xtim

    report = Report(p=p, shots=shots, seed=seed)

    # Gate 1 - format.
    report.format_errors = check_format(circuit_text, n_qubits=n_qubits)
    data_qubits = list(data_qubits)
    integers = all(isinstance(q, numbers.Integral) and not isinstance(q, bool) for q in data_qubits)
    if not (len(data_qubits) == N_CODE and integers and len(set(data_qubits)) == N_CODE
            and all(0 <= q < n_qubits for q in data_qubits)):
        report.format_errors.append(
            f"data_qubits must be exactly {N_CODE} distinct integers in 0..{n_qubits - 1}, got {data_qubits}")
    else:
        data_qubits = [int(q) for q in data_qubits]
    if noise.gates is not None:
        from .qasm import GATES
        allowed = {GATES[g] for g in noise.gates if g in GATES} | {"MX", "RX", "DETECTOR", "TICK"}
        report.format_errors += [f"line {ins.line_no}: {ins.name} is not allowed in this noise model"
                                 for ins in parse(circuit_text) if ins.name not in allowed]
    if not report.format_ok:
        return report
    report.n_2q = count_multiqubit_gates(circuit_text)

    n_meas = count_measurements(circuit_text)
    try:
        frame, post = declarations(circuit_text)
    except ValueError as e:
        report.format_errors.append(str(e))
        return report
    for _, recs in (list(frame.values()) if frame else []) + post:
        if any(not 0 <= r < n_meas for r in recs):
            report.format_errors.append(f"a FRAME or POSTSELECT record is outside 0..{n_meas - 1}")
            return report

    # Every noiseless analysis below runs on the circuit with its noise lines in place at strength 0.
    # xtim's preprocessing of a circuit with no noise lines at all can misjudge which records fix
    # the output's sign, and which detectors are deterministic, when an X-basis readout feeding
    # the Pauli frame is written as h + measure; with the lines present it is correct.
    noiseless = apply_noise(circuit_text, 0.0, noise)

    # Gate 2 - simulable class.
    try:
        clean = xtim.Circuit(noiseless)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            report.diagnosis = xtim.diagnose(clean)
    except Exception as e:                       # XtimRejectError, XtimParseError, ...
        report.class_error = f"{type(e).__name__}: {e}"
        return report

    # Gate 3 - correctness, exact over every pattern of measurement outcomes (tprep.exact).
    from .exact import ExactError, verify

    found = frame is None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if found:
                frame = auto_frame(clean, data_qubits)
            exact = verify(circuit_text, data_qubits, frame=frame, postselect=post, n_qubits=n_qubits)
    except (ValueError, ExactError) as e:
        report.correctness_errors.append(str(e))
        return report
    report.frame, report.ideal_acceptance = frame, exact.get("acceptance")
    report.correctness_errors = exact["errors"]
    if found and any("does not prepare" in e for e in exact["errors"]):
        report.correctness_errors.append(
            "(no FRAME_X/FRAME_Z declared: the judge used the Pauli frame it found itself, which can miss "
            "one; if the output is |T_L> up to a logical Pauli fixed by measurement results, declare it)")
    if not report.correct_ok:
        return report

    circuit_text, post = conditional(circuit_text, post, exact["conditional_detectors"])
    report.conditional_detectors = [k for k, _, _ in exact["conditional_detectors"]]
    noiseless = apply_noise(circuit_text, 0.0, noise)
    plan = Plan(data_qubits, frame, post, n_meas)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Cross-check of the two engines: xtim's noiseless runs must agree with the exact analysis.
        ideal = _sample(xtim.Circuit(noiseless + "\n" + plan.tail), 64, seed, plan)
        if ideal["acceptance"] == 0 or ideal["corrected"] or ideal["score"] > 1e-12:
            report.correctness_errors.append(
                "internal: noiseless runs disagree with the exact analysis; please report this circuit")
            return report

        # Gate 4 and the score.
        result = _sample(xtim.Circuit(apply_noise(circuit_text, p, noise) + "\n" + plan.tail), shots, seed, plan)
    report.score, report.sem, report.acceptance = result["score"], result["sem"], result["acceptance"]
    return report


def conditional(circuit_text, postselect, detectors):
    """Read conditional detectors (tprep.exact: deterministic only among POSTSELECT-ed runs) as
    POSTSELECT conditions on their accepted value. Returns the circuit with those DETECTOR lines
    turned into TICKs (line numbers and the other detectors' rec[-k] stay put; xtim's reference
    for a noiselessly random detector means nothing) and the extended POSTSELECT list."""
    if not detectors:
        return circuit_text, postselect
    drop, k, lines = {d for d, _, _ in detectors}, 0, circuit_text.splitlines()
    for i, line in enumerate(lines):
        if line.split("#")[0].split()[:1] == ["DETECTOR"]:
            if k in drop:
                lines[i] = "TICK"
            k += 1
    return "\n".join(lines) + "\n", list(postselect) + [(v, recs) for _, v, recs in detectors]


def auto_frame(clean, data_qubits):
    """The logical frame xtim finds when none is declared: records whose XOR fixes the signs of
    X^7 and Y^7, as {"X": (0, records), "Z": (0, records)}. `clean` carries its noise lines at
    strength 0 (see check_stim). tprep.exact then verifies it.

    Raises ValueError when xtim finds no such records (the output is then not |T_L>, or its frame
    needs a declaration).
    """
    import xtim

    n = clean.num_measurements
    signs = {}
    for P in TARGET:
        op = "*".join(f"{P}{q}" for q in data_qubits)
        try:
            signs[P] = {n - k for k in xtim.extract_frame(clean, op)}
        except xtim.XtimError:
            # |T_L> has a definite <X_L> and <Y_L> given the measurement record; a state where one
            # of them is random (e.g. a stabilizer state, or the wrong data_qubits) is not |T_L>.
            raise ValueError(f"the output on data_qubits is not |T_L>: the expectation of {P}^7 on it "
                             "is not fixed by the circuit (wrong state or wrong data_qubits)") from None
    # sign(X^7) flips with f_Z, sign(Y^7) with f_X XOR f_Z.
    return {"X": (0, sorted(signs["X"] ^ signs["Y"])), "Z": (0, sorted(signs["X"]))}


class Plan:
    """What the judge appends and how it reads the runs: the final round, X^7 and Y^7 with the
    frame's records, the frame's constants, and the POSTSELECT conditions (absolute records)."""

    def __init__(self, data_qubits, frame, postselect, n_meas):
        n_tail = 2 * len(SUPPORTS)
        (cx, sx), (cz, sz) = frame["X"], frame["Z"]
        recs = {"X": set(sz), "Y": set(sx) ^ set(sz)}     # X^7 flips with f_Z, Y^7 with f_X ^ f_Z
        self.frame = frame
        self.signs = np.array([(-1) ** cz, (-1) ** (cx ^ cz)], dtype=float)
        self.postselect = [(v, list(r)) for v, r in postselect]
        self.tail = final_round(data_qubits)
        for k, P in enumerate(TARGET):
            op = "*".join(f"{P}{q}" for q in data_qubits)
            self.tail += f"\nPAULI_EXPECTATION({k}) {op} " + " ".join(
                f"rec[-{n_meas + n_tail - r}]" for r in sorted(recs[P]))


def judge_tail(clean, data_qubits):
    """The final round and X-bar, Y-bar with the frame xtim finds (no declarations). Kept for tools
    written against judge 2.x; check() builds the tail with Plan."""
    return Plan(data_qubits, auto_frame(clean, data_qubits), (), clean.num_measurements).tail


def _sample(circuit, shots, seed, plan=None):
    """1 - F over the accepted shots, after a lookup-table correction from the final round.

    A shot is accepted when every DETECTOR reads 0 and every POSTSELECT parity of `plan` has its
    value; the final round never rejects. The frame's constants flip the signs of X^7 / Y^7. Any
    single-qubit Z correction flips X-bar and Y-bar, any single-qubit X correction flips Y-bar,
    so the correction is applied to the signs of the per-shot expectations.
    """
    n_syn = 2 * len(SUPPORTS)
    n_meas = circuit.num_measurements
    beta = np.array(list(TARGET.values()))
    kept = corrected = 0
    total = total_sq = 0.0
    for k, start in enumerate(range(0, shots, CHUNK)):
        # Bit-packed records: a circuit with thousands of measurements would otherwise take a byte
        # per record per shot (5 GB per million shots at the 5,000-statement limit).
        dets, exps, meas = circuit.compile_detector_sampler(seed=seed + k, engine="exact").sample(
            min(CHUNK, shots - start), return_expectations=True, return_measurements=True, bit_packed=True)
        keep = ~dets.any(axis=1)
        for v, recs in (plan.postselect if plan else ()):
            parity = np.zeros(len(keep), dtype=np.uint8)
            for r in recs:
                parity ^= (meas[:, r >> 3] >> (r & 7)) & 1
            keep &= parity == v
        syn = np.unpackbits(meas[keep], axis=1, bitorder="little")[:, n_meas - n_syn:n_meas].astype(bool)
        z_fix, x_fix = syn[:, :n_syn // 2].any(axis=1), syn[:, n_syn // 2:].any(axis=1)
        sign = np.stack([np.where(z_fix, -1, 1), np.where(z_fix ^ x_fix, -1, 1)], axis=1)
        if plan is not None:
            sign = sign * plan.signs
        infid = (1 - (exps[keep] * sign) @ beta) / 2
        kept += int(keep.sum())
        corrected += int((z_fix | x_fix).sum())
        total += infid.sum()
        total_sq += (infid**2).sum()
    if not kept:
        return {"score": float("nan"), "sem": float("nan"), "acceptance": 0.0, "corrected": corrected}
    mean = max(total / kept, 0.0)  # rounding leaves -1e-16 on a perfect run
    var = max(total_sq / kept - mean**2, 0.0)
    return {"score": mean, "sem": (var / kept) ** 0.5, "acceptance": kept / shots, "corrected": corrected}
