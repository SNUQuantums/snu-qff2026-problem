"""The submission format: OpenQASM 2.0 from Qiskit, plus `// DETECTOR` lines (problem.md section 1.3).

The judge turns a submission into the Stim text xtim reads. Each `measure` becomes one measurement
record, in program order, and each `// DETECTOR c[i] c[j] ...` line becomes a Stim DETECTOR over the
latest measurement written to each named bit.

Optional comment lines declare more, each over the latest measurement written to the named bits:
  // POSTSELECT(v) c[i] ...   keep only runs where the bits' XOR is v (may be random noiselessly)
  // FRAME_X(c) c[i] ...      f_X = c XOR the bits; likewise FRAME_Z. The output is read after the
  // FRAME_Z(c) c[j] ...      logical correction X_L^f_X Z_L^f_Z. Give both or neither; with
                              neither the judge finds the Pauli frame itself.
They reach the judge's Stim text as `#!FRAME_X(c) r ...` lines at the end, r the absolute record
index (parse.declarations reads them).
"""

import re

from .grid import N_QUBITS, adjacent
from .noise import STANDARD

# Qiskit name -> Stim name. These are the only operations a submission may use.
GATES = {
    "id": "I", "x": "X", "y": "Y", "z": "Z", "h": "H", "s": "S", "sdg": "S_DAG",
    "sx": "SQRT_X", "sxdg": "SQRT_X_DAG", "t": "T", "tdg": "T_DAG",
    "cx": "CX", "cz": "CZ", "swap": "SWAP", "ch": "CH", "cs": "CS", "csdg": "CS_DAG",
    "measure": "M", "reset": "R",
}
MAX_BYTES = 200_000       # UTF-8 bytes of the submission text (problem.md section 3)
MAX_STATEMENTS = 5_000    # QASM statements, i.e. `;` outside `//` comments
# measure + reset. xtim gives every one a fresh wire, and its analysis grows faster than the square
# of the wire count: 1,000 take about 20 s, 2,000 resets over 100 s, 4,900 measurements >10 min.
MAX_MEASURE_RESET = 1_000
DETECTOR = re.compile(r"^\s*//\s*DETECTOR\b(.*)$", re.M)
DECLARATION = re.compile(r"^\s*//\s*(FRAME_X|FRAME_Z|POSTSELECT)\b(.*)$", re.M)
CONSTANT = re.compile(r"^\s*\(\s*([01])\s*\)(.*)$")
BIT = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)\[(\d+)\]$")


def _custom_instructions():
    from qiskit import qasm2
    from qiskit.circuit.library import CSdgGate, CSGate

    # qelib1.inc has no cs/csdg; qasm2.dumps writes a definition for them, which this replaces.
    return list(qasm2.LEGACY_CUSTOM_INSTRUCTIONS) + [
        qasm2.CustomInstruction("cs", 0, 2, CSGate, builtin=True),
        qasm2.CustomInstruction("csdg", 0, 2, CSdgGate, builtin=True),
    ]


def load(text):
    """Parse submission text into a QuantumCircuit. `include` cannot reach the judge's disk."""
    from qiskit import qasm2

    return qasm2.loads(text, include_path=(), custom_instructions=_custom_instructions())


def _bit_names(qc, bits):
    names = []
    for b in bits:
        if isinstance(b, tuple):
            names.append(f"{b[0]}[{b[1]}]")
        else:
            reg, idx = qc.find_bit(b).registers[0]
            names.append(f"{reg.name}[{idx}]")
    return names


def detector_lines(qc, detectors):
    """`// DETECTOR` lines for `detectors`, each a list of Clbits or (register name, index) pairs."""
    return ["// DETECTOR " + " ".join(_bit_names(qc, det)) for det in detectors]


def declaration_lines(qc, frame=None, postselect=()):
    """`// FRAME_X`, `// FRAME_Z` and `// POSTSELECT` lines. frame = {"X": (c, bits), "Z": (c, bits)},
    postselect = [(v, bits), ...]; bits as in detector_lines."""
    lines = [f"// POSTSELECT({v}) " + " ".join(_bit_names(qc, bits)) for v, bits in postselect]
    if frame is not None:
        lines += [f"// FRAME_{P}({frame[P][0]}) " + " ".join(_bit_names(qc, frame[P][1])) for P in "XZ"]
    return [line.rstrip() for line in lines]


def submission_text(qc, detectors=(), frame=None, postselect=()):
    """The text to submit: `qasm2.dumps(qc)`, then the `// DETECTOR` lines, then any `// POSTSELECT`
    and `// FRAME_X`/`// FRAME_Z` lines (see declaration_lines)."""
    from qiskit import qasm2

    return "\n".join([qasm2.dumps(qc).rstrip(), *detector_lines(qc, detectors),
                      *declaration_lines(qc, frame, postselect)]) + "\n"


def _records(tokens, by_name):
    """Absolute record indices of `c[i]` tokens, or (None, the first bad token)."""
    recs = []
    for tok in tokens:
        b = BIT.match(tok)
        if not b or (b.group(1), int(b.group(2))) not in by_name:
            return None, tok
        recs.append(by_name[(b.group(1), int(b.group(2)))])
    return recs, None


def _x_readouts(lines, drop_h=False):
    """Run `h q` then `measure q` as `I q` then `MX q` when q is not used again before a reset, and
    `h q; measure q; h q` as `I q; MX q; I q`. With `drop_h` (native X-basis readout) the h's are
    removed instead: their lines become None.

    Both are the same channel: the h's depolarizing noise lands on the I, and the measurement flip
    lands as Z_ERROR before MX instead of X_ERROR after the h. Only the qubit left behind differs
    (|+>/|-> instead of |0>/|1>), so the rewrite is made only when nothing reads it afterwards.
    xtim samples circuits with plain H + M readouts that feed the Pauli frame about 20 times
    slower than with MX, with the same result.
    """
    ops = [(i, line.split()[0], [int(q) for q in line.split()[1:]]) for i, line in enumerate(lines)]
    on = {}
    for i, name, qs in ops:
        for q in qs:
            on.setdefault(q, []).append((i, name))
    for seq in on.values():
        for k in range(1, len(seq)):
            (j, prev), (i, name) = seq[k - 1], seq[k]
            after = seq[k + 1][1] if k + 1 < len(seq) else None
            if name != "M" or prev != "H" or lines[j] is None or not lines[j].startswith("H "):
                continue
            q = lines[i].split()[1]
            if after in (None, "R"):
                lines[j], lines[i] = (None if drop_h else f"I {q}"), f"MX {q}"
            elif after == "H":
                # h, measure, h is an X readout that leaves |+>/|-> for the next use: exactly MX
                # between two identities, noise included (conjugate the middle by H).
                h2 = seq[k + 1][0]
                lines[j], lines[i], lines[h2] = (None if drop_h else f"I {q}"), f"MX {q}", (None if drop_h else f"I {q}")


def _x_preparations(lines):
    """Run `reset q` then `h q` as one `RX q`: the same state, prepared natively (lines[j] = None)."""
    nxt = {}
    for i in range(len(lines) - 1, -1, -1):
        if lines[i] is None:
            continue
        name, *qs = lines[i].split()
        if name == "R" and nxt.get(qs[0], (None, ""))[1] == "H":
            j = nxt[qs[0]][0]
            lines[i], lines[j] = f"RX {qs[0]}", None
        for q in qs:
            nxt[q] = (i, name)


def to_stim(text, n_qubits=N_QUBITS, noise=STANDARD):
    """(stim text, format errors). The stim text is None when there are errors."""
    stim, _, errors = convert(text, n_qubits, noise)
    return stim, errors


def convert(text, n_qubits=N_QUBITS, noise=STANDARD):
    """(stim text, labels, format errors): labels[i] names the submission instruction behind
    line i of the stim text, for messages and drawings. `noise` may narrow the allowed operations
    and make X-basis readout and preparation native (NoiseModel.gates, .native_x)."""
    errors = []
    size = len(text.encode("utf-8"))
    statements = sum(line.split("//")[0].count(";") for line in text.splitlines())
    if size > MAX_BYTES:
        return None, [], [f"the submission is {size} bytes; the limit is {MAX_BYTES}"]
    if statements > MAX_STATEMENTS:
        return None, [], [f"the submission has {statements} statements; the limit is {MAX_STATEMENTS}"]
    first = next((l.strip() for l in text.splitlines() if l.strip() and not l.strip().startswith("//")), "")
    if not first.startswith("OPENQASM 2"):
        return None, [], ["submit OpenQASM 2.0 text starting with `OPENQASM 2.0;` (problem.md section 1.3)"]
    try:
        qc = load(text)
    except Exception as e:  # qiskit.qasm2.QASM2ParseError and friends
        return None, [], [f"OpenQASM 2.0 could not be parsed: {str(e)[:300]}"]

    if len(qc.qregs) != 1:
        return None, [], ["declare exactly one quantum register; qubit q[i] is chip qubit i"]
    if qc.num_qubits > n_qubits:
        errors.append(f"the quantum register has {qc.num_qubits} qubits; the chip has {n_qubits}")

    lines, labels, last_write, n_meas, meas_line = [], [], {}, 0, []
    barriers = {}  # stim line -> qubit groups of the barriers written just before it
    for k, inst in enumerate(qc.data, start=1):
        op = inst.operation
        qs = [qc.find_bit(q).index for q in inst.qubits]
        where = f"instruction {k} ({op.name} {', '.join(f'q[{q}]' for q in qs)})"
        if op.name == "barrier":
            barriers.setdefault(len(lines), []).append(" ".join(map(str, qs)))
            continue
        if op.name not in GATES:  # includes `if (c==1) ...`, which Qiskit loads as if_else
            errors.append(f"{where}: {op.name} is not an allowed operation (problem.md section 1.2)")
            continue
        if noise.gates is not None and op.name not in noise.gates:
            errors.append(f"{where}: {op.name} is not allowed in this noise model")
            continue
        if len(qs) == 2 and not adjacent(*qs):
            errors.append(f"{where}: q[{qs[0]}] and q[{qs[1]}] are not neighbours on the chip")
            continue
        lines.append(f"{GATES[op.name]} {' '.join(map(str, qs))}")
        labels.append(where)
        if op.name == "measure":
            last_write[qc.find_bit(inst.clbits[0]).registers[0]] = n_meas
            meas_line.append(len(lines) - 1)
            n_meas += 1

    _x_readouts(lines, drop_h=noise.native_x)
    if noise.native_x:
        _x_preparations(lines)
    n_mr = sum(1 for line in lines if line and line.split()[0] in ("M", "R", "MX", "RX"))
    if n_mr > MAX_MEASURE_RESET:
        errors.append(f"the circuit has {n_mr} measure and reset operations; the limit is {MAX_MEASURE_RESET}")
    by_name = {(reg.name, idx): rec for (reg, idx), rec in last_write.items()}
    after = {}  # stim line -> detectors declared right after it
    for i, m in enumerate(DETECTOR.finditer(text), start=1):
        recs, bad = _records(m.group(1).split(), by_name)
        if bad:
            errors.append(f"DETECTOR line {i}: {bad} is not a measured classical bit")
        else:
            if not recs:
                errors.append(f"DETECTOR line {i} names no bits")
            # Declared right after its last measurement, as Stim circuits do: xtim evaluates a
            # detector where it stands, and one placed after later non-Clifford gates can read as
            # random once the judge's final round is appended.
            last = max(recs, default=0)
            after.setdefault(meas_line[last] if recs else len(lines) - 1, []).append(
                ("DETECTOR " + " ".join(f"rec[-{last + 1 - r}]" for r in recs), f"DETECTOR line {i}"))
    directives = []
    seen = {"FRAME_X": 0, "FRAME_Z": 0}
    for m in DECLARATION.finditer(text):
        kind, rest = m.group(1), CONSTANT.match(m.group(2))
        where = f"{kind} line `{m.group(0).strip()}`"
        if not rest:
            errors.append(f"{where}: write {kind}(0) or {kind}(1) followed by classical bits")
            continue
        recs, bad = _records(rest.group(2).split(), by_name)
        if bad:
            errors.append(f"{where}: {bad} is not a measured classical bit")
        elif kind == "POSTSELECT" and not recs:
            errors.append(f"{where}: names no bits")
        else:
            seen[kind] = seen.get(kind, 0) + 1
            directives.append(f"#!{kind}({rest.group(1)}) " + " ".join(map(str, recs)))
    if seen["FRAME_X"] > 1 or seen["FRAME_Z"] > 1 or seen["FRAME_X"] != seen["FRAME_Z"]:
        errors.append("declare FRAME_X and FRAME_Z exactly once each, or neither "
                      "(then the judge finds the Pauli frame itself)")
    if not lines:
        errors.append("the circuit is empty")
    out, out_labels, held = [], [], []
    for k, (line, label) in enumerate(zip(lines, labels)):
        held += barriers.get(k, [])
        if line is None:  # an h folded into an X-basis readout or preparation
            continue
        # A barrier rides on the next instruction as a comment, so line numbers stay put; only
        # idle noise reads it (noise.barriers). A barrier after the last instruction does nothing.
        out.append(line + (f"  # barrier {' | '.join(held)}" if held else ""))
        held = []
        out_labels.append(label)
        for det, det_label in after.get(k, []):
            out.append(det)
            out_labels.append(det_label)
    for d in directives:
        out.append(d.rstrip())
        out_labels.append(d.split("(")[0][2:] + " line")
    return (None if errors else "\n".join(out) + "\n"), out_labels, errors
