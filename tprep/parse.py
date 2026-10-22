"""Line-by-line reader for the Stim-superset text a submission is written in."""

import re

from .grid import ALLOWED, ARITY, NOISE_OPS, NON_UNITARY, N_QUBITS, adjacent

HEAD = re.compile(r"^\s*([A-Z_][A-Z_0-9]*)\s*(\([^)]*\))?\s*(.*?)\s*$")
DIRECTIVE = re.compile(r"^#!(FRAME_X|FRAME_Z|POSTSELECT)\(([01])\)((?:\s+\d+)*)\s*$")


class Instruction:
    """One parsed line. `qubits` holds plain qubit targets only (not rec[] or Pauli terms)."""

    __slots__ = ("line_no", "name", "arg", "targets", "qubits", "text")

    def __init__(self, line_no, name, arg, targets, text):
        self.line_no = line_no
        self.name = name
        self.arg = arg
        self.targets = targets
        self.text = text
        self.qubits = [int(t) for t in targets if t.isdigit()]

    def __repr__(self):
        return f"line {self.line_no}: {self.text}"


def parse(text):
    """Every non-blank, non-comment line of `text` as an Instruction."""
    out = []
    for i, raw in enumerate(text.splitlines(), start=1):
        line = raw.split("#")[0].strip()
        if not line:
            continue
        m = HEAD.match(line)
        if not m:
            raise ValueError(f"line {i}: cannot parse {raw!r}")
        name, arg, rest = m.group(1), m.group(2), m.group(3)
        out.append(Instruction(i, name, arg, rest.split(), line))
    return out


def check_format(text, n_qubits=N_QUBITS):
    """Validity gate 1. Returns a list of human-readable violations (empty if clean)."""
    errors = []
    try:
        program = parse(text)
    except ValueError as e:
        return [str(e)]

    if not program:
        return ["the circuit is empty"]

    for ins in program:
        if ins.name in NOISE_OPS:
            errors.append(f"line {ins.line_no}: {ins.name} is not allowed; the judge inserts noise")
            continue
        if ins.name not in ALLOWED:
            errors.append(f"line {ins.line_no}: {ins.name} is not in the allowed instruction set")
            continue

        for q in ins.qubits:
            if not 0 <= q < n_qubits:
                errors.append(f"line {ins.line_no}: qubit {q} is outside the chip (0..{n_qubits - 1})")

        arity = ARITY.get(ins.name)
        if ins.name in NON_UNITARY and len(ins.targets) != len(ins.qubits):
            # `M !5` and the like: the judge charges noise to plain qubit targets only.
            errors.append(f"line {ins.line_no}: {ins.name} takes plain qubit targets only")
            continue
        if arity is None:
            continue
        if len(ins.targets) != len(ins.qubits):
            errors.append(f"line {ins.line_no}: {ins.name} takes plain qubit targets only")
            continue
        if arity > 1 and len(ins.qubits) % arity:
            errors.append(
                f"line {ins.line_no}: {ins.name} needs targets in groups of {arity}, got {len(ins.qubits)}")
            continue
        if arity > 1:
            for i in range(0, len(ins.qubits), arity):
                group = ins.qubits[i:i + arity]
                for a, b in zip(group, group[1:]):
                    if not adjacent(a, b):
                        errors.append(
                            f"line {ins.line_no}: {ins.name} {a} {b} — qubits are not neighbours on the chip")
    return errors


def count_multiqubit_gates(text):
    """Two-qubit-and-larger gate count, the first tie-break in problem.md section 2.3."""
    n = 0
    for ins in parse(text):
        arity = ARITY.get(ins.name, 1)
        if arity > 1 and ins.qubits:
            n += len(ins.qubits) // arity
    return n


def declarations(text):
    """(frame, postselect) declared by `#!` lines (qasm.convert writes them from `// FRAME_X`, ...).

    frame is {"X": (c, records), "Z": (c, records)} or None when not declared; postselect is a list
    of (v, records). Records are absolute indices (0 = first measurement). Raises ValueError when
    only one of FRAME_X and FRAME_Z is given, or one is given twice.
    """
    frame, post = {}, []
    for line in text.splitlines():
        m = DIRECTIVE.match(line.strip())
        if not m:
            continue
        kind, c, recs = m.group(1), int(m.group(2)), [int(r) for r in m.group(3).split()]
        if kind == "POSTSELECT":
            post.append((c, recs))
        elif kind[-1] in frame:
            raise ValueError(f"{kind} is declared twice")
        else:
            frame[kind[-1]] = (c, recs)
    if len(frame) == 1:
        raise ValueError("declare FRAME_X and FRAME_Z exactly once each, or neither")
    return (frame or None), post


def count_measurements(text):
    return sum(len(i.qubits) for i in parse(text) if i.name in ("M", "MZ", "MX", "MY"))
