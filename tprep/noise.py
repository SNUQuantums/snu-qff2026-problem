"""The judge's noise model (problem.md section 1.5) and the schedule its idle noise follows."""

from dataclasses import dataclass

from .grid import ARITY
from .parse import parse


@dataclass(frozen=True)
class NoiseModel:
    """Channel strengths as multiples of p. Set a field to 0 to switch that channel off. The
    defaults are judge 2.x's rules (LEGACY); the official model since 3.0 is STANDARD."""

    gate1: float = 1.0   # DEPOLARIZE1(p) after every 1-qubit gate
    gate2: float = 1.0   # DEPOLARIZE2(p) after every 2-qubit gate
    meas: float = 1.0    # flip of every measurement result (X_ERROR before M, Z_ERROR before MX)
    reset: float = 1.0   # the orthogonal state after every reset (X_ERROR after R, Z_ERROR after RX)
    idle: float = 0.0    # DEPOLARIZE1(p) per time step on every waiting qubit (see schedule)
    meas_after: float = 0.0   # DEPOLARIZE1(p) after every measurement
    initial: float = 0.0      # a qubit's first use, unless it is a reset, starts with a noisy |0>
    # X-basis readout and preparation are native: `h` then `measure`, and `reset` then `h`, run as
    # one MX / RX (one time step, no h noise). A readout is folded only when the qubit is reset or
    # unused afterwards (MX leaves |+>/|->, h + measure leaves |0>/|1>), or when an h follows it.
    native_x: bool = False
    gates: frozenset | None = None   # QASM operations allowed; None means every one in qasm.GATES


# Uniform depolarizing circuit noise: CNOT is the only two-qubit
# gate, waiting qubits depolarize, X and Z preparations and measurements are native, a measurement
# also depolarizes its qubit, and every qubit is prepared with noise.
STANDARD = NoiseModel(idle=1.0, meas_after=1.0, initial=1.0, native_x=True, gates=frozenset(
    {"id", "x", "y", "z", "h", "s", "sdg", "sx", "sxdg", "t", "tdg", "cx", "measure", "reset", "barrier"}))
UNIFORM = STANDARD
LEGACY = NoiseModel()             # judge 2.x: all of cx, cz, swap, cs, csdg, ch; no idle noise
IDLE = NoiseModel(idle=1.0)       # judge 2.4's comparison option: LEGACY plus idle noise

_MEAS = {"M", "MZ", "MX", "MY"}
_RESET = {"R", "RZ", "RX"}
# The Pauli that flips each basis: X_ERROR does nothing to an X-basis measurement or reset.
_FLIP = {"M": "X", "MZ": "X", "MY": "X", "R": "X", "RZ": "X", "MX": "Z", "RX": "Z"}


def _repeated(q, n, first_order):
    """DEPOLARIZE1(q) applied n times, as one channel (or its first-order weight n q)."""
    return n * q if first_order else 0.75 * (1 - (1 - 4 * q / 3) ** n)


def barriers(line):
    """Qubit groups of the `barrier`s written just before this Stim line (see qasm.convert)."""
    _, _, note = line.partition("# barrier")
    return [[int(q) for q in g.split()] for g in note.split("|")] if note else []


def schedule(text):
    """When each operation runs, and when each qubit waits.

    Operations run in the order written, each at the earliest time step its qubits are free (TICK
    is ignored). A `barrier` holds its qubits until all of them are free, as in Qiskit. A qubit
    waits between two of its operations, and from its last operation to the end of the circuit,
    except before a reset and after a measurement it is not used again (it no longer holds a
    state anyone reads). Returns (prog, steps, gaps): steps[i] lists (qubits, time step) of
    instruction i, one entry per gate it applies; gaps lists (i, qubit, first step, n): the qubit
    waits n steps from `first step`, after instruction i.
    """
    prog = parse(text)
    raw = text.splitlines()
    free, last, steps, gaps = {}, {}, [], []   # last: qubit -> (instruction, step, was a measurement)
    for i, ins in enumerate(prog):
        for group in barriers(raw[ins.line_no - 1]):
            step = max(free.get(q, 0) for q in group)
            free.update((q, step) for q in group)
        groups = []
        arity = ARITY.get(ins.name, 1)
        for g in range(0, len(ins.qubits), arity):
            group = ins.qubits[g:g + arity]
            step = max(free.get(q, 0) for q in group)
            for q in group:
                if q in last and step - last[q][1] > 1 and ins.name not in _RESET:
                    gaps.append((last[q][0], q, last[q][1] + 1, step - last[q][1] - 1))
                free[q] = step + 1
                last[q] = (i, step, ins.name in _MEAS)
            groups.append((group, step))
        steps.append(groups)
    end = max(free.values(), default=0)
    for q, (i, step, measured) in last.items():
        if not measured and end - step > 1:
            gaps.append((i, q, step + 1, end - step - 1))
    return prog, steps, sorted(gaps)


def noise_blocks(text, p, noise=STANDARD, first_order=False):
    """(instructions, blocks): blocks[i] is instruction i with its noise lines around it.

    A run of n waiting steps (see schedule) is one DEPOLARIZE1 equal to n repetitions, placed right
    after the operation that starts it, so the instruction order and every rec[-k] stay as written.
    `first_order` writes that channel's strength as n * idle * p instead, for the fault scan.
    """
    prog, _, gaps = schedule(text)
    blocks, seen = [], set()
    for ins in prog:
        qs = " ".join(str(q) for q in ins.qubits)
        arity = ARITY.get(ins.name)
        block = []
        fresh = [q for q in dict.fromkeys(ins.qubits) if q not in seen]
        seen.update(ins.qubits)
        if noise.initial and fresh and ins.name not in _RESET:
            block.append(f"X_ERROR({p * noise.initial}) " + " ".join(map(str, fresh)))
        if ins.name in _MEAS and noise.meas and ins.qubits:
            block.append(f"{_FLIP[ins.name]}_ERROR({p * noise.meas}) {qs}")
        block.append(ins.text)
        if ins.name in _MEAS and noise.meas_after and ins.qubits:
            block.append(f"DEPOLARIZE1({p * noise.meas_after}) {qs}")
        if ins.name in _RESET and noise.reset and ins.qubits:
            block.append(f"{_FLIP[ins.name]}_ERROR({p * noise.reset}) {qs}")
        elif arity == 1 and noise.gate1 and ins.qubits:
            block.append(f"DEPOLARIZE1({p * noise.gate1}) {qs}")
        elif arity == 2 and noise.gate2 and ins.qubits:
            block.append(f"DEPOLARIZE2({p * noise.gate2}) {qs}")
        blocks.append(block)
    if noise.idle:
        for i, q, _, n in gaps:
            blocks[i].append(f"DEPOLARIZE1({_repeated(p * noise.idle, n, first_order)}) {q}")
    return prog, blocks


def apply_noise(text, p, noise=STANDARD):
    """Return `text` with noise instructions inserted. Participant circuits carry none."""
    return "\n".join(line for block in noise_blocks(text, p, noise)[1] for line in block)
