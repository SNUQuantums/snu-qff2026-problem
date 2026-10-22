"""The chip: qubit numbering, adjacency, and the instruction set a submission may use."""

ROWS, COLS = 6, 6
N_QUBITS = ROWS * COLS

# Target state: Steane [[7,1,3]], code qubit i <-> binary i+1.
SUPPORTS = ({0, 2, 4, 6}, {1, 2, 5, 6}, {3, 4, 5, 6})
N_CODE = 7

ACCEPTANCE_FLOOR = 0.5

ONE_QUBIT = frozenset({"H", "S", "S_DAG", "SQRT_X", "SQRT_X_DAG", "X", "Y", "Z", "I", "T", "T_DAG"})
TWO_QUBIT = frozenset({"CX", "CNOT", "CZ", "SWAP", "CS", "CS_DAG", "CH"})
# No MPP: a Pauli-product measurement would reach across the chip in one noiseless step.
NON_UNITARY = frozenset({"M", "MZ", "MX", "MY", "R", "RZ", "RX"})
ANNOTATION = frozenset({"DETECTOR", "TICK", "QUBIT_COORDS", "SHIFT_COORDS"})

ALLOWED = ONE_QUBIT | TWO_QUBIT | NON_UNITARY | ANNOTATION

# Participants must not bake noise in; the judge inserts it (problem.md section 1.5).
NOISE_OPS = frozenset({
    "DEPOLARIZE1", "DEPOLARIZE2", "X_ERROR", "Y_ERROR", "Z_ERROR",
    "PAULI_CHANNEL_1", "PAULI_CHANNEL_2", "HERALDED_ERASE", "E", "ELSE_CORRELATED_ERROR",
})

ARITY = {}
ARITY.update({op: 1 for op in ONE_QUBIT})
ARITY.update({op: 2 for op in TWO_QUBIT})


def position(q):
    return divmod(q, COLS)


def adjacent(a, b):
    (ra, ca), (rb, cb) = position(a), position(b)
    return abs(ra - rb) + abs(ca - cb) == 1


def couplers():
    """Every nearest-neighbour pair, as sorted tuples."""
    out = []
    for q in range(N_QUBITS):
        r, c = position(q)
        if c + 1 < COLS:
            out.append((q, q + 1))
        if r + 1 < ROWS:
            out.append((q, q + COLS))
    return out


def grid_diagram():
    """The chip as ASCII, the way problem.md prints it."""
    lines = []
    for r in range(ROWS):
        lines.append(" -- ".join(f"{COLS * r + c:2d}" for c in range(COLS)))
        if r + 1 < ROWS:
            lines.append("  |    " * COLS)
    return "\n".join(lines)
