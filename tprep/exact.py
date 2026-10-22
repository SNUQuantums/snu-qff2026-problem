"""Exact, shot-free check that the accepted runs of a circuit leave |T_L> (after the logical frame).

The noiseless circuit is rewritten as one unitary on fresh wires (deferred measurement): `M q`
copies q onto a new wire with CX and the record is that wire's Z; `R q` moves q onto a new wire
in |0>. Every non-Clifford gate is a product of rotations exp(i theta P) with theta = +-pi/8 (T,
CS, CH). Pushing each rotation back through the Cliffords before it gives

    |Psi> = C R_k ... R_1 |0...0>,   R_j = exp(i theta_j P_j),

with C the whole circuit's Clifford part (a stim tableau) and P_j the rotation axes moved to the
start. |phi> = R_k ... R_1 |0...0> is kept as a sparse vector over computational basis states;
its size is 2^(rank of the axes' X parts), so the method is exact whenever that rank is small.
Any Pauli expectation is then <phi| C^dag Q C |phi>, and every quantity the judge needs is a sum
of Pauli expectations:

    acceptance  A = <Pi_acc>                       Pi_acc = prod_k (1 + (-1)^v_k Z_{S_k}) / 2
    fidelity    F = <Pi_acc  P_frame^dag |T_L><T_L| P_frame>

where Z_S is the product of Z over the record wires in S (a parity of records), and
|T_L><T_L| = Pi_code (1 + (X^7 - Y^7) / sqrt 2) / 2 on the data qubits. The accepted runs all
leave |T_L> exactly when F = A (> 0). Floating-point arithmetic only: "exact" means no sampling
and no approximation beyond rounding; the check is (A - F) / A <= 1e-9, per accepted run.
"""

import math
from itertools import product

import numpy as np

from .grid import N_QUBITS, SUPPORTS
from .parse import parse

TOL = 1e-9          # on quantities normalized by the acceptance (per accepted run)
A_ZERO = 1e-12      # an acceptance this small is rounding of 0
MAX_RANK = 22   # 2^22 amplitudes; above this the check refuses instead of running for hours

# Non-Clifford gates as rotations exp(i theta P), in time order, P on (first, second) target.
_P8 = math.pi / 8
ROTATIONS = {
    "T": [("Z_", -_P8)], "T_DAG": [("Z_", _P8)],
    # CS = diag(1,1,1,i) = exp(i pi/8 (1 - Z_a - Z_b + Z_a Z_b))
    "CS": [("Z_", -_P8), ("_Z", -_P8), ("ZZ", _P8)],
    "CS_DAG": [("Z_", _P8), ("_Z", _P8), ("ZZ", -_P8)],
}
CLIFFORD = {"I", "X", "Y", "Z", "H", "S", "S_DAG", "SQRT_X", "SQRT_X_DAG", "CX", "CZ", "SWAP"}


class ExactError(Exception):
    """The circuit is outside what the exact check can handle (reported, never guessed)."""


class _State:
    def __init__(self, n_wires):
        import stim

        self.stim = stim
        self.n = n_wires
        self.sim = stim.TableauSimulator()
        self.sim.set_num_qubits(n_wires)
        self.phi = {0: 1 + 0j}
        self.span = []          # GF(2) basis of the X parts seen so far (for the rank guard)

    def _masks(self, ps):
        xs, zs = ps.to_numpy()
        x = int.from_bytes(np.packbits(xs, bitorder="little").tobytes(), "little")
        z = int.from_bytes(np.packbits(zs, bitorder="little").tobytes(), "little")
        return x, z, complex(ps.sign) * 1j ** (x & z).bit_count()

    def _add_to_span(self, x):
        for v in self.span:
            x = min(x, x ^ v)
        if x:
            self.span.append(x)
            if len(self.span) > MAX_RANK:
                raise ExactError(f"the circuit's magic and postselection span more than 2^{MAX_RANK} "
                                 "branches, beyond the exact check")

    def rotate(self, paulis, theta):
        """Apply exp(i theta P) now; P maps qubit -> 'X'/'Y'/'Z' on current wires."""
        ps = self.stim.PauliString(self.n)
        for q, c in paulis.items():
            ps[q] = c
        x, z, coef = self._masks(self.sim.current_inverse_tableau()(ps))
        self._add_to_span(x)
        c, s = math.cos(theta), 1j * math.sin(theta) * coef
        new = {}
        for b, a in self.phi.items():
            new[b] = new.get(b, 0) + c * a
            nb = b ^ x
            new[nb] = new.get(nb, 0) + s * (-1) ** (b & z).bit_count() * a
        self.phi = {b: a for b, a in new.items() if abs(a) > 1e-15}

    def apply(self, vec, ps, inverse):
        """ps|vec> for a stim PauliString on the wires, |Psi> = C|vec>."""
        x, z, coef = self._masks(inverse(ps))
        self._add_to_span(x)
        out = {}
        for b, a in vec.items():
            out[b ^ x] = out.get(b ^ x, 0) + coef * (-1) ** (b & z).bit_count() * a
        return out

    def expectation(self, ps, inverse, vec=None):
        """<vec| ps |vec> (vec defaults to the circuit's state)."""
        vec = self.phi if vec is None else vec
        x, z, coef = self._masks(inverse(ps))
        total = 0j
        for b, a in vec.items():
            a2 = vec.get(b ^ x)
            if a2 is not None:
                total += a2.conjugate() * coef * (-1) ** (b & z).bit_count() * a
        return total.real


def run(stim_text, n_qubits):
    """Build |Psi>. Returns (state, final inverse tableau, records: list of wires, detectors:
    list of record-index lists, current wire of each physical qubit)."""
    prog = parse(stim_text)
    fresh = sum(len(i.qubits) for i in prog if i.name in ("M", "MX", "R", "RX"))
    st = _State(n_qubits + fresh)
    wire = list(range(n_qubits))
    nxt = n_qubits
    records, detectors = [], []
    for ins in prog:
        name, qs = ins.name, ins.qubits
        if name == "DETECTOR":
            detectors.append([len(records) - int(t[5:-1]) for t in ins.targets])
            continue
        if name in ("TICK", "QUBIT_COORDS", "SHIFT_COORDS") or not qs and name not in ROTATIONS:
            continue
        if name in ("M", "MX"):
            for q in qs:
                w, nxt = nxt, nxt + 1
                if name == "MX":
                    st.sim.h(wire[q])
                st.sim.cx(wire[q], w)
                if name == "MX":
                    st.sim.h(wire[q])
                records.append(w)
        elif name in ("R", "RX"):
            for q in qs:
                wire[q], nxt = nxt, nxt + 1
                if name == "RX":
                    st.sim.h(wire[q])
        elif name in CLIFFORD:
            if name != "I":
                st.sim.do(st.stim.CircuitInstruction(name, [wire[q] for q in qs]))
        elif name in ROTATIONS:
            arity = 2 if name.startswith("CS") else 1
            for g in range(0, len(qs), arity):
                tgt = [wire[q] for q in qs[g:g + arity]]
                for axis, theta in ROTATIONS[name]:
                    st.rotate({t: c for t, c in zip(tgt, axis) if c != "_"}, theta)
        elif name == "CH":
            # CH = exp(-i pi/8 Y_b) CZ exp(i pi/8 Y_b): the target's H is Y-rotation conjugate of Z.
            for g in range(0, len(qs), 2):
                a, b = wire[qs[g]], wire[qs[g + 1]]
                st.rotate({b: "Y"}, _P8)
                st.sim.cz(a, b)
                st.rotate({b: "Y"}, -_P8)
        else:
            raise ExactError(f"line {ins.line_no}: {name} is not supported by the exact check")
    return st, st.sim.current_inverse_tableau(), records, detectors, wire


def _parity(st, wires_z, extra=None):
    ps = st.stim.PauliString(st.n)
    for w in wires_z:
        ps[w] = "Z" if ps[w] == 0 else "I"   # parity: a repeated record cancels
    if extra:
        ps = ps * extra
    return ps


def verify(stim_text, data_qubits, *, frame=None, postselect=(), n_qubits=N_QUBITS):
    """Exact acceptance and correctness of the accepted noiseless runs.

    frame: {"X": (const, [record indices]), "Z": (const, [record indices])} — the logical
    correction X_L^f_X Z_L^f_Z with f = const XOR records. postselect: [(value, [record indices])].
    Record indices are absolute (0 = first measurement). Returns a dict with acceptance,
    infidelity (of the accepted, corrected output, before the judge's syndrome round),
    detector_values (exact <(-1)^parity> per DETECTOR among the accepted runs, each must be +-1),
    conditional_detectors and errors. POSTSELECT compares the raw parity with v; a DETECTOR compares
    with its value in the accepted noiseless runs.
    """
    st, inv, records, detectors, wire = run(stim_text, n_qubits)
    E = lambda ps: st.expectation(ps, inv)
    rec_wires = lambda idx: [records[i] for i in idx]
    out = {"errors": []}

    data = [wire[q] for q in data_qubits]
    frame = frame or {"X": (0, []), "Z": (0, [])}
    cx, sx = frame["X"]
    cz, sz = frame["Z"]

    def logical(pauli):
        ps = st.stim.PauliString(st.n)
        for w in data:
            ps[w] = pauli
        return ps

    stab_gens = []
    for P in "XZ":
        for s in SUPPORTS:
            g = st.stim.PauliString(st.n)
            for i in s:
                g[data[i]] = P
            stab_gens.append(g)
    group = []
    for bits in product((0, 1), repeat=len(stab_gens)):
        g = st.stim.PauliString(st.n)
        for on, s in zip(bits, stab_gens):
            if on:
                g = g * s
        group.append(g)

    # The accepted part of the state: each POSTSELECT is the projector (1 + (-1)^v Z_S) / 2. They
    # commute with each other and with every operator below, so <Pi O Pi> = <acc| O |acc>.
    acc = st.phi
    for v, idx in postselect:
        flipped = st.apply(acc, _parity(st, rec_wires(idx)), inv)
        sign = (-1) ** v
        acc = {b: (acc.get(b, 0) + sign * flipped.get(b, 0)) / 2 for b in set(acc) | set(flipped)}
        acc = {b: a for b, a in acc.items() if abs(a) > 1e-15}
    A = sum(abs(a) ** 2 for a in acc.values())
    out["acceptance"] = A
    if A <= A_ZERO:
        out["errors"].append("no noiseless run is accepted")
        return out

    # A DETECTOR compares its parity with the value it has in every accepted noiseless run, so it
    # must be deterministic there; a parity that is always 1 is fine. Deterministic without
    # POSTSELECT is the usual case (xtim fires it against its noiseless reference). Deterministic only
    # among the POSTSELECT-ed runs makes it "conditional": the judge then reads it like a POSTSELECT
    # of that value (out["conditional_detectors"] = [(detector index, value, records)]).
    parities = [_parity(st, rec_wires(d)) for d in detectors]
    out["detector_values"] = [st.expectation(P, inv, acc) / A for P in parities]
    random = [k for k, v in enumerate(out["detector_values"]) if abs(abs(v) - 1) > TOL]
    if random:
        out["errors"].append("detector(s) " + ", ".join(f"#{k + 1}" for k in random) + " (in the order "
                             "written) are not deterministic: noiselessly they read 1 for some measurement "
                             "outcomes" + (" among the runs POSTSELECT keeps" if postselect else "") +
                             " (use POSTSELECT for a condition that is random by design)")
        return out
    out["conditional_detectors"] = [
        (k, 0 if v > 0 else 1, detectors[k]) for k, (P, v) in enumerate(zip(parities, out["detector_values"]))
        if postselect and abs(abs(E(P)) - 1) > TOL]
    # Signed parities: (-1)^(c + records) as (sign, PauliString on record wires).
    sZ = ((-1) ** cz, _parity(st, rec_wires(sz)))
    sXZ = ((-1) ** (cx ^ cz), _parity(st, rec_wires(list(sx) + list(sz))))
    # Corrected target: P^dag (X^7 - Y^7) P = s_Z X^7 - s_XZ Y^7 (X_L flips Y-bar, Z_L flips both).
    X7, Y7 = logical("X"), logical("Y")
    target = [(1.0, st.stim.PauliString(st.n)),
              (sZ[0] / math.sqrt(2), X7 * sZ[1]),
              (-sXZ[0] / math.sqrt(2), Y7 * sXZ[1])]

    code = sum(st.expectation(g, inv, acc) for g in group) / len(group)          # <Pi_acc Pi_code>
    F = sum(t * st.expectation(g * T, inv, acc) for g in group for t, T in target) / (2 * len(group))
    out["infidelity"] = (A - F) / A
    # Compared per accepted run, (A - x) / A: an absolute A - x would let a tiny acceptance hide
    # a wrong output.
    if (A - code) / A > TOL:
        out["errors"].append("noiselessly the output is not in the code space in "
                             f"{(A - code) / A:.6g} of the accepted runs: a stabilizer of the output is -1 "
                             "or depends on the measurement outcomes")
    elif out["infidelity"] > TOL:
        out["errors"].append("the circuit does not prepare |T_L>: after the logical frame, the accepted "
                             f"noiseless runs have infidelity {out['infidelity']:.6g}")
    out["rank"] = len(st.span)
    return out
