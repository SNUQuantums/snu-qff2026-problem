"""One picture of a submission: where it sits on the chip, and its gates in time.

    from tprep import draw, fault_scan
    fig = draw(circuit, data_qubits)                          # chip + timeline
    fig = draw(circuit, data_qubits, scan=fault_scan(circuit, data_qubits))   # + first-order culprits
    fig.savefig("circuit.svg")   # vector: zoom in without blur

or from a terminal: `python -m tprep.vis circuit.qasm 1,2,0,13,6,7,12 --faults --out circuit.svg`.

Left, the 6x6 chip: data qubits are filled and labelled with their code-qubit number, and a
coupler's width grows with the two-qubit gates on it. Right, the circuit as the judge sees it,
each instruction placed at the time step the judge gives it. A small red dot marks each place the
judge inserts noise, including every time step a qubit waits. With a fault scan, the dots where a
single fault gets past every detector and breaks the output grow with their share of c1.
Labels are in English so that missing Korean fonts cannot break them.
"""

import sys

from .grid import COLS, ROWS, position
from .noise import STANDARD, schedule
from .qasm import convert, submission_text

TWO = {"CX", "CZ", "SWAP", "CH", "CS", "CS_DAG"}
MEASURES = {"M", "MX"}
RESETS = {"R", "RX"}
SHORT = {"S_DAG": "S†", "T_DAG": "T†", "SQRT_X": "√X", "SQRT_X_DAG": "√X†", "CS_DAG": "CS†"}
CONTROL_TARGET = {"CX": ("●", "⊕"), "CZ": ("●", "●"), "SWAP": ("×", "×"), "CH": ("●", "H"),
                  "CS": ("●", "S"), "CS_DAG": ("●", "S†")}


def _ops(stim):
    """(name, qubits, time step, instruction index) for each gate, as the judge schedules them,
    and the idle gaps (instruction, qubit, first step, n) of noise.schedule."""
    prog, steps, gaps = schedule(stim)
    return [(ins.name, group, step, i) for i, (ins, groups) in enumerate(zip(prog, steps))
            if ins.name != "DETECTOR" for group, step in groups], gaps


def draw(circuit, data_qubits, *, detectors=(), frame=None, postselect=(), scan=None, max_layers=None,
         noise=STANDARD):
    """Return a matplotlib Figure. `circuit` is submission text or a QuantumCircuit (then with its
    detectors, frame and postselect conditions, which the drawing checks but does not show)."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    if not isinstance(circuit, str):
        circuit = submission_text(circuit, detectors, frame, postselect)
    stim, _, errors = convert(circuit, noise=noise)
    if errors:
        raise ValueError("; ".join(errors))
    ops, gaps = _ops(stim)
    if max_layers:
        ops = [o for o in ops if o[2] < max_layers]
        gaps = [(i, q, a, min(n, max_layers - a)) for i, q, a, n in gaps if a < max_layers]
    used = sorted({q for _, qs, _, _ in ops for q in qs} | set(data_qubits))
    n_layers = max((o[2] for o in ops), default=0) + 1
    worst = max((scan.by_line.values() if scan else []), default=0) or 1

    fig, (chip, time) = plt.subplots(
        1, 2, figsize=(4 + min(0.32 * n_layers, 60), max(4.5, 0.32 * len(used) + 1)),
        gridspec_kw={"width_ratios": [4, min(0.32 * n_layers, 60)]})

    # Chip.
    load = {}
    for name, qs, _, _ in ops:
        if name in TWO:
            e = tuple(sorted(qs))
            load[e] = load.get(e, 0) + 1
    for q in range(ROWS * COLS):
        r, c = position(q)
        for nb in (q + 1 if c + 1 < COLS else None, q + COLS if r + 1 < ROWS else None):
            if nb is None:
                continue
            (r2, c2), n = position(nb), load.get((q, nb), 0)
            chip.plot([c, c2], [-r, -r2], color="#c0392b" if n else "#dddddd", lw=1 + min(n, 12) * 0.6,
                      zorder=1)
    for q in range(ROWS * COLS):
        r, c = position(q)
        if q in data_qubits:
            face, text = "#2c3e50", f"d{list(data_qubits).index(q)}"
        else:
            face, text = ("#ecf0f1" if q in used else "white"), str(q)
        chip.scatter(c, -r, s=420, color=face, edgecolor="#7f8c8d", zorder=2)
        chip.text(c, -r, text, ha="center", va="center", fontsize=7, zorder=3,
                  color="white" if q in data_qubits else "#2c3e50")
    chip.set_title("chip (d = data qubit, red = couplers used)", fontsize=9)
    chip.set_aspect("equal")
    chip.axis("off")

    # Timeline.
    row = {q: i for i, q in enumerate(used)}
    for q, i in row.items():
        time.plot([-0.5, n_layers - 0.5], [-i, -i], color="#bdc3c7", lw=0.6, zorder=0)
        time.text(-0.8, -i, f"q{q}" + (f" (d{list(data_qubits).index(q)})" if q in data_qubits else ""),
                  ha="right", va="center", fontsize=7)
    def dot(x, q, share=0.0):
        if share:
            time.scatter(x, -row[q], s=14 + 50 * share, color="#c0392b", edgecolor="black", lw=0.5, zorder=4)
        else:
            time.scatter(x, -row[q], s=6, color="#e74c3c", alpha=0.5, lw=0, zorder=4)

    seen = set()
    for name, qs, layer, line_no in ops:
        if name in TWO:
            a, b = qs
            time.plot([layer, layer], [-row[a], -row[b]], color="#2c3e50", lw=1, zorder=1)
        marks = CONTROL_TARGET[name] if name in TWO else [SHORT.get(name, name)] * len(qs)
        for q, mark in zip(qs, marks):
            time.add_patch(Rectangle((layer - 0.35, -row[q] - 0.35), 0.7, 0.7, facecolor="white",
                                     edgecolor="#2c3e50", lw=0.8, zorder=2))
            time.text(layer, -row[q], mark, ha="center", va="center", fontsize=7 if name in TWO else 6,
                      zorder=3)
        # Where the judge puts this instruction's noise: a measurement's flip before it, everything
        # else after; a first use that is not a reset is a noisy |0> just before it.
        share = scan.by_line.get(line_no, 0) / worst if scan else 0
        for q in qs:
            if noise.initial and q not in seen and name not in RESETS:
                dot(layer - 0.5, q)
            seen.add(q)
            if name in MEASURES:
                dot(layer - 0.5, q, share if noise.meas else 0)
                if noise.meas_after:
                    dot(layer + 0.5, q, share)
            else:
                dot(layer + 0.5, q, share)
    if noise.idle:   # one dot per time step a qubit waits
        for _, q, first, n in gaps:
            for x in range(first, first + n):
                dot(x, q)
    time.set_title(f"circuit as judged ({n_layers} time steps)", fontsize=9)
    legend = ("small red dot: where the judge inserts noise (after a gate or preparation, around a "
              "measurement, and each time step a qubit waits)")
    if scan:
        legend += (f"\nlarge dot: a single fault there gets past every detector and breaks the output; "
                   f"size ~ its share of c1 = {scan.c1:.3g}")
    time.text(-2.5, -len(used) - 0.2, legend, ha="left", va="top", fontsize=8)
    time.set_xlim(-2.5, n_layers)
    time.set_ylim(-len(used) - 1.2, 0.7)
    time.axis("off")
    fig.tight_layout()
    return fig


def main(argv=None):
    import argparse

    ap = argparse.ArgumentParser(description="Draw a submission on the chip and in time.")
    ap.add_argument("circuit", help="submission file (OpenQASM 2.0 with // DETECTOR lines)")
    ap.add_argument("data_qubits", help="comma-separated, code-qubit order, e.g. 1,2,0,13,6,7,12")
    ap.add_argument("--faults", action="store_true", help="run the first-order fault scan and mark the noise that breaks the output")
    ap.add_argument("--out", default="circuit.svg", help="file to write; .svg stays sharp at any zoom")
    args = ap.parse_args(argv)
    text = open(args.circuit).read()
    data = [int(q) for q in args.data_qubits.split(",")]
    scan = None
    if args.faults:
        from .faults import fault_scan
        scan = fault_scan(text, data)
        print(scan)
    draw(text, data, scan=scan).savefig(args.out)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
