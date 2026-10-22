"""First-order fault scan: which single faults get past your detectors and break the output.

Every noise channel the judge inserts (problem.md section 1.5) is replaced, one at a time, by each of
its Paulis applied with certainty, all other noise off, and run through the judge's own sampler,
final correction round included. Summing weight x acceptance x infidelity over all faults gives the
first-order coefficient c1 of the score, 1-F = c1 p + O(p^2). It is an estimate: each fault runs SHOTS
shots, and where measurement outcomes are random its acceptance and infidelity are averages over
them. c1 = 0 means no single fault broke the output in those shots; the scan names the instructions
responsible otherwise.

    from tprep import fault_scan
    scan = fault_scan(circuit, data_qubits)   # QASM text or QuantumCircuit (+ detectors=...)
    print(scan)                               # c1 and the worst instructions
"""

import os
import warnings
from dataclasses import dataclass, field
from multiprocessing import Pool

from .check import Plan, _sample, auto_frame, conditional
from .parse import count_measurements, declarations
from .noise import STANDARD, apply_noise, noise_blocks
from .qasm import convert, submission_text

SHOTS = 256  # per fault; measurement outcomes are random, and their average sets the contribution
CHANNELS = {"DEPOLARIZE1": (1, ["X", "Y", "Z"]),
            "DEPOLARIZE2": (2, [a + b for a in "IXYZ" for b in "IXYZ" if a + b != "II"]),
            "X_ERROR": (1, ["X"]), "Z_ERROR": (1, ["Z"])}


@dataclass
class FaultScan:
    c1: float
    n_faults: int
    by_instruction: list = field(default_factory=list)  # (contribution, label, worst fault)
    by_line: dict = field(default_factory=dict)          # judge's stim line -> contribution, for draw()

    def __str__(self):
        lines = [f"{self.n_faults} single faults, c1 ~ {self.c1:.4g}   (1-F ~ {self.c1:.3g} p at first order; "
                 f"estimated from {SHOTS} shots per fault)"]
        if not self.by_instruction:
            lines.append("  no single fault broke the output in these shots: fault tolerant at first order "
                         "as far as this scan can see")
        for contrib, label, worst in self.by_instruction[:15]:
            lines.append(f"  {contrib:7.4f}  {label}   e.g. {worst}")
        if len(self.by_instruction) > 15:
            lines.append(f"  ... {len(self.by_instruction) - 15} more")
        return "\n".join(lines)


def _run(args):
    body, plan, gates, w = args
    import xtim

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = _sample(xtim.Circuit(body.replace("@FAULT@", gates) + "\n" + plan.tail), SHOTS, 0, plan)
    return 0.0 if r["acceptance"] == 0 else w * r["acceptance"] * r["score"]


def fault_scan(circuit, data_qubits, *, detectors=(), frame=None, postselect=(), noise=STANDARD,
               processes=None):
    """The first-order fault scan of a submission, with its POSTSELECT conditions and logical frame
    (declared, or found as check() does). Raises ValueError if it does not pass format; run check()
    first, the scan assumes a correct circuit."""
    import xtim

    if not isinstance(circuit, str):
        circuit = submission_text(circuit, detectors, frame, postselect)
    elif frame is not None or postselect:
        raise TypeError("for submission text, write FRAME_X/FRAME_Z/POSTSELECT lines in the text")
    stim, labels, errors = convert(circuit, noise=noise)
    if errors:
        raise ValueError("; ".join(errors))
    from .exact import verify

    frame, post = declarations(stim)
    clean = xtim.Circuit(apply_noise(stim, 0.0, noise))
    frame = frame or auto_frame(clean, list(data_qubits))
    exact = verify(stim, list(data_qubits), frame=frame, postselect=post)
    if exact["errors"]:
        raise ValueError("the circuit fails the correctness check: " + exact["errors"][0])
    stim, post = conditional(stim, post, exact["conditional_detectors"])
    plan = Plan(list(data_qubits), frame, post, count_measurements(stim))

    # Noise per instruction, so each channel knows the instruction it belongs to; at p = 1 each
    # line's strength is its first-order weight in units of p. Every other channel stays in place
    # at strength 0 (see check_stim for why the lines must be present).
    blocks = noise_blocks(stim, 1.0, noise, first_order=True)[1]
    zeros = noise_blocks(stim, 0.0, noise)[1]
    jobs, owners = [], []
    for b, block in enumerate(blocks):
        for j, line in enumerate(block):
            name = line.split("(")[0]
            if name not in CHANNELS:
                continue
            arity, paulis = CHANNELS[name]
            weight = float(line.split("(", 1)[1].split(")", 1)[0])
            qs = line.split(")", 1)[1].split()
            body = "\n".join("@FAULT@" if (i, k) == (b, j) else l
                             for i, zl in enumerate(zeros) for k, l in enumerate(zl))
            for g in range(0, len(qs), arity):
                for P in paulis:
                    gates = "\n".join(f"{p}_ERROR(1) {q}" for p, q in zip(P, qs[g:g + arity]) if p != "I")
                    jobs.append((body, plan, gates, weight / len(paulis)))
                    owners.append((b, " ".join(f"{p}{q}" for p, q in zip(P, qs[g:g + arity]) if p != "I")))

    if processes == 1:
        res = list(map(_run, jobs))
    else:
        with Pool(processes or os.cpu_count()) as pool:
            res = pool.map(_run, jobs, chunksize=8)
    # The score averages over accepted runs: with POSTSELECT, divide by the noiseless acceptance.
    res = [r / exact["acceptance"] for r in res]

    per = {}
    for contrib, (b, fault) in zip(res, owners):
        if contrib > 1e-12:
            total, worst = per.get(b, (0.0, (0.0, "")))
            per[b] = (total + contrib, max(worst, (contrib, fault)))
    rows = sorted(((t, labels[b], w[1]) for b, (t, w) in per.items()), key=lambda r: -r[0])
    return FaultScan(sum(res), len(jobs), rows, {b: t for b, (t, _) in per.items()})
