"""tprep: the judge for the Qiskit Fall Fest 2026 @ SNU magic state preparation problem.

    from tprep import check

    report = check(circuit, data_qubits)
    print(report)

`circuit` is the submission: OpenQASM 2.0 from Qiskit with `// DETECTOR` lines (problem.md
section 1.3), or a Qiskit QuantumCircuit with `detectors=`. `data_qubits` names the seven physical
qubits holding |T_L>, in code-qubit order. The scoring server runs this same code, so
the score you get locally with `shots=10**8, seed=0` is the score you submit.
"""

from .check import Report, check, check_stim, final_round
from .faults import FaultScan, fault_scan
from .qasm import submission_text
from .grid import ACCEPTANCE_FLOOR, ALLOWED, COLS, N_CODE, N_QUBITS, ROWS, adjacent, couplers, grid_diagram
from .noise import IDLE, LEGACY, STANDARD, UNIFORM, NoiseModel, apply_noise
from .parse import check_format, count_multiqubit_gates, parse


def draw(*args, **kwargs):
    """See tprep.vis.draw. Imported on use so that `python -m tprep.vis` and matplotlib stay optional."""
    from .vis import draw as _draw

    return _draw(*args, **kwargs)


__version__ = "3.0.3"

__all__ = [
    "ACCEPTANCE_FLOOR", "ALLOWED", "COLS", "N_CODE", "N_QUBITS", "NoiseModel", "ROWS", "Report",
    "STANDARD", "IDLE", "UNIFORM", "LEGACY", "adjacent", "apply_noise", "check", "check_format", "count_multiqubit_gates",
    "couplers", "final_round", "grid_diagram", "parse", "__version__", "check_stim", "FaultScan",
    "draw", "fault_scan", "submission_text",
]
