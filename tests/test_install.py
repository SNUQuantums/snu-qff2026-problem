"""설치가 정상인지 확인합니다. 물리를 검증하는 것이 아니라 환경을 검증하는 테스트입니다.

    python -m pytest tests -q
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import pytest

ROOT = Path(__file__).resolve().parent.parent


def test_qiskit_version_is_the_one_we_score_with():
    """회로 작성과 라우팅에 쓰는 Qiskit 버전을 고정합니다."""
    import qiskit

    assert qiskit.__version__ == "2.5.2"


def test_xtim_is_installed():
    """xtim 은 PyPI 에 없습니다. 대회용 빌드를 problem.md 5절대로 설치하십시오."""
    xtim = pytest.importorskip(
        "xtim", reason="xtim 이 설치되지 않았습니다. problem.md 5절을 보십시오.")
    assert xtim.__version__ == "3.1.3+qff2"


def test_judge_imports_and_knows_the_chip():
    import tprep

    assert callable(tprep.check)
    assert tprep.N_QUBITS == tprep.ROWS * tprep.COLS
    assert len(tprep.couplers()) == tprep.ROWS * (tprep.COLS - 1) + tprep.COLS * (tprep.ROWS - 1)


def test_example_circuit_scores_the_documented_value():
    """problem.md 4.2절에 적힌 점수가 실제로 나오는지 확인합니다.

    순위는 유효숫자 3자리로 가릅니다. seed가 고정이라 같은 shots에서는 항상 같은 값이 나옵니다.
    """
    pytest.importorskip("xtim", reason="xtim 이 설치되지 않았습니다.")
    import tprep

    circuit = (ROOT / "examples" / "injection.qasm").read_text()
    report = tprep.check(circuit, [0, 6, 13, 2, 7, 8, 1], shots=10**6)
    assert report.passed, str(report)
    assert f"{report.score:.3g}" == "0.0252"


def test_format_gate_explains_itself():
    import tprep

    header = 'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[36];\n'
    assert any("not neighbours" in e for e in tprep.check(header + "cx q[0],q[35];", range(7)).format_errors)
    assert any("not an allowed" in e for e in tprep.check(header + "ccx q[0],q[1],q[2];", range(7)).format_errors)


def test_package_version_matches_pyproject():
    """태그를 낼 때 __version__ 과 pyproject 의 version 을 함께 올려야 합니다.

    채점 결과에 실려 나가는 judge 버전이 이 값이므로, 둘이 어긋나면 리더보드의 한 행이
    어느 코드로 채점되었는지 추적할 수 없습니다.
    """
    import tprep

    pyproject = (ROOT / "pyproject.toml").read_text()
    declared = re.search(r'^version = "([^"]+)"', pyproject, re.M).group(1)
    assert tprep.__version__ == declared
