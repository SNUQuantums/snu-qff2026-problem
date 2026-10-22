# Qiskit Fall Fest 2026 @ SNU

**Qiskit으로 원하는 quantum state를 만들어 봅시다.**
처음 배우는 분을 위한 튜토리얼 미션과, noise가 있는 환경에서 회로 성능을 겨루는 메인 해커톤을 진행합니다.

## 어디서 시작하나요?

| 프로그램 | 내용 | 자료 | 완료·제출 방법 |
|---|---|---|---|
| **튜토리얼 미션 1** | 세 qubit의 GHZ state 만들기 | [`tutorial.ipynb`](tutorial.ipynb) Part 1 | 노트북의 `TODO`를 채우고 결과 확인 |
| **튜토리얼 미션 2** | Noise가 있는 GHZ 회로 개선하기 | [`tutorial.ipynb`](tutorial.ipynb) Part 2 | 노트북의 `TODO`를 채우고 결과 비교 |
| **메인 해커톤** | 주어진 7-qubit state를 noise 아래에서 더 정확하게 준비하기 | [`problem.md`](problem.md) | 대회 사이트에 회로와 `data_qubits` 제출 |

**튜토리얼 미션 1·2는 입문 실습입니다. 메인 해커톤의 선행 자격이나 순위 점수가 아닙니다.**
두 미션을 완료하지 않아도 메인 해커톤에 바로 참여할 수 있습니다.

**AI 사용:** ChatGPT, Claude, GitHub Copilot 등 AI 도구를 자유롭게 사용해도 됩니다.
개념 학습, 코드 작성·수정, 디버깅, 회로 설계 등 모든 과정에서 사용할 수 있습니다.

## 튜토리얼 참여

시작일 Zoom 설명에서 노트북을 배포하고 함께 실습합니다. 미션 1·2는 설명 도중 차례로 안내합니다.
Zoom과 현장 설명이 끝난 뒤 레포를 공개하며, 설명에 참여하지 못한 분도 같은 노트북으로 실습할 수 있습니다.

1. 배포받은 `tutorial.ipynb`를 Google Colab에 업로드하거나 로컬 Jupyter / VS Code에서 엽니다.
2. 맨 위의 설치 셀부터 순서대로 실행합니다. 튜토리얼은 Python 3.11 이상을 기준으로 합니다.
3. 미션의 `TODO`를 채우고 셀을 실행해 회로와 출력 state, fidelity를 확인합니다.
4. 회로를 수정하고 결과를 비교해 보세요.

미션에는 IBM 계정이나 메인 채점기 설치가 필요하지 않습니다. 결과를 살펴보는 코드도 노트북 안에 있습니다.
노트북 마지막의 메인 채점기 실습은 선택 사항입니다.

## 메인 해커톤 참여

1. [`problem.md`](problem.md)에서 목표 state, 허용 연산, hardware, 채점 규칙을 확인합니다.
2. 아래 안내에 따라 레포와 채점기를 준비합니다.
3. 예시 회로를 실행하고 회로·배치를 수정해 점수를 비교합니다.
4. **OpenQASM 2.0 회로와 `data_qubits`**를 대회 사이트에 제출합니다.

| 안내 | 위치 |
|---|---|
| 최종 문제·채점 규칙 | [`problem.md`](problem.md) |
| 제출 가능한 시작용 회로 | [`examples/injection.qasm`](examples/injection.qasm) |
| 제출 사이트·순위표 | 행사 안내에 따라 공개 예정 |
| 질문 | Discord `#qff-질문` 채널 |

### 메인 채점기 설치

레포와 Release 공개 후 아래 명령을 실행합니다.

```bash
git clone https://github.com/SNUQuantums/snu-qff2026-problem.git
cd snu-qff2026-problem
pip install -r requirements.txt
```

대회용 `xtim==3.1.3+qff2`을 이 레포의 Release에서 설치합니다.
지원 환경은 **CPython 3.11–3.13**, Linux (x86_64, aarch64), macOS (Apple Silicon), Windows (64비트)입니다.
GPU는 필요하지 않습니다. 자세한 조건은 `problem.md` 5.1절에 있고, 설치가 되지 않으면 Google Colab을 사용하세요.
튜토리얼 미션 1·2에는 이 제약이 적용되지 않습니다.

### 첫 점수 확인

레포 최상위 폴더에서 실행하세요.

```python
from pathlib import Path
import tprep

circuit = Path("examples/injection.qasm").read_text(encoding="utf-8")
data_qubits = [0, 6, 13, 2, 7, 8, 1]
report = tprep.check(circuit, data_qubits)
print(report)
```

```text
p=0.001  shots=1,000,000  seed=0
  format:      PASS
  class:       PASS
  correctness: PASS
  acceptance:  PASS  (1.0000, floor 0.5)
  score 1-F = 0.0251615 ± 0.00014   (3 significant figures: 0.0252)
  2q+ gates = 28
  RESULT: PASS
```

점수는 **채점기의 최종 정정 절차를 거친 뒤의 infidelity $1-F$**이며 작을수록 좋습니다.
Gate 수와 배치를 개선하는 접근부터 시작해도 됩니다.

`check()`의 기본값은 $10^6$ shots입니다. 공식 설정은 `shots=10**8, seed=0`이며,
같은 채점기 빌드와 환경에서 재현할 수 있습니다. `±`는 Monte Carlo 추정의 표준오차입니다.
채점 코드와 자세한 규칙은 공개되어 있습니다.

### Qiskit 회로와 검출기

직접 만든 `QuantumCircuit`을 `check()`에 넘길 수도 있습니다.
`qc`, `data_qubits`, `detectors`를 준비한 뒤 다음처럼 사용합니다.
검출기가 없다면 `detectors = []`로 두세요.

```python
report = tprep.check(qc, data_qubits, detectors=detectors)
text = tprep.submission_text(qc, detectors=detectors)
```

검출기는 측정 결과로 실행을 받아들일지 정하는 조건입니다.
예를 들어 `[[qc.clbits[0]]]`은 해당 bit가 0인 실행만 받아들입니다.
정확한 문법과 acceptance 조건은 `problem.md` 1.3절과 1.6절에 있습니다.

### 회로와 오류 위치 보기

```python
scan = tprep.fault_scan(circuit, data_qubits)
print(scan)
fig = tprep.draw(circuit, data_qubits, scan=scan)
fig.savefig("my_circuit.svg")
```

`fault_scan`은 어떤 연산의 오류가 성능을 떨어뜨리는지 분석합니다.
`draw`는 qubit 배치와 회로를 그리고, noise가 들어가는 자리를 빨간 점으로, 출력을 깨는 자리를 큰 점으로
표시합니다. 자세한 설명은 `problem.md` 4.4절에 있습니다.

## 파일 구성

| 파일·폴더 | 내용 |
|---|---|
| `tutorial.ipynb` | Qiskit 입문, 튜토리얼 미션 1·2, 메인 문제 연결 실습 |
| `problem.md` | 메인 해커톤의 최종 문제와 규칙 |
| `tprep/` | 공식 채점과 로컬 검증에 사용하는 채점기 |
| `examples/injection.qasm` | 제출 가능한 시작용 회로 |
| `requirements.txt` | 메인 채점기 설치 목록 |
| `tests/` | 채점기 자체 검사 |

## 라이선스

코드(`tprep/`, `examples/`, `tests/`)는 MIT입니다. [`LICENSE`](LICENSE)를 참고하세요.
문서(`problem.md`, `README.md`)는 [`LICENSE-docs`](LICENSE-docs)를 따릅니다.
