# Tech Stack and Development Standards

## Architecture Overview
Nanomix is built as a high-performance hybrid application, combining the developer productivity of Python with the computational efficiency and safety of Rust.

### Core Languages
- **Python (>=3.7):** Used for high-level logic, CLI orchestration, data manipulation (Pandas/NumPy), and scientific computing (SciPy).
- **Rust (2021 Edition):** Used for performance-critical components, including the implementation of probabilistic models and heavy data processing loops.

### Bridging and Build System
- **PyO3:** Provides the bindings between Rust and Python.
- **Maturin:** The build backend for managing the compilation of Rust extensions and packaging them as Python wheels.

## Key Libraries & Tools
- **Data Science:** `pandas`, `numpy`, `scipy`.
- **Bioinformatics:** `pyranges` for genomic interval arithmetic.
- **Parallelism:** `ray` (Python) and `crossbeam` (Rust) for multi-threaded execution.
- **CLI:** `clap` (Rust) and standard Python entry points.

## Programming Paradigms

### Spec-Driven Design
This project adheres to a **Spec-Driven Design** paradigm. All significant features, architectural changes, and milestones are first documented as Markdown files in the `specs/` directory. This ensures that the codebase remains aligned with its documented intent and facilitates collaboration between technical and clinical contributors.

### Probabilistic Modeling
The core logic relies on statistical inference, specifically:
- **Maximum Likelihood Estimation (MLE)** for mixture proportion optimization.
- **Expectation-Maximization (EM)** for soft-assignment mixture models.
- **Binomial/Bernoulli Distributions** for modeling modification call success.

## Installation & Usage
### Developer Setup
```bash
pip install maturin
git clone https://github.com/simpsonlab/nanomix.git
cd nanomix
maturin develop
```
### CLI Usage
The primary entry point is the `nanomix` command, which provides subcommands for deconvolution (`deconvolute`), fragment assignment (`assign`), simulation (`simulate`), and evaluation (`evaluate`).
