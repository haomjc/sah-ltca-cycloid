# SAH-LTCA-Cycloid

**Loaded Tooth Contact Analysis (LTCA) for Cycloid-Pin-Housing Reducers with Semantic Artifact Harness (SAH) Formal Verification**

This repository contains:

1. **MATLAB solver** — a comprehensive LTCA engine for cycloid speed reducers.
2. **Formal verification artifacts** — JSON Schema, MATLAB runtime assertions, and Dafny contracts that demonstrate the Semantic Artifact Harness methodology.
3. **LaTeX manuscript source** — the complete paper source (`main.tex`, `references.bib`, `arxiv.sty`, figures, pgfplots data) designed as a **machine-readable artifact**: every claim is linked to evidence, every algorithm has a passport, every formula has dimensional bindings, and every figure has a regeneration path. The manuscript itself is the primary demonstration of the SAH framework.

## Quick Start

```matlab
% In MATLAB, navigate to src/ and run:
run_parametric_study
```

This executes 5 comparative cases (rigid / friction / stiffness / pin-displacement / full-coupled) and exports `.dat` files to `data/` for LaTeX pgfplots visualization.

## Repository Structure

```
src/
  run_parametric_study.m    Entry point: 5-case comparative study
  runTCA_analysis.m         Main engine: DCA -> stiffness -> LTCA (parfor)
  fun_force_analysis.m      Force balance solver with SAH assertion contracts
  fun_findtc.m              Initial guess for contact angle TC
  fun_output_mod.m          DCA contact geometry (requires Symbolic Toolbox)
  fun_tooth_stiffness.m     Energy method tooth body stiffness
  export_analysis_data.m    Export .dat files for LaTeX pgfplots

verification/
  sah-ltca.schema.json      JSON Schema defining SAH harness structure
  ltca_contracts.dfy        Dafny contract sketch with loop invariants

manuscript/                 LaTeX source — machine-readable paper artifact
  main.tex                  Full paper with claim-evidence maps, algorithm passports,
                            formula passports, code-artifact maps, and pgfplots figures
  references.bib            BibTeX records with dense evidence annotations
  arxiv.sty                 Document style
  graphicalAbstract.png     Graphical abstract
  data/                     pgfplots .dat files exported from MATLAB solver

data/                       Output directory (generated at runtime)
```

## Dependencies

### MATLAB (for solver)

- **MATLAB R2020b** or later
- **Symbolic Math Toolbox** (required by `fun_output_mod.m`)
- **Optimization Toolbox** (`fsolve`, `fzero`)
- **Parallel Computing Toolbox** (`parfor`, `parpool`)

### LaTeX (for manuscript)

- **TeX Live** or **MiKTeX** with `pdflatex` and `bibtex`
- Packages: `pgfplots`, `tikz`, `tcolorbox`, `natbib`, `algorithm`, `algpseudocode`, `geometry`, `fancyhdr`

## Call Chain

```
run_parametric_study.m
  +-- runTCA_analysis.m
  |     +-- fun_output_mod.m      (DCA: fsolve, symbolic)
  |     +-- fun_findtc.m          (TC initial guess: fsolve)
  |     +-- fun_tooth_stiffness.m (energy method: integral)
  |     +-- fun_force_analysis.m  (LTCA: scan + fzero + Hertz iteration)
  +-- export_analysis_data.m      (export .dat for pgfplots)
```

## Manuscript as Machine-Readable Artifact

The `manuscript/` directory contains the complete LaTeX source of the paper *"Semantic Artifact Harness for AI-Native Research"*. The manuscript is not a conventional paper—it is designed as a **typed, auditable artifact** where:

- Every **claim** is linked to a code path, figure, table, or verification artifact via a claim-evidence map.
- Every **algorithm** has a structured passport (inputs, outputs, assumptions, failure modes, convergence criteria).
- Every **formula** has a passport with variable bindings and dimensional analysis.
- Every **figure** has a regeneration path (MATLAB script → `.dat` → pgfplots → TikZ).
- Every **AI interaction** is recorded as a compressed correction trace.
- The **literature** is compressed into a dense evidence matrix mapping source → role → boundary.

This structure makes the paper inspectable by both humans and AI systems, which is the core thesis of the SAH framework. The LaTeX source can be compiled with `pdflatex` + `bibtex`; figures are generated from solver output in `data/`.

## Formal Verification

### Layer 1: MATLAB Runtime Assertions (implemented)

`fun_force_analysis.m` includes SAH-branded assertion contracts:

- **SAH_PRE** (preconditions): parameter range checks at function entry
- **SAH_POST** (postconditions): force non-negativity, backside zeroing, residual finiteness
- **SAH_CONV** (convergence warning): alerts when pin-hole iteration exhausts max_iter

### Layer 2: Dafny Contracts (sketch)

`verification/ltca_contracts.dfy` provides:

- `PinState` datatype with force, deformation, and working-side fields
- `ForceBalanceInvariant` predicate (non-negativity + backside zero)
- `ComputePinForces` method with pre/postconditions
- `ConvergenceGuarantee` lemma (axiom, with proof sketch)

### Layer 3: JSON Schema

`verification/sah-ltca.schema.json` defines:

- `ResearchCard`: hypothesis, methodology, expected outcomes
- `ClaimEvidenceMap`: typed claims with evidence bindings
- `AlgorithmPassport`: function interfaces with parameters, invariants, convergence criteria
- `FormulaPassport`: equations with variable bindings and dimensional analysis
- `VerificationReport`: pass/fail results with timestamps

## Citation

If you use this code in your research, please cite:

```bibtex
@article{sah-ltca-2026,
  title   = {Semantic Artifact Harness for AI-Native Research},
  author  = {OpenCode and Jiacheng Miao},
  journal = {aiXiv},
  year    = {2026},
  url     = {https://aixiv.science/aisc2026/}
}
```

## License

MIT License. See [LICENSE](LICENSE).
