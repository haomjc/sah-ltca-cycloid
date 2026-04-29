# SAH-LTCA-Cycloid

**Loaded Tooth Contact Analysis (LTCA) for Cycloid-Pin-Housing Reducers with Semantic Artifact Harness (SAH) Formal Verification**

This repository contains the MATLAB implementation of a comprehensive LTCA solver for cycloid speed reducers, along with formal verification artifacts (JSON Schema, Dafny contracts) that demonstrate the Semantic Artifact Harness methodology.

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

data/                       Output directory (generated at runtime)
```

## Dependencies

- **MATLAB R2020b** or later
- **Symbolic Math Toolbox** (required by `fun_output_mod.m`)
- **Optimization Toolbox** (`fsolve`, `fzero`)
- **Parallel Computing Toolbox** (`parfor`, `parpool`)

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
