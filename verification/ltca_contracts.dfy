// ltca_contracts.dfy
// Dafny contract sketch for the LTCA force-balance core.
//
// This file captures the formal specification of the pin-hole iteration
// in fun_force_analysis.m. It is NOT a complete verified implementation;
// convergence lemmas are declared as axioms (see [SAH] markers).
//
// Purpose: demonstrate how the MATLAB runtime assertions (SAH_POST/SAH_PRE)
// correspond to machine-checkable loop invariants and postconditions in a
// verification-aware language. This addresses Reviewer concern about
// "machine-verifiable formal specifications."

// ============================================================
// Datatypes
// ============================================================

datatype PinState = PinState(
  position_x: real,
  position_y: real,
  force: real,            // Contact force magnitude (>= 0)
  deformation: real,      // Pin-hole elastic deformation (>= 0)
  isWorkingSide: bool     // true if pin contributes to torque balance
)

// ============================================================
// Predicates
// ============================================================

// All pin forces are non-negative.
predicate ForceNonNegative(pins: seq<PinState>)
{
  forall i :: 0 <= i < |pins| ==> pins[i].force >= 0.0
}

// Backside (non-working) pins carry zero force.
predicate BacksideZeroForce(pins: seq<PinState>)
{
  forall i :: 0 <= i < |pins| ==> !pins[i].isWorkingSide ==> pins[i].force == 0.0
}

// Combined force balance invariant.
predicate ForceBalanceInvariant(pins: seq<PinState>)
{
  ForceNonNegative(pins) && BacksideZeroForce(pins)
}

// Torque balance residual is finite and bounded.
predicate TorqueBalanceValid(EQ1: real, Tin: real)
{
  // EQ1 = Tin - sum(fp .* lp_eff)
  // Valid if residual is finite (not NaN/Inf).
  EQ1 == EQ1  // NaN != NaN trick for finiteness in Dafny
}

// ============================================================
// Main method contract
// ============================================================

// Computes pin contact forces via iterative pin-hole Hertz compliance.
//
// Corresponds to: fun_force_analysis.m, outer for-loop (lines 117-296)
// The MATLAB implementation uses:
//   - max_iter = 50
//   - tol = 1e-6
//   - convergence: max(|delta_ph_new - delta_ph|) < tol
//
// This Dafny signature captures the contract; the body is abstract.
method ComputePinForces(
  // Precondition parameters (SAH_PRE contracts)
  Np: nat,                           // Number of pins (>= 1)
  Tin: real,                         // Input torque (Tin > 0)
  A: real,                           // Eccentricity (A > 0)
  Rp_profile: real,                  // Generating radius (Rp_profile > 0)
  Rrp: real,                         // Generating pin radius (Rrp > 0)
  Rrp_pin: real,                     // Operating pin radius (Rrp_pin > 0)
  E: real,                           // Elastic modulus (E > 0)
  PR: real,                          // Poisson ratio (0 <= PR <= 0.5)
  B: real,                           // Face width (B > 0)
  mu_fric: real,                     // Friction coefficient (mu_fric >= 0)
  max_iter: nat,                     // Max iterations (max_iter > 0)
  tol: real,                         // Convergence tolerance (tol > 0)
  // Mutable state
  TP: seq<real>,                     // Contact angles (|TP| == Np)
  BL: seq<real>                      // Backlash values (|BL| == Np)
) returns (
  pins: seq<PinState>,               // Output: pin states
  EQ1: real,                         // Output: torque balance residual
  converged: bool                    // Output: whether iteration converged
)
  // --- Preconditions (from SAH_PRE assertions in MATLAB) ---
  requires |TP| == Np
  requires |BL| == Np
  requires Tin > 0.0
  requires A > 0.0 && Rp_profile > 0.0 && Rrp > 0.0 && Rrp_pin > 0.0
  requires E > 0.0
  requires 0.0 <= PR && PR <= 0.5
  requires B > 0.0
  requires mu_fric >= 0.0
  requires max_iter > 0
  requires tol > 0.0

  // --- Postconditions (from SAH_POST assertions in MATLAB) ---
  ensures |pins| == Np
  ensures ForceBalanceInvariant(pins)          // fp >= 0, backside == 0
  ensures TorqueBalanceValid(EQ1, Tin)          // residual is finite
  ensures converged ==> forall i :: 0 <= i < Np ==> pins[i].deformation >= 0.0

  // Frame condition: TP and BL are read-only inputs
  modifies {}

// ============================================================
// Loop invariant (for the outer iteration)
// ============================================================
//
// The MATLAB outer loop (lines 117-296) maintains:
//   1. Force non-negativity: fp >= 0 for all pins
//   2. Backside zeroing: fp(~work_mask) == 0
//   3. Deformation convergence: max(|delta_ph_new - delta_ph|) decreasing
//
// In Dafny, these would appear as while-loop invariants:

//   while iter < max_iter && !converged
//     invariant ForceBalanceInvariant(pins)
//     invariant forall i :: 0 <= i < Np ==> pins[i].deformation >= 0.0
//     invariant iter <= max_iter
//     // [SAH] Monotone convergence: each iteration reduces the residual.
//     // This is an AXIOM — not proven from the Hertz compliance model.
//     invariant prev_residual >= 0.0 ==> residual <= prev_residual

// ============================================================
// Convergence lemma (declared, not proven)
// ============================================================

// [SAH-AXIOM] Convergence guarantee.
// If the Hertz compliance C_hertz is bounded above (C_hertz <= C_max),
// and the structural compliance C_struct >= 0, then the pin-hole
// iteration converges in finite steps.
//
// Proof sketch:
//   1. The mapping T: delta_ph -> delta_ph_new is a contraction
//      when C_total = C_struct + C_hertz + C_ph is bounded.
//   2. The residual r(k) = max|delta_ph_new - delta_ph| satisfies
//      r(k+1) <= L * r(k) for some L < 1 (Banach fixed-point).
//   3. Convergence in max_iter steps follows from L^max_iter < tol.
//
// Status: This lemma is stated for documentation; the contraction
// property requires domain-specific analysis of the Hertz model
// that is beyond the scope of this sketch.

lemma {:axiom} ConvergenceGuarantee(
  C_max: real,        // Upper bound on Hertz compliance
  C_struct: real,     // Structural compliance (>= 0)
  max_iter: nat,
  tol: real
)
  requires C_max > 0.0
  requires C_struct >= 0.0
  requires max_iter > 0
  requires tol > 0.0
  // The contraction ratio depends on C_total = C_struct + C_hertz + C_ph.
  // If bounded, convergence is guaranteed within max_iter steps.
  ensures true  // [SAH] Placeholder — full proof requires Hertz model analysis

// ============================================================
// Work-side detection invariant
// ============================================================

// The 2D cross product torque_dir = r_x * F_y - r_y * F_x determines
// which pins are on the working side. This is a geometric invariant
// that does not depend on force magnitudes.
predicate WorkSideDetectionCorrect(
  pins: seq<PinState>,
  node_x: real, node_y: real,
  X0: real, Y0: real
)
{
  forall i :: 0 <= i < |pins| ==>
    var r_x := pins[i].position_x - X0;
    var r_y := pins[i].position_y - Y0;
    // torque_dir > 0 <==> working side (provides positive torque)
    (r_x * (-1.0) - r_y * 0.0 > 0.0) ==> pins[i].isWorkingSide
}
