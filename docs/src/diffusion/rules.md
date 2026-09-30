# Placement and Continuation Rules

These rules state what the diffusion simulation guarantees when a pair forms, moves while bound and dissociates,
and when a run continues from an earlier one. The code cites them by section, and the tests assert each clause.

## Pair placement

### Distance

Formation and dissociation measure one distance between two partners, from their stored coordinates in Float64:
- per axis, the difference of the two coordinates, rounded once. Under periodic boundaries it is the minimum image of
  that difference (less `box_size·round(d/box_size)`), still rounded once: the subtraction's rounding error is carried
  into the result;
- then the Euclidean norm of those differences, scaled so that no normal difference squares to zero.

It is the same in either partner order. A pair forms when this distance is below `r_react`, so two monomers across a
periodic edge form a pair as anywhere else (0.7.1 measured the plain difference, and such a pair never formed).
Dissociation measures the same distance, so a pair it leaves in place is not within formation distance, unless its
placement cannot fit (Dissociation, below).

**Orientation.** Formation, a bound step and a split each orient the pair along this
displacement, in Float64 (the per-axis differences above, never a position rounded to the
coordinate type), and each leaves both partners finite and inside `[0, box_size]` in the
coordinate type.
- **Formation** places the partners `d_dimer` apart by minimum image, to within rounding,
  unless the placement cannot fit (A bound pair and the box, below): an anchored pair when
  some axis fits neither way, a mobile pair in a reflecting box narrower than `d_dimer`.
- **A bound step** keeps the bond `d_dimer` by minimum image, to within rounding, turned by
  the rotation draws from the displacement's direction; a mobile pair in a reflecting box
  narrower than `d_dimer` has each end reflected on its own.
- **A split** of partners closer than `r_react` leaves them about `r_react` apart (the
  separation `s` of Dissociation, below), not `d_dimer`, unless the placement cannot fit.

The direction is exact to within rounding whenever the displacement's components are normal
numbers of the coordinate type; for a subnormal displacement any finite unit direction is
acceptable. A zero displacement (partners at exactly the same position) is taken as +x, with
no random draw of its own.

**Range.** Coordinates are micrometres. Displacements, positions and box sizes beyond about
`1e150` in magnitude are outside the model (their squares overflow Float64), and these rules
promise nothing there.

### A bound pair and the box

**An anchored pair** (`pair_mobility = :min`, one partner immobile; if both are immobile, the lower `track_id` is the
anchor):
- The anchor never moves.
- The other partner is placed inside the box in its coordinate type, at the bond length `d_dimer` from the anchor (the
  minimum-image distance under periodic boundaries), to within rounding.
- It lies along the anchor-to-partner axis, with every axis that would leave a reflecting box mirrored across the
  anchor.
- When some axis fits neither way (possible only when `box_size < 2·d_dimer`), the partner keeps its current position
  and the pair forms at its current separation.

**A mobile pair in a reflecting box** is a rigid body.
- At formation, when an end would leave the box, both partners move to `d_dimer/2` either side of their midpoint,
  along their axis, with the midpoint shifted inward just enough to fit.
- While bound, the pair's center reflects off the walls, moved in by each end's half-extent, keeping the pair's
  orientation and bond length. The center folds into its interval as many times as it crosses a wall (a triangle
  wave).
- So both partners stay inside the box, `d_dimer` apart to within rounding.
- When the pair cannot fit (possible only when `box_size < d_dimer`), each partner is reflected on its own, as in
  0.7.1.

**Every fit** is decided against the physical box `[0, box_size]`. Formation and a bound
step compute positions in Float64 from the stored coordinates and the partner's Float64
minimum image, apply the boundary there, and convert each placed coordinate to the
coordinate type once, clamped into the box; a split places its partners the same way
(Dissociation, below). A bound step carries its orientation through the boundary: the center
folds or wraps, and the direction is never rebuilt from rounded ends.

**Under periodic boundaries**, a bound pair's center and axis are taken from its partner's
minimum image. The moved center is wrapped into the box and each end is then wrapped, so a
pair straddling the boundary moves by one step. A pair forming under periodic boundaries is
placed from its partner's minimum image too (the anchor-to-partner axis, or the midpoint and axis), and each partner
is wrapped into `[0, box_size]` in its coordinate type on the formation step. An anchored pair already is; a mobile
pair keeps its 0.7.1 placement (`d_dimer/2` either side of its midpoint, the minimum image's) up to that wrap.

**Every boundary step**, a monomer's included, ends inside the box in the coordinate type.

**Under the default `:fixed`**, forming a pair places both partners `d_dimer` apart about their midpoint, and a bound
pair moves with `diff_dimer` and rotates with `diff_dimer_rot`, even when a member is immobile (monomer D = 0).
`simulate` warns once per run when any step moves an immobile member (at formation, while bound or at dissociation,
recorded by the camera or not), pointing to `:min`, which keeps such a pair in place.

### Dissociation

**Which pairs move.** The partners' distance is the one above, with the minimum image under periodic boundaries.
Partners at least `r_react` apart by that distance keep their positions. Partners closer than `r_react` end inside the
box in their coordinate type and at least `r_react` apart by that distance, unless the construction cannot fit
(below); the same minimum image gives the pair's midpoint and axis.

So a pair keeps its positions whenever it is bound at `d_dimer` and `d_dimer` exceeds `r_react` by more than rounding,
which holds in every box of at least `2·d_dimer`. In a smaller box a pair can be closer than `d_dimer` (a reflecting
anchored formation that kept its current separation, or a periodic bond whose minimum image is shorter), and such a
pair closer than `r_react` moves at the split.

**Where they go.** The separation `s` is `r_react` plus a margin that survives rounding to the coordinate type.
- Under `pair_mobility = :min`, an immobile partner never moves, and the other is placed `s` from it by the anchored
  construction above. If both are immobile, the lower `track_id` stays and the higher one moves.
- Otherwise, including every pair under the default `:fixed`, both move to `s/2` either side of their midpoint along
  their axis, with the midpoint shifted inward just enough to fit a reflecting box, or each end wrapped under
  periodic boundaries.

**When it cannot fit**, both partners stay where they are (the 0.7.1 behaviour), and `simulate` warns once per run.
The construction cannot fit when its fit test fails:
- anchored: some axis fits neither way, possible only when `box_size < 2s`;
- about the midpoint in a reflecting box: some `s·|u_k| > box_size`;
- periodic: some `s·|u_k| > box_size/2`.

**Under `:fixed`**, a split that moves an immobile member counts toward the one warning per run above.

### Definitions

- **Inside the box in the coordinate type T:** every coordinate `c` satisfies `0 <= c <= box_size`, compared exactly
  (T promoted to Float64). A computed Float64 coordinate is converted to T and clamped into `[0, t]`, where `t` is the
  largest value of T not above `box_size`; that moves it by at most one unit in the last place of T.
- **Fits** (the only fallback trigger):
  - reflecting: each axis k has `a_k + d·u_k` or `a_k − d·u_k` in `[0, box_size]` (anchored), or `d·|u_k| <= box_size`
    (midpoint and rigid-pair forms);
  - periodic: `d·|u_k| <= box_size/2` on every axis, so the minimum image of the placed offset is the offset itself.
- **Rounding tolerance** of "at the bond length": `|dist − d_dimer| <= 4·√N·eps(T)·max(1, box_size)`. The margin in
  `s` is `8·eps(T)·max(1, box_size + r_react)`, and at least `1e-9·r_react`. It is larger than the rounding of both
  partners, so the separation after conversion is still at least `r_react`.

## Continuation

These rules apply when `starting_conditions` is an SMLD or the output of `extract_end_state`.

**Brightness.** Without an explicit `γ`, a continuation restamps every resumed molecule to `γ·dt` at the current `dt`
only when both hold:
- the source's metadata says its rate was set with γ (`rate_source == "γ"`, with its `γ` and `dt`);
- every resumed molecule carries `γ·dt_saved` photons, to a relative 1e-6.

In every other case (a default or `photons` source, a 0.7.1 SMLD, a Vector, or a γ source with any molecule off that
rate because it was edited or lacks a saved `dt`), each molecule keeps its photons per record, as in 0.7.1. There is
one warning when the source claimed a γ rate that some molecule does not carry. An explicit `γ` restamps every
molecule to `γ·dt`.

**D.** A resumed track keeps a saved D when the source saved one for it, including in a filtered, time-cut or
photon-edited subset of the run. Every other track takes the current setting at run time:
- a fresh draw from a non-empty `monomer_mobility`, saved as for a new molecule; or
- `diff_monomer`, never saved (and absent from `metadata["monomer_D"]`), so a later change applies.

A kept saved D warns once when this config's `monomer_mobility`, weights included, differs from the source run's saved
copy.

**dt.** The source's saved `dt` is read only by the brightness check (`photons ≈ γ·dt_saved`, rtol 1e-6). A restamp
uses the current `dt`, and a source that keeps its photons per record ignores both.

**Scope.** Continuation assumes the SMLD comes from one simulation run, or a filtered subset of one; continuing a
concatenation of different runs is unsupported. The γ check is a consistency check, not proof of provenance: it cannot
detect a molecule from another run that happens to carry `γ·dt_saved`.

**Provenance.** An SMLD has unknown provenance when either holds:
- its metadata carries `concatenated_from` or `merged_from` (written by SMLMData's `cat_smld` and `merge_smld`);
- its last frame holds two records of one track at the same timestamp (a backstop for a hand-built concatenation,
  which no single run produces).

Then `extract_end_state` resumes the latest record per track (the first of tied ones), carries no `γ`, `rate_source`,
`monomer_D` or `monomer_class`, and warns once, so every molecule keeps its photons per record and takes the current D
setting. The tie is only a backstop: in a Float32 SMLD from one run, timestamps `k·dt` can collide after about 8.4e6
steps. Such a run then shows a false tie, which only drops its γ rate and saved D, with the warning.

### Definitions

- **Rate source:** `metadata["rate_source"]` is always one of the Strings `"γ"`, `"photons"`, `"default"`, set by
  `simulate` from its own keywords. A continuation that does not restamp stores `"default"` for a default source and
  `"photons"` otherwise; any other or missing saved value, including a non-String, reads as `"photons"`.
- **Unknown provenance** is decided once, in `extract_end_state`, which every SMLD `starting_conditions` passes
  through. Its output then lacks the four keys above, so a second extraction and `simulate` take it at face value. Its
  output's metadata holds copies (including `simulation_parameters`), so editing the extract does not change the
  source.
- **Saved mixture:** `create_smld` stores `metadata["monomer_mobility"]`, a copy of the config's mixture at run time,
  so an in-place edit of the mutable config after the run cannot hide a change.
- **Warnings:** each names what was kept or dropped and why, and is shown at most once per Julia session
  (`maxlog=1`).
