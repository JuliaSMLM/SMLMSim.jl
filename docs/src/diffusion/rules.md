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

**Range.** Coordinates are micrometres and times seconds. Every number of a `DiffusionSMLMConfig`
but `r_react` (mixture fractions included), every velocity component (the config's, `track_velocity`
or a saved one), every per-molecule D (`track_D` or a saved one) and `simulate`'s `γ` is 0 or has a
magnitude from `1e-50` to `1e50`, as given, before any conversion (the deprecated `photons` keyword is
not checked, as in 0.7, until its removal in 0.8.0); `k_off` may also be
infinite. A run has at most `2^53` steps and `2^53` molecules, and the default camera at most `2^53` pixels
per side (counts compared exactly). Construction and `simulate` check all of this, the config
again on its fields as they are when `simulate` runs, and throw `ArgumentError` otherwise. In that
range the products, squares and step scales the simulation computes, such as `sqrt(2·D·dt)`, `v·dt`,
`γ·dt` and the `:min` pair coefficient, are normal Float64 numbers.

`r_react` may be any value `>= 0`, infinite included. Binding compares it with a distance that is
scaled against underflow; a dissociation spacing it makes too large to place (infinite included)
fails the fit test, and the partners stay where they are (Dissociation, below).

Positions given in `starting_conditions` must be finite (`ArgumentError` otherwise), and a coordinate in
`[0, box_size]` is used as given. Under periodic boundaries a position outside the box names its image in the
box: each coordinate outside is wrapped at entry, as the first step would wrap it, so no record lies outside the
box. Under reflecting boundaries a position beyond a wall throws `ArgumentError`, except a coordinate up to the box
written in its own type (`Float32(box_size)` can round above `box_size`), which is clamped inside. Coordinates, photons and
timestamps are stored in each emitter's own type: with Float32 starting conditions (all or some),
`simulate` throws `ArgumentError` when the box, the run's duration or the photons per record a `γ`
restamp gives does not fit Float32 (about `3e38`); photons kept from the source, or set by the
deprecated `photons` keyword, are not checked, as in 0.7. A diffusion step beyond that range overflows
the coordinate. Molecules of different types may pair: a pair forms, moves and dissociates with each partner
placed in Float64 and stored in its own type, and partners split closer than `r_react` are spaced with the
margin of the coarser type. A coordinate that nothing moves keeps its stored value: a diffusion step with
D = 0 and the boundary leave a coordinate inside the box as it is, and an anchor keeps its position, so an
immobile BigFloat molecule keeps every bit, even beyond the current `setprecision`. A step that moves a BigFloat
coordinate computes it at the current precision.

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
- When the pair cannot fit (possible only when `box_size < d_dimer`), each partner is reflected on its own, as a
  monomer is (Directed motion, A wall).

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
pair moves with `diff_dimer` and rotates with `diff_dimer_rot`, even when a member is immobile (monomer D = 0 and no velocity).
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

## Directed motion

A molecule can carry a constant velocity `v` in μm/s. The config's `monomer_velocity` gives one velocity to every
monomer, or one per `monomer_mobility` population, drawn with the population's D. `simulate`'s `track_velocity` gives
one to single molecules.
- **A monomer step** is the diffusion step (Gaussian, variance `2·D·dt` per axis, the same draws as without a
  velocity), then the displacement `v·dt`, then the boundary. A zero velocity takes exactly the step without one, so
  the output and the random stream are unchanged.
- **A wall.** In a reflecting box the position is reflected off the wall it crossed, as for any monomer step, and the
  velocity is not reversed, so a molecule driven into a wall stays within one step's displacement of it: it piles up
  at the wall, as a constant force against a wall that passes nothing does. Near a wall the drift is resolved to about
  `sqrt(2·D·dt)`, so the layer at the wall is resolved better at smaller `dt`. A step longer than the box folds back as
  often as it crosses a wall (before 0.7.4 such a step stopped at the wall, with or without a velocity). Under
  periodic boundaries it wraps. The drift and the boundary are applied in Float64 to the stored position after the
  diffusion step, which is converted once more to its coordinate type (so a Float32 molecule with D > 0 and a
  velocity is rounded twice, by at most one unit in the last place), and a drift within the model's range (Range,
  above) leaves the molecule finite and inside the box, with Float32 coordinates too.
- **Binding.** The pair search runs at the start of each step, on the positions the previous step left, as without a
  velocity: two monomers closer than `r_react` there form a pair in that step instead of moving. So molecules that
  come within `r_react` at time `t` are still monomers at `t` and a pair from `t + dt`; a record shows the pair only
  if one is taken while it lasts (a pair that forms and dissociates between records never shows). Two molecules whose
  relative step `|v1 - v2|·dt` exceeds `2·r_react` can pass through each other without binding.
  `r_react = 0` turns binding off: no pair forms, so molecules pass through each other even at one point, while
  pairs given in `starting_conditions` still dissociate at `k_off`.
- **Float32 coordinates** move in whole units of their precision: a drift step `v·dt` below the coordinate's
  resolution is rounded, and a molecule with D 0 and a very small velocity can move slower than `v`, or not at all.
- **A bound pair** drifts at the pair velocity, from overdamped force balance (a partner without a velocity counts as
  zero). Under `:min`, where the pair's diffusion follows the slower partner, the partners' velocities are weighted
  by their drag `1/D`: `v = (v1·D2 + v2·D1)/(D1 + D2)`. So a partner with D 0 carries the pair at its own velocity,
  a still partner holds it (the anchor), and a partner of very small D nearly does; with both at D 0, the plain
  mean. Under `:fixed`, whose pair diffusion does not depend on the partners, the plain mean (equal drag). The
  velocities move only the pair's center and apply no torque: partners moving in opposite directions across the bond
  neither move nor turn the pair. The drift is added to the pair's center after its diffusion and rotation draws. On
  the formation step a mobile pair with a velocity is placed about its partners' midpoint moved by `v·dt` (no
  diffusion on that step), so re-binding does not cost the pair its drift. An anchored pair does not move.
- **Immobile** means monomer D 0 and no velocity, everywhere these rules or `pair_mobility` say immobile: the `:min`
  anchor at formation, while bound and at dissociation, the `:fixed` warning, and `frame_dimer_truth`'s `mixed`. A
  molecule with D 0 and a velocity is mobile. Under `:min` its pair has D 0 (the pair rule is by D) and moves at its
  velocity.
- **Single molecules.** `track_velocity` and `track_D` replace the velocity and D of the molecules they name.
  Populations are still drawn for every molecule, so the other molecules' draws are unchanged.
- **Saved.** When the run has any velocity (the config's, a keyword's or a continued run's),
  `metadata["monomer_velocity"]` (`track_id => velocity`) holds the velocity every molecule moved with, zero included,
  and `metadata["monomer_velocity_config"]` a copy of the config's `monomer_velocity`. A run without any velocity has
  neither key, and its molecules moved with none.

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

**Velocity.** A velocity is a property of the molecule, as a saved D is. A resumed track keeps its saved velocity, with
one warning when this config's `monomer_velocity` differs from the source run's saved copy
(`metadata["monomer_velocity_config"]`), as for D. A track without a saved velocity (its run had none) takes this
config's: its population's entry, by the class it drew or kept with its saved D (none when that class is not one of
this config's populations), or the one velocity of an empty mixture, else none. A kept D keeps its saved class, or
takes the source run's population with that D; a D that is no population's (one set by `track_D`, or a source without
a mixture) gives no class, so no population's velocity. `track_velocity` and `track_D`
replace any of these for the molecules they name.

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
`monomer_D`, `monomer_class` or `monomer_velocity`, and warns once, so every molecule keeps its photons per record and
takes the current D and velocity settings. The tie is only a backstop: in a Float32 SMLD from one run, timestamps `k·dt` can collide after about 8.4e6
steps. Such a run then shows a false tie, which only drops its γ rate, saved D and saved velocity, with the warning.

### Definitions

- **Rate source:** `metadata["rate_source"]` is always one of the Strings `"γ"`, `"photons"`, `"default"`, set by
  `simulate` from its own keywords. A continuation that does not restamp stores `"default"` for a default source and
  `"photons"` otherwise; any other or missing saved value, including a non-String, reads as `"photons"`.
- **Unknown provenance** is decided once, in `extract_end_state`, which every SMLD `starting_conditions` passes
  through. Its output then lacks the five keys above, so a second extraction and `simulate` take it at face value. Its
  output's metadata holds copies (including `simulation_parameters`), so editing the extract does not change the
  source.
- **Saved mixture:** `create_smld` stores `metadata["monomer_mobility"]`, a copy of the config's mixture at run time,
  so an in-place edit of the mutable config after the run cannot hide a change.
- **Warnings:** each names what was kept or dropped and why, and is shown at most once per Julia session
  (`maxlog=1`).
