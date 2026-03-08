# Regular Polytopes in Minkowski 2+1 Spacetime

**Date:** 2026-03-08  
**Repository:** hyperbolic_tiling  
**Related discussion:** [StackExchange: "Regular polytopes in Minkowski spacetime"](https://math.stackexchange.com/questions/1496720/regular-polytopes-in-minkowski-spacetime)

---

## 1. Minkowski ℝ²·¹ as the Primary Space

Working directly in flat Minkowski space ℝ²·¹ (signature +,+,−), there
are three types of "spheres" (loci at constant Minkowski distance from
the origin):

| Locus | Equation | Topology | Induced signature | Name |
|-------|----------|----------|-------------------|------|
| Two-sheeted hyperboloid | t²−x²−y² = +1 | Two copies of ℝ² | (+,+) Riemannian | H² (hyperbolic plane) |
| Light cone | t²−x²−y² = 0 | Cone | Degenerate | Ideal boundary |
| One-sheeted hyperboloid | x²+y²−t² = +1 | ℝ × S¹ | (+,−) Lorentzian | dS₂ (de Sitter 1+1) |

The stabilizer of a point on each locus determines how many edges can
meet symmetrically at a vertex:

| Vertex causal type | Stabilizer | Compact? | Symmetric edges | Vertex figure |
|-------------------|------------|----------|-----------------|---------------|
| Timelike (on H²)  | SO(2) (rotation) | Yes | Finite      | Polygon {q}   |
| Lightlike (on cone)| ℝ (parabolic)   | No  | Infinite    | Apeirogon {∞} |
| Spacelike (on dS₂) | SO(1,1) (boost) | No  | Infinite    | Apeirogon {∞} |

**Key insight:** An infinite vertex figure at a spacelike vertex is not
pathological — it is forced by the non-compact boost stabilizer, just as
{3,∞} has infinitely many triangles at each ideal vertex because the
parabolic stabilizer is non-compact.

---

## 2. Mirrors and Reflections: The 1+1 Picture

The simplest illustration of Minkowski mirror geometry is in ℝ¹·¹
(one space + one time dimension). A "mirror" is a line through the origin.

**Reflection through a spacelike line** (e.g., the x-axis, normal = (1,0,0)):
- Preserves each branch of t²−x² = 1 separately
- Vertices stay on one sheet of the two-sheeted hyperboloid

**Reflection through a timelike line** (e.g., the t-axis, normal = (0,0,1)):
- Swaps the two branches: (t,x) → (−t,x)
- Vertices appear on both sheets

```
        t
        |      /
        |    /          Light cone: t = ±x
        |  /
   -----+-----→ x       Spacelike mirror ↔ preserves sheets
        |  \             Timelike mirror  ↔ swaps sheets
        |    \
        |      \
```

The **orthogonal complement** (in Minkowski metric) of a spacelike subspace
is timelike, and vice versa. This is the mechanism behind the StackExchange
2-to-1 correspondence.

---

## 3. The StackExchange Result (Corrected)

The SE answer establishes that for every regular hyperbolic tiling {p,q},
there are **two** Minkowski polytopes, differing in the choice of the face
mirror R₀:

| R₀ choice    | Mirror normal type | Effect | Vertices on |
|--------------|--------------------|--------|-------------|
| Direct       | Spacelike normal   | Preserves each sheet | One sheet of two-sheeted hyperboloid |
| Complement   | Timelike normal    | Swaps sheets | Both sheets of two-sheeted hyperboloid |

**Important:** For compact tilings like {7,3}, the R₀-complement does
**NOT** produce vertices on the one-sheeted hyperboloid. All Coxeter
mirror normals for {7,3} are spacelike, and the vertex is at a timelike
(ordinary) point. The complement merely swaps sheets: (1,0,0) → (−1,0,0).
Vertices remain on t²−x²−y² = +1.

Computation confirms (see `great_heptagons_mirrors.wls`):
```
m1 = (0, 0, 1),      HNorm² = −1  (spacelike normal)
m2 = (0, −sin π/7, cos π/7), HNorm² = −1  (spacelike normal)
m3 = (derived),       HNorm² = −1.31 (spacelike normal)
Vertex p2 = (1, 0, 0), HNorm² = +1  (timelike, on two-sheeted hyperboloid)
```

The SE answer mentions the one-sheeted hyperboloid as a theoretical vertex
locus for ultraideal points, but provides **no concrete example** — all its
GIF visualizations ({∞,3}, {5,4}, {3,∞}) are on the two-sheeted hyperboloid
or light cone.

---

## 4. Stellation: A New Path to the One-Sheeted Hyperboloid

### 4.1 The stellation depth gradation

For a regular p-gon with p > 5, stellation at depth k replaces each face
{p} with a star polygon {p/k}. The vertices move outward from the face
center by the ratio tan(kπ/p)/tan(π/p):

**The {7,3} family:**

| k | Face | Ratio | Vertex locus | Tiling |
|---|------|-------|--------------|--------|
| 1 | {7}  | 1.00  | Two-sheeted hyperboloid (timelike) | {7,3} |
| 2 | {7/2}| 2.25  | Two-sheeted hyperboloid (timelike) | {7/2,7} |
| 3 | {7/3}| 5.05  | **One-sheeted hyperboloid (spacelike)** | {7/3, ∞} |

As k increases, vertices cross the light cone and emerge on the one-sheeted
hyperboloid. This is a **constructive geometric mechanism** for producing
vertices on dS₂ — fundamentally different from the SE's R₀-complement swap
(which only shuffles sheets of the two-sheeted hyperboloid).

### 4.2 The mirror signature change

When reconstructing the level-3 stellated tiling from a seed face, one of
the Coxeter mirrors changes causal type:

- **Unstellated {7,3}:** All three mirrors have spacelike normals. Spacelike
  mirrors alone reconstruct only half the faces.
- **Stellated {7/3}:** The third mirror acquires a timelike normal. Adding
  this timelike mirror recovers the full tiling.

This is discovered computationally in `great_heptagons_mirrors.wls`. The
timelike mirror is essential — without it, only half the structure is
generated.

### 4.3 The convex hull: "great heptagons"

Taking the convex hull of each {7/3} face yields a larger {7} (a "great
heptagon"), computed in `stellate.wls`:

```mathematica
convexHullOrder = Table[Mod[-6/stellationLevel * k, 7] + 1, {k, 0, 6}];
```

These great heptagons are ~5× larger than the original faces. They are
**entirely ultraideal** — all vertices lie on the one-sheeted hyperboloid.
The tiling by great heptagons corresponds to the known {7,7/2} (great
heptagonal tiling) in the Polytope Wiki, with central density 3 and
imaginary circumradius ≈ 0.87306i.

---

## 5. Properties of the {7/3} Tiling as a Minkowski Object

### 5.1 Infinite vertex figure

The vertex figure grows without bound as more of the tiling is included:

| Faces included | Vertex figure segments |
|---------------|----------------------|
| ~10           | 5                    |
| ~30           | 9                    |
| ~170          | 13                   |
| ~340          | 17                   |

This is consistent with the non-compact stabilizer SO(1,1) at spacelike
vertices (see Section 1).

### 5.2 Edges cross the light cone

Each edge connects two spacelike vertices (on the one-sheeted hyperboloid)
via a straight line through the interior:

```
vertex (spacelike, outside light cone)
  │
  ├── light cone crossing
  │
  ├── middle segment (timelike, inside light cone)
  │
  ├── light cone crossing
  │
vertex (spacelike, outside light cone)
```

The edges are neither "in H²" nor "in dS₂" — they are geodesics of flat
Minkowski space. From the hyperbolic viewpoint, edges "cross the ideal
boundary." From the Minkowski viewpoint, they are ordinary straight lines.

This mixed causal profile is why the {7/3} tiling is best understood as
a **Minkowski object**, not as a hyperbolic tiling with ultraideal
decorations.

### 5.3 The Schläfli symbol and terminology

| Tiling | Face | Vertex figure | Coxeter group |
|--------|------|---------------|---------------|
| {7,3}  | {7}  | {3}           | [7,3]         |
| {7/2,7}| {7/2}| {7}           | [7,3]         |
| {7,7/2}| {7}  | {7/2}         | [7,3]         |
| {7/3, ∞} | {7/3} | apeirogon  | [7,3] (extended) |
| {7, ∞} (convex hull) | {7} | apeirogon | [7,3] (extended) |

In McMullen-Schulte terminology, these are **apeirohedra** (infinite
polyhedra). Every hyperbolic tiling is already an apeirohedron; {7/3, ∞}
is additionally non-locally-finite. But from the Minkowski-native viewpoint,
"infinite vertex figure" is just the geometry of non-compact stabilizers —
there is nothing pathological about it.

---

## 6. Complete Landscape of Regular Structures in ℝ²·¹

The Polytope Wiki states: **only two infinite families of non-dense regular
star tilings** exist in the hyperbolic plane: {p/2, p} and {p, p/2} for
odd p > 5. All deeper stellations produce **dense** arrangements (see
[Dense polytope](https://polytope.miraheze.org/wiki/Dense_polytope),
[Regular polyhedron](https://polytope.miraheze.org/wiki/Regular_polyhedron)).

### Category A: Ordinary hyperbolic tilings (well-understood)

Every {p,q} with (p−2)(q−2) > 4. Each gives 2 Minkowski tilings via the
SE correspondence (R₀ direct and complement). Finite faces, finite vertex
figures. Infinitely many families.

### Category B: Non-dense star tilings (well-understood)

{p/2, p} and {p, p/2} for odd p ≥ 7. Central density 3. Each also gives
2 Minkowski tilings. Finite faces, finite vertex figures, self-intersecting
but discrete.

| Tiling | Face | Vertex fig. | Coxeter group |
|--------|------|-------------|---------------|
| {7/2, 7} | {7/2} | {7} | [7,3] |
| {7, 7/2} | {7}   | {7/2} | [7,3] |
| {9/2, 9} | {9/2} | {9} | [9,3] |
| {11/2, 11} | {11/2} | {11} | [11,3] |
| ... | | | |

For even p, these degenerate into compounds of three simpler tilings.

### Category C: Deep stellations — this repository (new)

Stellations {p/k} with k ≥ 3 push vertices onto the one-sheeted hyperboloid
(dS₂). The Polytope Wiki classifies these as "dense" from the hyperbolic
viewpoint. From the Minkowski viewpoint, they are discrete tilings with
infinite vertex figures forced by non-compact stabilizers.

This repository explores {7/3} as the first example. Other candidates:
- {9/4, ?}: stellation of {9,3} at depth 4
- {11/3, ?}, {11/4, ?}, {11/5, ?}: deeper stellations of {11,3}
- In general: {p/k, ∞} for p ≥ 7, gcd(p,k) = 1, k ≥ 3

### Category D: Timelike faces (open)

Faces with Lorentzian signature (+,−). Vertices on hyperbolas, giving
infinitely many vertices per face. No examples are known. Whether
flag-transitive tilings with timelike faces exist is an open question.

### Summary

| Category | Vertex figure | Vertex locus | Status |
|----------|---------------|--------------|--------|
| A. {p,q} | Finite | Two-sheeted hyperboloid | Well-understood |
| B. {p/2,p}, {p,p/2} | Finite | Two-sheeted hyperboloid | Well-understood |
| C. {p/k, ∞}, k≥3 | Infinite | One-sheeted hyperboloid | **This repo** |
| D. Timelike faces | Infinite | ??? | Open |

---

## 7. Convex Counterparts in Minkowski Space

### 7.1 The Euclidean pattern

In Euclidean/hyperbolic geometry, every star tiling has an "army" — the
convex tiling sharing the same vertices:

| Star tiling | Army (convex hull) |
|-------------|--------------------|
| {5/2, 5}    | {3, 5} (icosahedron) |
| {7/2, 7}    | {3, 7} (order-7 triangular tiling) |
| {7, 7/2}    | {3, 7} |

### 7.2 The Minkowski problem

In Minkowski space, "convex hull" is not well-defined — the light cone
breaks convexity. We can still ask for the **simplest non-self-intersecting
tiling** with a given vertex set, but this requires a specific construction
(Delaunay, Epstein-Penner, etc.) computed case by case.

| Structure | Vertex set | Vertex fig. | Convex counterpart |
|-----------|------------|-------------|-------------------|
| SE Type 2 of {7,3} | Both sheets, two-sheeted hyp. | {3} (finite) | Not convex (faces cross light cone) |
| {7, ∞} from stellation | One-sheeted hyp. | ∞ (non-convex) | Open (Delaunay triangulation?) |

For the stellation-derived vertices on the one-sheeted hyperboloid, the
Delaunay triangulation would likely yield a {3, q} tiling for some q.
Computing this is an open problem.

---

## 8. Summary: What's New

### What the StackExchange provides
- 2-to-1 correspondence (R₀ vs complement, sheet-preserving vs sheet-swapping)
- Framework: hyperbolic ↔ Minkowski via Coxeter mirrors
- Examples: {∞,3}, {5,4}, {3,∞} on two-sheeted hyperboloid

### What this repository adds
1. **Vertices on the one-sheeted hyperboloid** via stellation — a mechanism
   absent from the SE answer, which only swaps sheets of the two-sheeted
   hyperboloid.
2. **Stellation depth gradation**: k=1→2→3 traces a path from ordinary
   (timelike) through ideal (lightlike) to ultraideal (spacelike) vertices.
3. **Timelike mirror discovery**: The Coxeter mirror that becomes timelike
   during stellation is essential for reconstructing the full tiling.
4. **Infinite vertex figures** from non-compact stabilizers SO(1,1) at
   spacelike vertices — natural from the Minkowski viewpoint.
5. **Edges crossing the light cone**: the {7/3} tiling is intrinsically a
   Minkowski object, not reducible to hyperbolic geometry.
6. **Category C** in the landscape: deep stellations {p/k, ∞} for k ≥ 3,
   excluded by the Polytope Wiki as "dense" but well-defined in Minkowski.

---

## 9. References

### Primary
- [StackExchange: "Regular polytopes in Minkowski spacetime"](https://math.stackexchange.com/questions/1496720/regular-polytopes-in-minkowski-spacetime)
- [Polytope Wiki: Regular polyhedron (star tiling classification)](https://polytope.miraheze.org/wiki/Regular_polyhedron)
- [Polytope Wiki: Dense polytope](https://polytope.miraheze.org/wiki/Dense_polytope)
- [Polytope Wiki: Great heptagonal tiling {7,7/2}](https://polytope.miraheze.org/wiki/Great_heptagonal_tiling)
- [Polytope Wiki: Stellated heptagonal tiling {7/2,7}](https://polytope.miraheze.org/wiki/Stellated_heptagonal_tiling)

### Repository code
- `stellate.wls` — Multi-level stellation (k=1,2,3) with convex hull
- `great_heptagons.wls` — Vertex figure analysis
- `great_heptagons_mirrors.wls` — Mirror causal type analysis
- `great_heptagons_vfch.wls` — Vertex figure convex hull
- `great_heptagons_vfch_mirror_analysis.wls` — Seed-face reconstruction

### Additional
- [Nan Ma: {7/2,7}](http://nan.ma/star/#/hyperbolic_72_7) and [{7,7/2}](http://nan.ma/star/#/hyperbolic_7_72)
- [Barrett, "Minkowski Space-Time and Hyperbolic Geometry"](https://www.eprints.soton.ac.uk/397637/2/J_F_Barrett_MICOM_2015_2018_revision_.pdf)
- McMullen & Schulte, *Abstract Regular Polytopes* (2002)
- [Schulte, "Symmetry of Polytopes and Polyhedra"](https://www.csun.edu/~ctoth/Handbook/chap18.pdf)
