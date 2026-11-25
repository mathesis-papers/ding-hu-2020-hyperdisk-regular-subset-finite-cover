# Negative answers to Ding–Hu's two open finite-cover questions

Haobo Ma — ChronoAI Pte Ltd, Singapore  
Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore

The note answers both open questions of Ding and Hu, Section V, p. 8,
arXiv:1909.11256v3, with a regular subset of a four-dimensional hyperdisk
whose contained hyperdisks are all rays and whose least cover cardinality
is continuum. It gives complete maskable sets with the same cover
cardinality for every positive marginal spectrum of rank d >= 4, at every
input dimension 3 <= n < d. These are exactly the dimension pairs admitting
nondegenerate counterexamples. The span-dimension bound proposed in the 2019
version and both degenerate examples also have explicit coverage.

For Zhu's distinct convex-hull maximality conjecture (arXiv:2010.07843v2,
Section VII, p. 9, Conjecture 1), the note proves the assertion for d=2,3
and gives a rank-two counterexample for every d >= 4. The full
fixed-diagonal density-matrix set is the corrected maximal maskable set.

## Files and build

- `note.tex`: self-contained LaTeX, including the bibliography.
- `note.pdf`: readable mathematical note.
- `metadata.json`: bibliographic description, license, problem identifier and source link.
- `verify.py`: standard-library Python 3 exact finite checks.
- `checker-output.json`: the checker's mathematical output.

Build with Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, run `pdflatex note.tex` twice. The bibliography is embedded;
BibTeX and external figures are unnecessary. Run the finite checks with:

```sh
python3 verify.py
```

## Mathematical coverage

The support obstruction and all-parameter construction have complete written
proofs in the note. The finite checker enumerates the 15 support partitions,
checks identities at 32 exact rational parameters, proves independence of the
three displayed span vectors by a nonzero exact minor, and checks the marginal
spectrum. It also checks the spectrum extension in dimensions 4 through 10,
17 exact parameter values for each degenerate example, and the four-row
rank-two correlation witness and its Hermitian constraint determinant.
The finite examples do not replace the all-parameter and all-dimension proofs.

The declaration `result` in
[`HyperdiskRegularSubsetCover.lean`](https://github.com/the-omega-institute/trureturing/blob/01b573fd9b691a4be3a8b0ecbb41b3311eda19b5/D5/S3/Quantum/Masking/HyperdiskRegularSubsetCover.lean)
proves the negation of the universal finite-cover assertion for the formal
hyperdisk and regular-subset definitions. Its local lemmas establish
phase membership on the interval, a countable family of distinct rays and
the single-ray obstruction. The continuum cardinality, span equality,
minimal dimensions, arbitrary-spectrum masking realization, degenerate
examples and density-matrix maximality have written proofs in the note.

The posing source is Feng Ding and Xueyuan Hu, *Masking quantum information
on hyperdisks*, Phys. Rev. A **102**, 042404 (2020),
[arXiv:1909.11256v3](https://arxiv.org/abs/1909.11256v3), Section V, p. 8.

## License

CC BY 4.0 (Creative Commons Attribution 4.0 International).
