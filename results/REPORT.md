# PTP1B pose-ranking results

Retrospective known-pocket holo redocking; no MD performed.

| Fragment | Method | Top-1 RMSD (Å) | Top-1 success | Top-3 success | Oracle RMSD (Å) |
|---|---|---:|---:|---:|---:|
| 8G65 | Vina | 0.232 | 1 | 1 | 0.232 |
| 8G65 | ML | 0.232 | 1 | 1 | 0.232 |
| 8G65 | Shuffled-label ML | 3.852 | 0 | 0 | 0.232 |
| 8G68 | Vina | 0.478 | 1 | 1 | 0.478 |
| 8G68 | ML | 0.478 | 1 | 1 | 0.478 |
| 8G68 | Shuffled-label ML | 6.626 | 0 | 1 | 0.478 |
| 8G69 | Vina | 0.515 | 1 | 1 | 0.515 |
| 8G69 | ML | 0.515 | 1 | 1 | 0.515 |
| 8G69 | Shuffled-label ML | 6.237 | 0 | 0 | 0.515 |
| 8G67 | Vina | 5.814 | 0 | 0 | 1.175 |
| 8G67 | ML | 4.318 | 0 | 1 | 1.175 |
| 8G67 | Shuffled-label ML | 4.769 | 0 | 0 | 1.175 |
| 8G6A | Vina | 0.933 | 1 | 1 | 0.689 |
| 8G6A | ML | 0.933 | 1 | 1 | 0.689 |
| 8G6A | Shuffled-label ML | 6.603 | 0 | 0 | 0.689 |

Mean paired ML − Vina top-1 RMSD: -0.299 Å. Negative favors ML.
Descriptive fragment-bootstrap 95% interval: [-0.8978184035392152, 0.0]. Two independent chemical series only.

![Ranking](ranking.png)

![Pose overlays](poses.png)

## Theory and interpretation

RMSD is the square root of the mean squared heavy-atom displacement, minimized over graph automorphisms.
Ligands are never independently aligned to the crystal ligand: that would erase placement errors.
Vina samples conformations inside a known pocket and scores approximate interactions. Its scores are not measured binding affinities.
A standardized ridge model (alpha=10) predicts log(1+RMSD) from Vina score and receptor-contact features.
No native coordinates, native contacts, atom displacement, pocket ID, or fragment identity enter its feature vector.
The entire held-out chemical series is excluded from training. The two folds test transfer between different allosteric pockets as well as chemistry.
The candidate oracle separates sampling failure from ranking failure; shuffled training labels provide a negative control.
Polar proximity is a coarse contact descriptor, not an angle-aware hydrogen-bond calculation.

## Limits

Five structures are insufficient to establish generalizable AI performance. Poses are not independent experiments.
Holo receptors and crystal-derived box centers make this an optimistic retrospective benchmark; it does not discover cryptic pockets.
One selected monomer is used. Alternate conformations, crystal packing, missing residues, protonation, waters, and receptor flexibility can affect results.
Inspect data-quality.json and the raw preparation logs before interpreting a result. Default ligand protonation follows the deposited CCD graph.
This workflow does not reproduce the Anton study, run MD, establish affinity, or identify a new drug.

Paper: https://doi.org/10.1021/acs.jcim.3c00236
