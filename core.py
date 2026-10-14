"""Geometry and grouped evaluation. Reference coordinates never enter features."""
import numpy as np
from scipy.spatial.distance import cdist
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

FEATURES = ['vina_score', 'clash', 'contact', 'polar_contact', 'hydrophobic_contact',
            'nearest_distance', 'burial_fraction']
RADII = {'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'P': 1.80,
         'F': 1.47, 'CL': 1.75, 'BR': 1.85, 'I': 1.98}


def validate_xyz(x):
    x = np.asarray(x, dtype=float)
    if x.ndim != 2 or x.shape[1] != 3 or not len(x) or not np.isfinite(x).all():
        raise ValueError('Expected finite, nonempty N x 3 coordinates')
    return x


def fixed_rmsd(a, b):
    a, b = validate_xyz(a), validate_xyz(b)
    if a.shape != b.shape:
        raise ValueError('Atom counts differ')
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def contact_features(ligand, ligand_elements, protein, protein_elements, score):
    ligand, protein = validate_xyz(ligand), validate_xyz(protein)
    if len(ligand_elements) != len(ligand) or len(protein_elements) != len(protein):
        raise ValueError('Element counts differ from coordinates')
    if not np.isfinite(score):
        raise ValueError('Nonfinite docking score')
    le, pe = np.asarray(ligand_elements), np.asarray(protein_elements)
    radii = np.array([RADII.get(e.upper(), 1.7) for e in le])[:, None]
    radii = radii + np.array([RADII.get(e.upper(), 1.7) for e in pe])[None, :]
    d = cdist(ligand, protein)
    gap = d - radii
    n = len(ligand)
    polar = np.isin(le, ['N', 'O', 'S'])[:, None] & np.isin(pe, ['N', 'O', 'S'])[None, :]
    hydro = (le == 'C')[:, None] & (pe == 'C')[None, :]
    return np.array([score, np.sum(np.maximum(-gap, 0) ** 2) / n,
                     np.sum(np.exp(-(gap / 0.8) ** 2)) / n,
                     np.sum(polar * ((d >= 2.5) & (d <= 3.5))) / n,
                     np.sum(hydro * ((d >= 3) & (d <= 4.5))) / n,
                     np.mean(d.min(axis=1)), np.mean((d < 4.5).sum(axis=1) >= 4)])


def grouped_predictions(x, rmsd, groups, seed=17, shuffle=False):
    """Outer leave-series-out; fixed hyperparameters, equal fragment weights."""
    x, rmsd, groups = np.asarray(x), np.asarray(rmsd), np.asarray(groups)
    if not np.isfinite(x).all() or not np.isfinite(rmsd).all():
        raise ValueError('Nonfinite training data')
    series = np.array([g.split(':')[0] for g in groups])
    if len(set(series)) < 2:
        raise ValueError('At least two independent series are required')
    predictions = np.full(len(x), np.nan)
    rng = np.random.default_rng(seed)
    folds = []
    for held in sorted(set(series)):
        train, test = series != held, series == held
        y = np.log1p(rmsd[train])
        training_groups = groups[train]
        if shuffle:
            # Negative control: preserve each fragment's label distribution.
            for g in sorted(set(training_groups)):
                sel = training_groups == g
                y[sel] = rng.permutation(y[sel])
        weights = np.array([1 / np.sum(training_groups == g) for g in training_groups])
        weights *= len(weights) / weights.sum()
        model = make_pipeline(StandardScaler(), Ridge(alpha=10))
        model.fit(x[train], y, ridge__sample_weight=weights)
        predictions[test] = model.predict(x[test])
        folds.append({'held_out_series': held, 'train_fragments': sorted(set(groups[train])),
                      'test_fragments': sorted(set(groups[test]))})
    return predictions, folds


def rank_metrics(rmsd, scores, threshold=2.0):
    rmsd, scores = np.asarray(rmsd), np.asarray(scores)
    if not len(rmsd) or len(rmsd) != len(scores) or not np.isfinite(scores).all():
        raise ValueError('Invalid ranking arrays')
    order = np.argsort(scores, kind='stable')
    return {'top1_rmsd_A': float(rmsd[order[0]]),
            'top1_success': int(rmsd[order[0]] <= threshold),
            'top3_success': int(np.any(rmsd[order[:3]] <= threshold)),
            'oracle_success': int(np.any(rmsd <= threshold)),
            'oracle_rmsd_A': float(rmsd.min())}


def paired_bootstrap(a, b, seed=17, n=10000):
    """Fragment-level paired interval, never a pose-level pseudo-replicate CI."""
    delta = np.asarray(a, float) - np.asarray(b, float)
    if not len(delta):
        raise ValueError('Empty bootstrap sample')
    rng = np.random.default_rng(seed)
    draws = delta[rng.integers(0, len(delta), (n, len(delta)))].mean(axis=1)
    return {'mean_difference': float(delta.mean()),
            'percentile_95_CI': np.quantile(draws, [.025, .975]).tolist(),
            'n_fragments': len(delta), 'resamples': n}
