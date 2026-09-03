#!/usr/bin/env python3
"""Which factors co-move with a dependent variable over time, and when.

Each entity -- a country, a patient, a cell line -- is a trajectory through a
space whose axes are measured indicators. For a factor F and a dependent D the
signed contribution at time t is

    c(t) = (D_t - b_D) * (F_t - b_F)

with baselines b chosen by ``baseline_mode``. Positive means the two deviate
from their baselines in the same direction at that moment; negative means they
diverge. Summed over time it gives one number per (factor, dependent, entity).

A permutation test then asks a deliberately narrow question: is *this* entity's
pairing special compared with the other entities in the panel? It answers it by
swapping in another entity's dependent trajectory, keeping time order intact.

The data here is SYNTHETIC and generated below, because the repository ships no
time-series data. It is built so the answer is known: in entity 1 the dependent
really is driven by factor 1, every other entity is unrelated noise. That makes
the output checkable -- see how the planted pair separates from the rest.

Run from the repository root, after ``./build.sh``::

    PYTHONPATH=python python3 python/examples/trajectory_contributions.py
"""

import argparse
import sys

import numpy as np

from tensor_omics import (compute_all_contributions, compute_contributions,
                          compute_p_values, compute_velocity_acceleration_contributions,
                          normalize_all_trajectories, perform_permutation_test)

#: index (1-based) of the dependent variable in the synthetic panel
DEPENDENT = 4
#: index (1-based) of the factor that genuinely drives it, in entity 1 only
DRIVER = 1


def make_panel(n_factors, n_entities, n_timepoints, seed):
    """Build the synthetic panel described in the module docstring.

    Returns
    -------
    np.ndarray[np.float64] of shape (n_factors, n_entities, n_timepoints), order='F'
        The layout every routine in this module expects. Note the axis order:
        the time axis is *last*.
    """
    rng = np.random.default_rng(seed)
    trajectories = np.zeros((n_factors, n_entities, n_timepoints), order="F")

    # every entity, every indicator: an unrelated series
    for i_entity in range(n_entities):
        for i_factor in range(n_factors):
            trajectories[i_factor, i_entity, :] = rng.standard_normal(n_timepoints)

    # Entity 1 alone: its dependent follows its own driving factor. Note what is
    # NOT done here -- the dependent keeps the same distribution as every other
    # entity's, so only the *pairing* is special. Giving entity 1 a dependent that
    # also looks different from the rest would make every factor in it significant,
    # because the test compares against other entities' dependents.
    trajectories[DEPENDENT - 1, 0, :] = (
        trajectories[DRIVER - 1, 0, :] + 0.05 * rng.standard_normal(n_timepoints)
    )
    return trajectories


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--baseline", default="baseline_mean",
                        choices=["baseline_raw", "baseline_mean", "baseline_min"],
                        help="where deviations are measured from (default: baseline_mean)")
    parser.add_argument("--permutations", type=int, default=2000,
                        help="permutations per tested pair (default: 2000)")
    parser.add_argument("--normalize", action="store_true",
                        help="min-max each series to [0,1] before analysing")
    parser.add_argument("--entities", type=int, default=120)
    parser.add_argument("--timepoints", type=int, default=21)
    parser.add_argument("--seed", type=int, default=11)
    args = parser.parse_args(argv)

    n_factors = 4
    trajectories = make_panel(n_factors, args.entities, args.timepoints, args.seed)
    print(f"panel: {n_factors} indicators x {args.entities} entities "
          f"x {args.timepoints} timepoints, baseline {args.baseline}")

    if args.normalize:
        # Each series is scaled independently, so shape is preserved and only
        # magnitude is removed. `status` is per (factor, entity) and is NOT
        # raised: a constant series cannot be scaled and comes back all zeros.
        result = normalize_all_trajectories(trajectories)
        n_flat = int(np.count_nonzero(np.asarray(result["status"])))
        trajectories = np.asfortranarray(result["trajectories_norm"])
        print(f"normalized to [0,1]; {n_flat} constant series could not be scaled")

    # --- descriptive: one total per (factor, dependent, entity) -------------
    factors = np.array([f for f in range(1, n_factors + 1) if f != DEPENDENT], dtype=np.int32)
    contributions = compute_all_contributions(
        trajectories, factors, np.array([DEPENDENT], dtype=np.int32), args.baseline
    )
    totals = contributions["total_contributions"]      # (n_factors, n_dependents, n_entities)
    print(f"\ntotal contribution to indicator {DEPENDENT}, entity 1:")
    for row, i_factor in enumerate(factors):
        mark = "  <- planted driver" if i_factor == DRIVER else ""
        print(f"  factor {i_factor}: {totals[row, 0, 0]:8.3f}{mark}")

    # --- inferential: is entity 1's pairing special in this panel? ----------
    print(f"\npermutation test, entity 1, {args.permutations} permutations:")
    for i_factor in factors:
        observed = compute_contributions(
            trajectories[i_factor - 1, 0, :].copy(),
            trajectories[DEPENDENT - 1, 0, :].copy(),
            args.baseline,
        )
        permuted = perform_permutation_test(
            trajectories, int(i_factor), DEPENDENT, 1,
            args.baseline, args.permutations, random_seed=5,
        )
        p = compute_p_values(
            observed["local_contributions"], observed["total_contribution"],
            permuted["local_contributions"], permuted["total_contributions"],
        )
        # p counts (permuted >= observed) / n_permutations, with no correction,
        # so it can be exactly 0. Report that as "< 1/n" rather than as zero.
        total_p = p["total_p_value"]
        shown = f"{total_p:.4f}" if total_p > 0 else f"<{1.0 / args.permutations:.4g}"
        n_local = int(np.count_nonzero(np.asarray(p["local_p_values"]) < 0.05))
        print(f"  factor {i_factor}: C={observed['total_contribution']:8.3f}  "
              f"p={shown}  ({n_local} of {args.timepoints} timepoints below 0.05)")

    # --- the derivatives, without transposing anything by hand -------------
    # compute_velocity_trajectories returns (n_timepoints-1, n_factors, n_entities),
    # which is NOT the layout the contribution routines take. This one starts from
    # positions and stays in the original layout.
    derivatives = compute_velocity_acceleration_contributions(trajectories, args.baseline)
    vel = derivatives["contrib_velocity"][DRIVER - 1, DEPENDENT - 1, 0]
    acc = derivatives["contrib_acceleration"][DRIVER - 1, DEPENDENT - 1, 0]
    print(f"\nentity 1, factor {DRIVER} vs indicator {DEPENDENT}: "
          f"velocity contribution {vel:.3f}, acceleration contribution {acc:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
