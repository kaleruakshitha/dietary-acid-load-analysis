"""Reusable statistics for the two-day PRAL reliability analysis."""

from __future__ import annotations

import numpy as np


def _finite_vector(values: np.ndarray | list[float], name: str) -> np.ndarray:
    result = np.asarray(values, dtype=float)
    if result.ndim != 1 or result.size == 0:
        raise ValueError(f"{name} must be a nonempty one-dimensional vector")
    if not np.isfinite(result).all():
        raise ValueError(f"{name} must contain only finite values")
    return result


def _validate_weights(weights: np.ndarray | list[float], size: int) -> np.ndarray:
    result = _finite_vector(weights, "weights")
    if result.size != size:
        raise ValueError("weights must have the same length as the values")
    if np.any(result <= 0):
        raise ValueError("weights must be strictly positive")
    return result


def weighted_mean(
    values: np.ndarray | list[float], weights: np.ndarray | list[float]
) -> float:
    """Return a positive-weighted arithmetic mean."""

    x = _finite_vector(values, "values")
    w = _validate_weights(weights, x.size)
    return float(np.average(x, weights=w))


def weighted_standard_deviation(
    values: np.ndarray | list[float], weights: np.ndarray | list[float]
) -> float:
    """Return the descriptive population SD under normalized survey weights."""

    x = _finite_vector(values, "values")
    w = _validate_weights(weights, x.size)
    mean = np.average(x, weights=w)
    return float(np.sqrt(np.average((x - mean) ** 2, weights=w)))


def weighted_quantile(
    values: np.ndarray | list[float],
    weights: np.ndarray | list[float],
    probabilities: np.ndarray | list[float],
) -> np.ndarray:
    """Return inverse-CDF quantiles from a weighted empirical distribution."""

    x = _finite_vector(values, "values")
    w = _validate_weights(weights, x.size)
    q = np.asarray(probabilities, dtype=float)
    if q.ndim != 1 or not np.isfinite(q).all() or np.any((q < 0) | (q > 1)):
        raise ValueError("probabilities must be a finite vector between zero and one")
    order = np.argsort(x, kind="mergesort")
    sorted_x = x[order]
    cumulative = np.cumsum(w[order]) / np.sum(w)
    indices = np.searchsorted(cumulative, q, side="left")
    indices = np.clip(indices, 0, sorted_x.size - 1)
    return sorted_x[indices]


def weighted_pearson(
    first: np.ndarray | list[float],
    second: np.ndarray | list[float],
    weights: np.ndarray | list[float],
) -> float:
    """Return the normalized-weight Pearson correlation point estimate."""

    x = _finite_vector(first, "first")
    y = _finite_vector(second, "second")
    if y.size != x.size:
        raise ValueError("first and second must have the same length")
    w = _validate_weights(weights, x.size)
    x_centered = x - np.average(x, weights=w)
    y_centered = y - np.average(y, weights=w)
    covariance = np.average(x_centered * y_centered, weights=w)
    x_variance = np.average(x_centered**2, weights=w)
    y_variance = np.average(y_centered**2, weights=w)
    denominator = np.sqrt(x_variance * y_variance)
    if denominator == 0:
        raise ValueError("correlation is undefined for a constant vector")
    return float(covariance / denominator)


def icc_two_way_random_absolute(
    measurements: np.ndarray | list[list[float]],
) -> tuple[float, float]:
    """Return ICC(A,1) and ICC(C,1) for a complete subject-by-occasion matrix.

    ICC(A,1) is a two-way random-effects, absolute-agreement, single-measure
    coefficient. ICC(C,1), the consistency counterpart, is returned as a
    sensitivity statistic.
    """

    values = np.asarray(measurements, dtype=float)
    if values.ndim != 2 or values.shape[0] < 2 or values.shape[1] < 2:
        raise ValueError("measurements must contain at least two subjects and occasions")
    if not np.isfinite(values).all():
        raise ValueError("measurements must contain only finite values")

    n_subjects, n_occasions = values.shape
    grand_mean = values.mean()
    subject_means = values.mean(axis=1)
    occasion_means = values.mean(axis=0)

    ss_subject = n_occasions * np.sum((subject_means - grand_mean) ** 2)
    ss_occasion = n_subjects * np.sum((occasion_means - grand_mean) ** 2)
    residuals = (
        values
        - subject_means[:, None]
        - occasion_means[None, :]
        + grand_mean
    )
    ss_error = np.sum(residuals**2)

    ms_subject = ss_subject / (n_subjects - 1)
    ms_occasion = ss_occasion / (n_occasions - 1)
    ms_error = ss_error / ((n_subjects - 1) * (n_occasions - 1))

    absolute_denominator = (
        ms_subject
        + (n_occasions - 1) * ms_error
        + n_occasions * (ms_occasion - ms_error) / n_subjects
    )
    consistency_denominator = ms_subject + (n_occasions - 1) * ms_error
    if absolute_denominator == 0 or consistency_denominator == 0:
        raise ValueError("ICC is undefined because there is no measurable variance")

    icc_absolute = (ms_subject - ms_error) / absolute_denominator
    icc_consistency = (ms_subject - ms_error) / consistency_denominator
    return float(icc_absolute), float(icc_consistency)


def cohen_kappa_from_matrix(
    matrix: np.ndarray | list[list[float]], scheme: str = "linear"
) -> float:
    """Return Cohen's kappa from a square count or weight matrix."""

    counts = np.asarray(matrix, dtype=float)
    if counts.ndim != 2 or counts.shape[0] != counts.shape[1] or counts.shape[0] < 2:
        raise ValueError("matrix must be square with at least two categories")
    if not np.isfinite(counts).all() or np.any(counts < 0) or counts.sum() <= 0:
        raise ValueError("matrix must contain finite nonnegative values with a positive sum")
    if scheme not in {"unweighted", "linear", "quadratic"}:
        raise ValueError("scheme must be unweighted, linear, or quadratic")

    categories = counts.shape[0]
    distance = np.abs(np.subtract.outer(np.arange(categories), np.arange(categories)))
    if scheme == "unweighted":
        agreement_weights = np.eye(categories)
    elif scheme == "linear":
        agreement_weights = 1 - distance / (categories - 1)
    else:
        agreement_weights = 1 - (distance / (categories - 1)) ** 2

    observed = counts / counts.sum()
    expected = np.outer(observed.sum(axis=1), observed.sum(axis=0))
    observed_agreement = np.sum(agreement_weights * observed)
    expected_agreement = np.sum(agreement_weights * expected)
    if np.isclose(expected_agreement, 1):
        raise ValueError("kappa is undefined because expected agreement is one")
    return float((observed_agreement - expected_agreement) / (1 - expected_agreement))
