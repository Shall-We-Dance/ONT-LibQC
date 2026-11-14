"""Recover existing deepTools SEM bounds from retained numeric SVG geometry.

This only supports the two linear mean-to-bound polygons produced by
``plotProfile --plotType se``.  Every mean vertex is checked against the numeric
profile table; unsupported or inconsistent figures simply have no recovered
band.  No SVG markup, styling or executable content is returned to the report.
"""
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET


_NS = '{http://www.w3.org/2000/svg}'
_NUMBER = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
_PAIR = rf'({_NUMBER})[\s,]+({_NUMBER})'
_POLYGON = re.compile(rf'\s*M\s*{_PAIR}(?:\s*L\s*{_PAIR})+\s*[zZ]\s*')
_VERTEX = re.compile(rf'[ML]\s*{_PAIR}')


def _close(a, b, tolerance=3e-6):
    return abs(a - b) <= tolerance


def _tick_positions(axis, prefix):
    ticks = []
    for group in axis.findall(_NS + 'g'):
        if not group.get('id', '').startswith(prefix):
            # An offset/multiplier makes the displayed tick numbers incomplete.
            # The ordinary vertical axis title is harmless.
            for text in group.iter(_NS + 'text'):
                if text.text and not text.get('transform', '').startswith('rotate(-90 '):
                    raise ValueError('Unsupported axis offset')
            continue
        labels = list(group.iter(_NS + 'text'))
        marks = list(group.iter(_NS + 'use'))
        if len(labels) != 1 or len(marks) != 1 or len(labels[0]):
            raise ValueError('Unsupported tick representation')
        ticks.append((labels[0].text or '', float(marks[0].get('x')), float(marks[0].get('y'))))
    if len(ticks) < 2:
        raise ValueError('Insufficient axis ticks')
    return ticks


def _polygon(path):
    d = path.get('d', '')
    if not _POLYGON.fullmatch(d):
        raise ValueError('Unsupported polygon commands')
    points = [(float(x), float(y)) for x, y in _VERTEX.findall(d)]
    if any(not math.isfinite(v) for point in points for v in point):
        raise ValueError('Nonfinite polygon coordinates')
    if len(points) < 4 or len(points) % 2:
        raise ValueError('Unsupported polygon shape')
    size = (len(points) - 2) // 2
    means, bounds = points[1:size + 1], list(reversed(points[size + 2:]))
    if points[0] != bounds[0] or points[size + 1] != bounds[-1]:
        raise ValueError('Unsupported polygon closure')
    if any(not _close(p[0], q[0]) for p, q in zip(means, bounds)):
        raise ValueError('Band and mean coordinates differ')
    if any(b[0] <= a[0] for a, b in zip(means, means[1:])):
        raise ValueError('Nonascending polygon bins')
    return means, bounds


def read_profile_band(svg_path, gene_tuple):
    """Return original ``(lower, upper)`` values, or ``None`` if unsupported.

    Missing profile bins remain NaN.  SVG coordinates are rounded by matplotlib,
    so validation permits at most a few millionths of a plotting point.  The
    returned numbers retain that source precision; no band is estimated from
    the means and no missing observations are interpolated.
    """
    try:
        if gene_tuple is None:
            return None
        labels, bins, means = gene_tuple
        count = len(bins)
        if count < 2 or len(means) != count or len(labels) != count:
            return None
        if any(not math.isfinite(x) for x in bins) or any(b <= a for a, b in zip(bins, bins[1:])):
            return None
        if any(math.isinf(v) for v in means) or not any(math.isfinite(v) for v in means):
            return None
        path = Path(svg_path)
        if path.stat().st_size > 2_000_000:
            return None
        source = path.read_bytes()
        if b'<!ENTITY' in source.upper():
            return None
        root = ET.fromstring(source)
        if root.tag != _NS + 'svg':
            return None
        # Geometry transforms would require a different coordinate conversion.
        if any(node.get('transform') for node in root.iter() if node.tag != _NS + 'text'):
            return None
        axes = [g for g in root.iter(_NS + 'g') if g.get('id', '').startswith('axes_')]
        if len(axes) != 1:
            return None
        axis = axes[0]
        x_axis = axis.find(f"{_NS}g[@id='matplotlib.axis_1']")
        y_axis = axis.find(f"{_NS}g[@id='matplotlib.axis_2']")
        if x_axis is None or y_axis is None:
            return None
        x_ticks, y_ticks = _tick_positions(x_axis, 'xtick_'), _tick_positions(y_axis, 'ytick_')
        left, right = x_ticks[0][1], x_ticks[-1][1]
        if not all(math.isfinite(v) for t in x_ticks + y_ticks for v in t[1:]) or right <= left:
            return None
        values = [(float(label.replace('\N{MINUS SIGN}', '-')), y) for label, _, y in y_ticks]
        (first_value, first_y), (last_value, last_y) = values[0], values[-1]
        slope = (last_y - first_y) / (last_value - first_value)
        intercept = first_y - slope * first_value
        if not math.isfinite(slope) or slope >= 0 or any(not math.isfinite(v) or not _close(y, intercept + slope * v) for v, y in values):
            return None
        expected_x = [left + (b - bins[0]) / (bins[-1] - bins[0]) * (right - left) for b in bins]
        groups = [g for g in axis.findall(_NS + 'g') if g.get('id', '').startswith('FillBetweenPolyCollection_')]
        if len(groups) != 2:
            return None
        recovered = []
        for group in groups:
            band = [math.nan] * count
            paths = group.findall(_NS + 'path')
            if not paths or len(paths) != len(group):
                return None
            for polygon in paths:
                vertices, bounds = _polygon(polygon)
                previous = None
                for (x, y), (_, bound_y) in zip(vertices, bounds):
                    # deepTools uses equally spaced bins in the SVG.  Validate
                    # against the entire retained axis before trusting a bin.
                    index = min(range(count), key=lambda i: abs(expected_x[i] - x))
                    if not _close(expected_x[index], x) or not math.isfinite(means[index]) or math.isfinite(band[index]):
                        return None
                    if previous is not None and index != previous + 1:
                        return None
                    if not _close(y, intercept + slope * means[index]):
                        return None
                    band[index] = (bound_y - intercept) / slope
                    previous = index
            if any(math.isfinite(mean) != math.isfinite(value) for mean, value in zip(means, band)):
                return None
            recovered.append(band)
        lower, upper = [], []
        tolerance = 4e-6 / abs(slope)
        for mean, a, b in zip(means, *recovered):
            if not math.isfinite(mean):
                lower.append(math.nan)
                upper.append(math.nan)
                continue
            lo, hi = min(a, b), max(a, b)
            if lo > mean + tolerance or hi < mean - tolerance or abs((lo + hi) / 2 - mean) > tolerance:
                return None
            lower.append(lo)
            upper.append(hi)
        return lower, upper
    except (OSError, ValueError, TypeError, ZeroDivisionError, OverflowError, ET.ParseError):
        return None
