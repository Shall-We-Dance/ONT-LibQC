"""Categorical colors shared by the HTML, SVG and printed reports."""
from zlib import crc32

SUBLIBRARY_COLORS = ('#2878c8', '#df7a14', '#219357', '#9356c0', '#ca4d87', '#1196aa', '#819610', '#d34d45')
POOLED = '#496576'
AMBIGUOUS_COLOR = '#a4a49c'
UNASSIGNED_COLOR = '#dddcd4'


def scope_color(scope, index=None):
    """Resolve by stable library ID, including calls that also pass its former list index."""
    if scope == 'total':
        return POOLED
    suffix = scope.removeprefix('lib')
    if suffix != scope and suffix.isascii() and suffix.isdecimal() and int(suffix) > 0:
        return SUBLIBRARY_COLORS[(int(suffix) - 1) % len(SUBLIBRARY_COLORS)]
    return SUBLIBRARY_COLORS[crc32(scope.encode()) % len(SUBLIBRARY_COLORS)]


def bar_ink(fill):
    """Prefer white when it meets the 4.5:1 contrast ratio for small percentage labels."""
    rgb = (int(fill[i:i + 2], 16) / 255 for i in (1, 3, 5))
    linear = (v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4 for v in rgb)
    luminance = sum(w * v for w, v in zip((0.2126, 0.7152, 0.0722), linear))
    white_contrast = 1.05 / (luminance + 0.05)
    return '#ffffff' if white_contrast >= 4.5 else '#000000'
