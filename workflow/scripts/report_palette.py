"""Categorical colors shared by the HTML, SVG and printed reports."""
from zlib import crc32

# ggsci pal_npg("nrc"), in Nature Reviews Cancer palette order.
# https://nanx.me/ggsci/reference/pal_npg.html
NPG = (
    '#e64b35', '#4dbbd5', '#00a087', '#3c5488', '#f39b7f',
    '#8491b4', '#91d1c2', '#dc0000', '#7e6148', '#b09c85',
)
POOLED = '#2a6cbd'


def scope_color(scope):
    if scope == 'total':
        return POOLED
    suffix = scope.removeprefix('lib')
    if suffix != scope and suffix.isascii() and suffix.isdecimal() and int(suffix) > 0:
        return NPG[(int(suffix) - 1) % len(NPG)]
    return NPG[crc32(scope.encode()) % len(NPG)]


def bar_ink(fill):
    """Choose the higher-contrast black or white label using sRGB luminance."""
    rgb = (int(fill[i:i + 2], 16) / 255 for i in (1, 3, 5))
    linear = (v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4 for v in rgb)
    luminance = sum(w * v for w, v in zip((0.2126, 0.7152, 0.0722), linear))
    black_contrast = (luminance + 0.05) / 0.05
    white_contrast = 1.05 / (luminance + 0.05)
    return '#000000' if black_contrast >= white_contrast else '#ffffff'
