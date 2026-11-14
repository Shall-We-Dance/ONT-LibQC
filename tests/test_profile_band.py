"""Numeric recovery of retained SEM polygons must never synthesize a band."""
import math
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow' / 'scripts'))
from report_profile_band import read_profile_band


def fixture(means=(.2, .4, .5, .7), bounds=((.1, .2, .3, .4), (.3, .6, .7, 1.0))):
    def polygon(indices, band):
        points = [(10 + 10 * i, 200 - 100 * means[i]) for i in indices]
        edge = [(10 + 10 * i, 200 - 100 * band[i]) for i in indices]
        vertices = [edge[0]] + points + [edge[-1]] + list(reversed(edge))
        d = ' '.join(('M' if i == 0 else 'L') + f' {x} {y}' for i, (x, y) in enumerate(vertices))
        return f'<path d="{d} z"/>'
    runs, run = [], []
    for i, mean in enumerate(means):
        if math.isfinite(mean):
            run.append(i)
        elif run:
            runs.append(run)
            run = []
    if run:
        runs.append(run)
    groups = ''.join(f'<g id="FillBetweenPolyCollection_{i + 1}">' + ''.join(polygon(r, band) for r in runs) + '</g>' for i, band in enumerate(bounds))
    x_axis = '<g id="matplotlib.axis_1"><g id="xtick_1"><use x="10" y="200"/><text>start</text></g><g id="xtick_2"><use x="40" y="200"/><text>end</text></g></g>'
    y_axis = '<g id="matplotlib.axis_2"><g id="ytick_1"><use x="10" y="200"/><text>0</text></g><g id="ytick_2"><use x="10" y="100"/><text>1</text></g></g>'
    svg = '<svg xmlns="http://www.w3.org/2000/svg"><g id="axes_1">' + groups + x_axis + y_axis + '</g></svg>'
    return svg, (['start', '', '', 'end'], [1, 2, 3, 4], list(means))


class ProfileBandTests(unittest.TestCase):
    def recover(self, svg, gene):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'profile.svg'
            path.write_text(svg)
            return read_profile_band(path, gene)

    def test_recovers_original_bounds_on_numeric_axis(self):
        svg, gene = fixture()
        lower, upper = self.recover(svg, gene)
        for actual, expected in zip(lower, [.1, .2, .3, .4]):
            self.assertAlmostEqual(actual, expected)
        for actual, expected in zip(upper, [.3, .6, .7, 1.]):
            self.assertAlmostEqual(actual, expected)

    def test_missing_bins_remain_missing(self):
        svg, gene = fixture(means=(.2, math.nan, .5, .7))
        lower, upper = self.recover(svg, gene)
        self.assertTrue(math.isnan(lower[1]))
        self.assertTrue(math.isnan(upper[1]))
        self.assertAlmostEqual(lower[2], .3)
        self.assertAlmostEqual(upper[3], 1.)

    def test_mismatched_means_reject_the_entire_band(self):
        svg, gene = fixture()
        gene[2][1] += .01
        self.assertIsNone(self.recover(svg, gene))

    def test_asymmetric_or_missing_bounds_are_not_sem(self):
        svg, gene = fixture(bounds=((.1, .2, .3, .4), (.3, .6, .8, 1.)))
        self.assertIsNone(self.recover(svg, gene))
        svg, gene = fixture()
        self.assertIsNone(self.recover(svg.replace('FillBetweenPolyCollection_2', 'other'), gene))

    def test_unsupported_geometry_and_axis_scaling_fail_closed(self):
        svg, gene = fixture()
        variants = [
            svg.replace('<g id="axes_1">', '<g id="axes_1" transform="scale(2)">'),
            svg.replace('<path d="M', '<path d="Q', 1),
            svg.replace('<g id="matplotlib.axis_2">', '<g id="matplotlib.axis_2"><g><text>1e6</text></g>'),
            svg.replace('<text>1</text>', '<text><tspan>1</tspan></text>'),
        ]
        for variant in variants:
            with self.subTest(variant=variant[:70]):
                self.assertIsNone(self.recover(variant, gene))

    def test_invalid_and_entity_documents_fail_closed(self):
        svg, gene = fixture()
        self.assertIsNone(self.recover('<svg>', gene))
        self.assertIsNone(self.recover('<!DOCTYPE svg [<!ENTITY number "1">]>' + svg, gene))
        self.assertIsNone(read_profile_band('/nonexistent/profile.svg', gene))


if __name__ == '__main__':
    unittest.main()
