import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/scripts"))
from custom_regions_ui import render_section, upgrade_html
from report_palette import NPG, bar_ink, scope_color
from report_refresh import refresh_series_paints
from report_sublibrary import mixture_counts, render_sublibrary_composition
from species_composition import AMBIG, NONE, NAMES, SPECIES_COLORS, render


class PaletteTests(unittest.TestCase):
    def test_sublibraries_keep_their_color_across_gaps_and_display_names(self):
        metric = dict(library_units=20, sub_libraries=[
            dict(scope_id='lib03', scope_name='Control', library_units=10),
            dict(scope_id='lib01', scope_name='Treatment', library_units=6),
        ])
        total, rows = mixture_counts(metric)
        self.assertEqual(total, 20)
        self.assertEqual(rows[:2], [('Control', 10, '#00a087'), ('Treatment', 6, '#e64b35')])
        self.assertEqual(scope_color('lib10'), '#b09c85')
        self.assertEqual(scope_color('lib11'), scope_color('lib01'))
        doc = render_sublibrary_composition(metric)
        for label, _, color in rows[:2]:
            self.assertIn(f'background:{color};color:#ffffff" title="{label}:', doc)
        self.assertIn('background:#8794a7;color:#000000" title="Unassigned / ambiguous:', doc)
        custom = render_section(dict(reference='hg38', scopes={'lib03': {'name': 'Control'}}))
        self.assertIn("value='lib03' data-color='#00a087'>Control</option>", custom)
        self.assertEqual(upgrade_html(custom), custom)

    def test_sublibrary_stacks_use_both_label_colors_across_the_npg_palette(self):
        metric = dict(library_units=100, sub_libraries=[
            dict(scope_id=f'lib{i:02d}', library_units=10) for i in range(1, 11)])
        doc = render_sublibrary_composition(metric)
        for i, fill in enumerate(NPG, 1):
            ink = '#ffffff' if i in (1, 2, 3, 4, 8, 9) else '#000000'
            self.assertIn(f'background:{fill};color:{ink}" title="lib{i:02d}: 10.0%', doc)

    def test_species_legend_and_stacks_use_reference_ids(self):
        for panel, expected in [(['mm39', 'hg38'], 'hg38'), (['hg38', 'mm39'], 'mm39'), (['hg38'], 'hg38')]:
            keys = panel + [AMBIG, NONE]
            doc = render(dict(panel=panel, reference=expected, rows=[
                dict(scope='lib01', n=len(keys), counts=dict.fromkeys(keys, 1), windows={})]))
            for key in keys:
                label = NAMES.get(key, key) + (' (expected)' if key == expected else '')
                color = SPECIES_COLORS[key]
                self.assertIn(f"background:{color}'></i>{label}</li>", doc)
                self.assertIn(f"background:{color};color:{bar_ink(color)}'", doc)
        bmj = {'#2a6ebb', '#f0ab00', '#c50084', '#7d5cc6', '#e37222',
               '#69be28', '#00b2a9', '#cd202c', '#747678'}
        self.assertEqual({SPECIES_COLORS[k] for k in NAMES}, bmj)

    def test_automatic_label_colors_have_normal_text_contrast_on_all_fills(self):
        for fill in (*NPG, *SPECIES_COLORS.values(), '#717d92', '#aab5c4', '#8794a7'):
            rgb = [int(fill[j:j + 2], 16) / 255 for j in (1, 3, 5)]
            rgb = [v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4 for v in rgb]
            luminance = sum(v * w for v, w in zip(rgb, (0.2126, 0.7152, 0.0722)))
            ink = {'#000000': 0, '#ffffff': 1}[bar_ink(fill)]
            ratio = (max(ink, luminance) + 0.05) / (min(ink, luminance) + 0.05)
            self.assertGreaterEqual(ratio, 4.5, fill)

    def test_white_is_preferred_when_both_colors_are_readable(self):
        self.assertEqual(bar_ink(SPECIES_COLORS['sacCer3']), '#ffffff')
        self.assertEqual(bar_ink('#767676'), '#ffffff')
        self.assertEqual(bar_ink('#777777'), '#000000')

    def test_retained_scoped_figures_change_only_series_paint(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            assets = root / 'assets'
            assets.mkdir()
            original = """<svg><title>#000080</title><path d="M0 1 L6 7" style="stroke: #000080; fill: none"/></svg>"""
            scoped = assets / 's.lib02.TSS_TES.svg'
            pooled = assets / 's.total.TSS_TES.svg'
            rna = root / 'data/rna/s/scope_0123456789ab.profiles.svg'
            rna.parent.mkdir(parents=True)
            for path in (scoped, pooled, rna):
                path.write_text(original)
            for _ in range(2):
                refresh_series_paints(root, 's', 'lib02', rna.name)
                self.assertEqual(pooled.read_text(), original)
                for path in (scoped, rna):
                    self.assertEqual(path.read_text(), original.replace('stroke: #000080', 'stroke: #4dbbd5'))
