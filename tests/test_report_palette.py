import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/scripts"))
from custom_regions_ui import render_section, upgrade_html
from report_palette import AMBIGUOUS_COLOR, SUBLIBRARY_COLORS, UNASSIGNED_COLOR, bar_ink, scope_color
from report_refresh import refresh_gene_profile_color
from report_sublibrary import mixture_counts, render_sublibrary_composition
from species_composition import AMBIG, NONE, NAMES, SPECIES_COLORS, SPECIES_STATUS_COLORS, render


class PaletteTests(unittest.TestCase):
    def test_sublibraries_keep_their_color_across_gaps_and_display_names(self):
        metric = dict(library_units=20, sub_libraries=[
            dict(scope_id='lib03', scope_name='Control', library_units=10),
            dict(scope_id='lib01', scope_name='Treatment', library_units=6),
        ])
        total, rows = mixture_counts(metric)
        self.assertEqual(total, 20)
        self.assertEqual(rows[:2], [('Control', 10, '#219357'), ('Treatment', 6, '#2878c8')])
        self.assertEqual(scope_color('lib08'), '#d34d45')
        self.assertEqual(scope_color('lib09'), scope_color('lib01'))
        doc = render_sublibrary_composition(metric)
        for label, _, color in rows[:2]:
            self.assertIn(f'background:{color};color:#ffffff" title="{label}:', doc)
        self.assertIn('background:#bcbab1;color:#000000" title="Unassigned / ambiguous:', doc)
        custom = render_section(dict(reference='hg38', scopes={'lib03': {'name': 'Control'}}))
        self.assertIn("value='lib03' data-color='#219357'>Control</option>", custom)
        self.assertEqual(upgrade_html(custom), custom)
        self.assertEqual(scope_color('named-scope', 0), scope_color('named-scope', 5))

    def test_sublibrary_stacks_use_both_label_colors_across_the_palette(self):
        metric = dict(library_units=100, sub_libraries=[
            dict(scope_id=f'lib{i:02d}', library_units=10) for i in range(1, 9)])
        doc = render_sublibrary_composition(metric)
        for i, fill in enumerate(SUBLIBRARY_COLORS, 1):
            ink = '#ffffff' if i in (1, 2, 3, 4) else '#000000'
            self.assertIn(f'background:{fill};color:{ink}" title="lib{i:02d}: 10.0%', doc)

    def test_species_legend_and_stacks_use_reference_ids(self):
        for panel, expected in [(['mm39', 'hg38'], 'hg38'), (['hg38', 'mm39'], 'mm39'), (['hg38'], 'hg38')]:
            keys = panel + [AMBIG, NONE]
            doc = render(dict(panel=panel, reference=expected, rows=[
                dict(scope='lib01', n=len(keys), counts=dict.fromkeys(keys, 1), windows={})]))
            for key in keys:
                label = NAMES.get(key, key) + (' (expected)' if key == expected else '')
                color = {**SPECIES_COLORS, **SPECIES_STATUS_COLORS}[key]
                self.assertIn(f"background:{color}'></i>{label}</li>", doc)
                self.assertIn(f"background:{color};color:{bar_ink(color)}'", doc)
        self.assertEqual(SPECIES_COLORS['hg38'], '#cd202c')
        self.assertEqual(SPECIES_COLORS['dm6'], '#00b2a9')
        self.assertEqual(SPECIES_STATUS_COLORS, {AMBIG: '#a7aab8', NONE: '#e2e7ec'})
        self.assertEqual(set(SPECIES_COLORS), set(NAMES) | {'danRer11'})
        self.assertEqual(SPECIES_COLORS['danRer11'], SPECIES_COLORS['GRCz12ab'])
        for fish in ('danRer11', 'GRCz12ab'):
            panel = [*sorted(set(NAMES) - {'GRCz12ab'}), fish, 'custom']
            colors = render(dict(panel=panel, reference=fish, rows=[
                dict(scope='lib01', n=len(panel), counts=dict.fromkeys(panel, 1), windows={})]))
            self.assertNotIn(f"background:{SPECIES_COLORS[fish]}'></i>custom</li>", colors)

    def test_automatic_label_colors_have_normal_text_contrast_on_all_fills(self):
        for fill in (*SUBLIBRARY_COLORS, *SPECIES_COLORS.values(), *SPECIES_STATUS_COLORS.values(), AMBIGUOUS_COLOR, UNASSIGNED_COLOR, '#bcbab1'):
            rgb = [int(fill[j:j + 2], 16) / 255 for j in (1, 3, 5)]
            rgb = [v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4 for v in rgb]
            luminance = sum(v * w for v, w in zip(rgb, (0.2126, 0.7152, 0.0722)))
            ink = {'#000000': 0, '#ffffff': 1}[bar_ink(fill)]
            ratio = (max(ink, luminance) + 0.05) / (min(ink, luminance) + 0.05)
            self.assertGreaterEqual(ratio, 4.5, fill)

    def test_white_is_preferred_when_both_colors_are_readable(self):
        self.assertEqual(bar_ink(SPECIES_COLORS['hg38']), '#ffffff')
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
            for old in ('#000080', '#4dbbd5'):
                source = original.replace('#000080', old)
                scoped.write_text(source)
                pooled.write_text(source)
                for _ in range(2):
                    refresh_gene_profile_color(scoped, 'lib02')
                    self.assertEqual(pooled.read_text(), source)
                    self.assertEqual(scoped.read_text(), source.replace('stroke: ' + old, 'stroke: #df7a14'))
