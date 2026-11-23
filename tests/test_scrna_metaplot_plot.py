import sys
import unittest
from pathlib import Path
from xml.etree import ElementTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow' / 'scripts'))
from report_palette import POOLED, scope_color
from scrna_metaplot_plot import figures
from scrna_report import metaplot_section


def segment(genes=2):
    return {'contributing_genes': genes, 'covered_genes': genes, 'supporting_units': genes * 2,
            'mean_relative_coverage': [0., 1., 2.] if genes else [None] * 3,
            'sem_across_genes': [0., .4, .8] if genes >= 2 else [None] * 3}


def profile(genes, coding):
    return {key: segment(genes if key == 'gene_body' else coding)
            for key in ('gene_body', 'five_prime_UTR', 'CDS', 'three_prime_UTR')}


class TranscriptPlots(unittest.TestCase):
    def test_coding_and_gene_body_fallback_have_accessible_downloads(self):
        metrics = {'scopes': {'all': profile(5, 3), 'lib01': profile(1, 0), 'lib02': profile(0, 0)}}
        images = figures(metrics)
        self.assertEqual(set(images), {'all', 'lib01'})
        for scope, axis in [('all', 'CDS'), ('lib01', 'spliced transcripts')]:
            svg = ElementTree.fromstring(images[scope])
            self.assertEqual(svg.get('role'), 'img')
            self.assertIn('transcript coverage', svg.find('{*}title').text)
            self.assertIn(axis, svg.find('{*}desc').text)
            color = POOLED if scope == 'all' else scope_color(scope)
            self.assertIn(f'stroke: {color}; stroke-width: 1.7', images[scope].decode('utf-8'))
        page = metaplot_section(metrics, images, 20)
        self.assertIn('Genes with coverage', page)
        self.assertIn('transcript_metaplot_lib01.svg', page)
        self.assertIn('1 gene · 2 cDNA fragments', page)
        self.assertNotIn('lib02', page)

    def test_empty_summary_and_duplicate_single_scope(self):
        empty = {'scopes': {'all': profile(0, 0)}}
        self.assertEqual(figures(empty), {})
        self.assertIn('No eligible', metaplot_section(empty, {}, 20))
        single = {'scopes': {'all': profile(2, 2), 'lib01': profile(2, 2)}}
        page = metaplot_section(single, {'all': b'<svg/>', 'lib01': b'<svg/>'}, 20)
        self.assertEqual(page.count('<figure '), 1)
