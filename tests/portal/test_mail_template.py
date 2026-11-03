"""Notification presentation and delivery identity regressions; no mail is sent."""
import json
import unittest
from html import escape
from html.parser import HTMLParser

from ont_portal.mail_template import build_message
from ont_portal.store import token_for


class MailMarkup(HTMLParser):
    def __init__(self, html):
        super().__init__()
        self.tags = []
        self.links = []
        self.images = []
        self.headings = []
        self._heading = None
        self.feed(html)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        self.tags.append(tag)
        if tag == 'a':
            self.links.append(attrs.get('href'))
        if tag == 'img':
            self.images.append(attrs.get('src'))
        if tag == 'h1':
            self._heading = []

    def handle_data(self, data):
        if self._heading is not None:
            self._heading.append(data)

    def handle_endtag(self, tag):
        if tag == 'h1' and self._heading is not None:
            self.headings.append(''.join(self._heading).strip())
            self._heading = None


class MailTemplateTests(unittest.TestCase):
    def setUp(self):
        self.config = {'public_url': 'https://example.org/'}
        self.smtp = {'from_address': 'reports@example.org', 'provider': 'resend'}
        self.secret = 'mail-template-test-secret'

    def job(self, jid='one', status='succeeded', **changes):
        row = {
            'id': jid,
            'batch_id': 'submission-one',
            'status': status,
            'email': 'user@example.org',
            'filename': jid + '.fastq',
            'options': json.dumps({'sample': 'Sample ' + jid, 'assay': 'microc', 'reference': 'dm6'}),
            'created': 1789470000,
            'expires': 1790074800,
        }
        row.update(changes)
        return row

    def message(self, row, rows=None):
        return build_message(self.config, self.smtp, row, self.secret, rows=rows)

    def bodies(self, message):
        return (
            message.get_body(preferencelist=('plain',)).get_content(),
            message.get_body(preferencelist=('html',)).get_content(),
        )

    def task_link(self, row):
        return 'https://example.org/jobs/' + row['id'] + '#' + token_for(self.secret, row['id'])

    def assert_presentation(self, message, heading, action):
        plain, html = self.bodies(message)
        self.assertEqual(plain.splitlines()[0], heading)
        self.assertEqual(MailMarkup(html).headings, [heading])
        self.assertIn(action + ':', plain)
        self.assertIn(action, html)
        self.assertNotIn('\u2014', message['Subject'])
        for body in (plain, html):
            for redundant in ('LIBRARY 01', 'YOUR LIBRARY', 'ANALYSIS COMPLETE'):
                self.assertNotIn(redundant, body)

    def test_single_success_uses_singular_copy_for_legacy_and_batch_calls(self):
        row = self.job()
        for rows in (None, [row]):
            with self.subTest(batch=rows is not None):
                message = self.message(row, rows)
                self.assertEqual(message['Subject'], 'ONT-LibQC - report ready')
                self.assert_presentation(message, 'Report ready', 'Open report')
                plain, html = self.bodies(message)
                self.assertEqual(set(MailMarkup(html).links), {self.task_link(row)})
                self.assertNotIn('/#jobs=', plain)
                for body in (plain, html):
                    for plural in ('1 files', 'Reports ready', 'View all results', 'compare your libraries'):
                        self.assertNotIn(plural, body)

    def test_single_failure_uses_failure_copy_for_legacy_and_batch_calls(self):
        row = self.job(status='failed')
        for rows in (None, [row]):
            with self.subTest(batch=rows is not None):
                message = self.message(row, rows)
                self.assertEqual(message['Subject'], 'ONT-LibQC - analysis failed')
                self.assert_presentation(message, 'Analysis failed', 'View analysis')
                for body in self.bodies(message):
                    self.assertIn('Failed', body)
                    for misleading in ('Report ready', 'View all results', 'compare your libraries'):
                        self.assertNotIn(misleading, body)

    def test_single_cancelled_analysis_is_not_presented_as_failed_or_ready(self):
        row = self.job(status='cancelled')
        for rows in (None, [row]):
            with self.subTest(batch=rows is not None):
                message = self.message(row, rows)
                self.assertEqual(message['Subject'], 'ONT-LibQC - analysis cancelled')
                self.assert_presentation(message, 'Analysis cancelled', 'View analysis')
                for body in self.bodies(message):
                    self.assertIn(self.task_link(row), body)
                    self.assertNotIn('Report ready', body)
                    self.assertNotIn('Analysis failed', body)

    def test_single_expired_report_distinguishes_manual_deletion_from_retention(self):
        for expires in (0, 1790074800):
            row = self.job(status='expired', expires=expires)
            for rows in (None, [row]):
                with self.subTest(manual=expires == 0, batch=rows is not None):
                    message = self.message(row, rows)
                    self.assertEqual(message['Subject'], 'ONT-LibQC - report unavailable')
                    self.assert_presentation(message, 'Report unavailable', 'View analysis')
                    for body in self.bodies(message):
                        self.assertIn(self.task_link(row), body)
                        self.assertNotIn('Available until', body)
                        self.assertNotIn('1970', body)
                        self.assertNotIn('Report ready', body)
                        if expires:
                            self.assertIn('Expired on', body)
                        else:
                            self.assertIn('deleted', body.lower())
                            self.assertNotIn('Expired on', body)

    def test_completed_batch_has_plural_subject_and_all_private_links(self):
        rows = [self.job(jid) for jid in ('one', 'two', 'three')]
        message = self.message(rows[0], rows)
        self.assertEqual(message['Subject'], 'ONT-LibQC - 3 reports ready')
        self.assert_presentation(message, 'Reports ready', 'View all results')
        batch_link = 'https://example.org/#jobs=' + ','.join(
            row['id'] + ':' + token_for(self.secret, row['id']) for row in rows
        )
        plain, html = self.bodies(message)
        markup = MailMarkup(html)
        self.assertIn(batch_link, markup.links)
        for body in (plain, html):
            self.assertIn('View all results', body)
            self.assertIn(batch_link, body)
            for row in rows:
                self.assertIn(self.task_link(row), body)
                self.assertIn(row['filename'], body)
        for row in rows:
            self.assertIn(self.task_link(row), markup.links)

    def test_mixed_batch_reports_each_terminal_status_without_claiming_all_ready(self):
        rows = [self.job(str(i), status) for i, status in enumerate(('succeeded', 'failed', 'cancelled', 'expired'))]
        message = self.message(rows[0], rows)
        self.assertEqual(message['Subject'], 'ONT-LibQC - 4 analyses finished')
        self.assert_presentation(message, 'Submission finished', 'View all results')
        for body in self.bodies(message):
            self.assertNotIn('Reports ready', body)
            for label in ('Report ready', 'Failed', 'Cancelled', 'Expired', 'Expired on'):
                self.assertIn(label, body)
            for row in rows:
                self.assertIn(self.task_link(row), body)

    def test_untrusted_sample_filename_and_reference_are_escaped_in_html(self):
        sample = '<script>alert("sample")</script> & Sample'
        filename = '<img src=x onerror="bad()">.fastq'
        reference = 'reference <b>dm6</b> & "custom"'
        row = self.job(filename=filename, options=json.dumps({'sample': sample, 'assay': 'microc', 'reference': reference}))
        plain, html = self.bodies(self.message(row, [row]))
        for value in (sample, filename, reference):
            self.assertIn(value, plain)
            self.assertIn(escape(value), html)
            self.assertNotIn(value, html)
        markup = MailMarkup(html)
        self.assertNotIn('script', markup.tags)
        self.assertEqual(markup.images, ['cid:ontlibqc-logo'])
        self.assertIn(self.task_link(row), markup.links)

    def test_single_batch_keeps_batch_delivery_identity_after_singular_presentation(self):
        row = self.job()
        message = self.message(row, [row])
        self.assertEqual(message['Message-ID'], '<ont-libqc-batch-submission-one-complete@example.org>')
        self.assertEqual(message['Resend-Idempotency-Key'], 'ont-libqc/batch-submission-one/complete')
        self.assertEqual(message.get_boundary(), 'ontlibqc_alternative_batch-submission-one')
        related = next(part for part in message.walk() if part.get_content_type() == 'multipart/related')
        self.assertEqual(related.get_boundary(), 'ontlibqc_related_batch-submission-one')
        self.assertEqual(message.as_bytes(), self.message(row, [row]).as_bytes())

    def test_legacy_delivery_identity_still_uses_job_and_status(self):
        for status in ('succeeded', 'failed', 'cancelled', 'expired'):
            with self.subTest(status=status):
                row = self.job(status=status)
                message = self.message(row)
                self.assertEqual(message['Message-ID'], f'<ont-libqc-one-{status}@example.org>')
                self.assertEqual(message['Resend-Idempotency-Key'], f'ont-libqc/one/{status}')
                self.assertEqual(message.get_boundary(), 'ontlibqc_alternative_one')
                related = next(part for part in message.walk() if part.get_content_type() == 'multipart/related')
                self.assertEqual(related.get_boundary(), 'ontlibqc_related_one')
                self.assertEqual(message.as_bytes(), self.message(row).as_bytes())

    def test_plain_html_and_embedded_logo_remain_available_without_external_images(self):
        row = self.job()
        message = self.message(row, [row])
        self.assertEqual(message.get_content_type(), 'multipart/alternative')
        self.assertEqual(message['From'], 'reports@example.org')
        self.assertEqual(message['To'], 'user@example.org')
        plain, html = self.bodies(message)
        self.assertIn(self.task_link(row), plain)
        self.assertIn('private', plain.lower())
        self.assertIn('private', html.lower())
        logos = [part for part in message.walk() if part.get_content_type() == 'image/png']
        self.assertEqual(len(logos), 1)
        self.assertEqual(logos[0]['Content-ID'], '<ontlibqc-logo>')
        self.assertEqual(logos[0].get_content_disposition(), 'inline')
        self.assertTrue(logos[0].get_payload(decode=True).startswith(b'\x89PNG\r\n\x1a\n'))
        self.assertEqual(MailMarkup(html).images, ['cid:ontlibqc-logo'])

    def test_smtp_provider_does_not_receive_resend_specific_header(self):
        self.smtp['provider'] = 'smtp'
        row = self.job()
        self.assertNotIn('Resend-Idempotency-Key', self.message(row, [row]))


if __name__ == '__main__':
    unittest.main()
