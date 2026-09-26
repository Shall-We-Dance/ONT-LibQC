"""FASTQ header compatibility and validation before analysis."""
import gzip
from pathlib import Path
import tempfile
import unittest

from ont_portal.run_job import InputError, normalize_fastq

READ_ID = b'00000000-0000-4000-8000-000000000001'


def record(header, seq=b'ACGTN', quality=b'IIIII', newline=b'\n'):
    return newline.join([b'@' + header, seq, b'+', quality, b''])


class FastqNormalizationTests(unittest.TestCase):
    def normalize(self, data, compressed=False, **limits):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, dest = root/'input.upload', root/'normalized.fastq.gz'
            source.write_bytes(gzip.compress(data) if compressed else data)
            count = normalize_fastq(source, dest, compressed,
                                    limits.get('max_bytes', 100000), limits.get('max_reads', 100))
            return count, gzip.decompress(dest.read_bytes())

    def test_plain_and_gzip_headers_with_space_or_tab_metadata(self):
        # SAM-style tags include punctuation that is valid metadata, but would
        # be rejected if accidentally parsed as part of the identifier.
        for compressed in [False, True]:
            for suffix in [b'', b' runid=example ch=12',
                           b'\tqs:f:16.8\tch:i:12\tRG:Z:run_dna_r10.4.1_e8.2_400bps_sup@v5.2.0_barcode16',
                           b' \tch:i:12\tSM:Z:barcode16']:
                for newline in [b'\n', b'\r\n']:
                    with self.subTest(compressed=compressed, suffix=suffix, newline=newline):
                        count, data = self.normalize(record(READ_ID + suffix, newline=newline), compressed)
                        self.assertEqual(count, 1)
                        self.assertEqual(data, record(READ_ID))

    def test_duplicate_ids_are_detected_across_metadata_formats(self):
        data = record(READ_ID + b' ch=12') + record(READ_ID + b'\tch:i:13')
        with self.assertRaisesRegex(InputError, 'Duplicate read identifier at read 2'):
            self.normalize(data)

    def test_invalid_identifiers_are_not_hidden_by_metadata(self):
        for identifier in [b'', b' ', b'\t', b'bad/name', b'bad;name', b'bad\x00name', b'x'*256]:
            with self.subTest(identifier=identifier):
                with self.assertRaisesRegex(InputError, 'Invalid FASTQ identifier at read 1'):
                    self.normalize(record(identifier + b'\tch:i:12'))
        count, data = self.normalize(record(b'x'*255 + b'\tch:i:12'))
        self.assertEqual((count, data), (1, record(b'x'*255)))

    def test_tabs_do_not_bypass_sequence_and_quality_validation(self):
        for seq, quality in [(b'ACGT!', b'IIIII'), (b'ACGTN', b'IIII'), (b'ACGTN', b'III I')]:
            with self.subTest(seq=seq, quality=quality):
                with self.assertRaisesRegex(InputError, 'Invalid sequence or quality values at read 1'):
                    self.normalize(record(READ_ID + b'\tch:i:12', seq, quality))

    def test_record_and_size_limits_still_apply(self):
        data = record(READ_ID + b'\tch:i:12')
        with self.assertRaisesRegex(InputError, 'exceeds 1 reads'):
            self.normalize(data + record(b'second\tch:i:13'), max_reads=1)
        with self.assertRaisesRegex(InputError, 'expanded FASTQ exceeds'):
            self.normalize(data, max_bytes=len(data)-1)
        with self.assertRaisesRegex(InputError, 'Invalid or excessively long FASTQ record'):
            self.normalize(data[:-1])
        with self.assertRaisesRegex(InputError, 'Invalid or excessively long FASTQ record'):
            self.normalize(record(READ_ID + b'\tcomment:Z:' + b'x'*4096))


if __name__ == '__main__':
    unittest.main()
