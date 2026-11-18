"""Stream the portable report with a snapshot of its saved display names."""
import zipfile
from starlette.responses import JSONResponse, StreamingResponse

from .sublibrary_names import MAP_FILE, annotate_page, name_map_tsv


class NamedArchiveResponse(StreamingResponse):
    def __init__(self, download, annotations, slots):
        self.slots = slots
        self.chunks = archive_chunks(download.path, annotations)
        super().__init__(self.chunks, media_type='application/zip',
                         headers={'Content-Disposition': download.headers['Content-Disposition']})

    async def __call__(self, scope, receive, send):
        if self.slots.locked():
            return await JSONResponse({'detail': 'Downloads are busy. Try again shortly.'}, 503,
                                      headers={'Retry-After': '3'})(scope, receive, send)
        async with self.slots:
            try:
                await super().__call__(scope, receive, send)
            finally:
                self.chunks.close()


class ChunkBuffer:
    """A non-seekable ZIP output. Drain it after each bounded input chunk."""
    def __init__(self):
        self.parts = []

    def write(self, data):
        if data:
            self.parts.append(data)
        return len(data)

    def flush(self):
        pass

    def drain(self):
        data = b''.join(self.parts)
        self.parts.clear()
        return data


def archive_chunks(path, annotations):
    sink = ChunkBuffer()
    pages = {f'report/samples/{sample}.html': payload for sample, payload in annotations.items()}
    map_path = 'report/' + MAP_FILE
    with zipfile.ZipFile(path) as source, zipfile.ZipFile(sink, 'w', zipfile.ZIP_DEFLATED, compresslevel=1) as target:
        for info in source.infolist():
            if info.filename == map_path:
                continue
            if info.filename in pages:
                catalog, names = pages[info.filename]
                body = source.read(info).decode('utf-8')
                target.writestr(info.filename, annotate_page(body, catalog, names))
            elif info.filename == 'report/downloads.html':
                body = source.read(info).decode('utf-8')
                link = '<p><a href="' + MAP_FILE + '">Sublibrary names (TSV)</a></p>'
                target.writestr(info.filename, body.replace('</main>', link + '</main>', 1))
            else:
                with source.open(info) as reader, target.open(info.filename, 'w') as writer:
                    while chunk := reader.read(128 * 1024):
                        writer.write(chunk)
                        if sink.parts:
                            yield sink.drain()
            if sink.parts:
                yield sink.drain()
        target.writestr(map_path, name_map_tsv(annotations))
    yield sink.drain()
