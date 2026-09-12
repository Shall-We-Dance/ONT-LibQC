"""Report asset types shared by the dispatcher, transfer agent and HTTP service."""

# Configured custom BED plots include PDF alongside the SVG displayed in HTML.
REPORT_FILE_SUFFIXES = frozenset({
    '.html', '.svg', '.css', '.js', '.json', '.tsv', '.csv', '.txt', '.md',
    '.bed', '.gz', '.png', '.pdf',
})
