"""Common labels and canvas for report figures."""
import html


def svg_start(width, height, title):
    return [f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' role='img' aria-label='{html.escape(title, quote=True)}' style='font-family:Arial,Helvetica,sans-serif'>",
            f'<title>{html.escape(title)}</title>', "<rect width='100%' height='100%' fill='white'/>",
            f"<text x='64' y='27' font-size='17' fill='#293f53' font-weight='600'>{html.escape(title)}</text>"]


def text(x, y, label, anchor='middle', size=12, color='#526579'):
    return f"<text x='{x:.3f}' y='{y:.3f}' text-anchor='{anchor}' font-size='{size}' fill='{color}'>{html.escape(str(label))}</text>"
