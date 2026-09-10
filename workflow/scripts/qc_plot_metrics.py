"""Depth-normalized concentration and ONT-derived electropherogram-style plots."""
import collections
import html
import math
import textwrap


def normalized_collision_stats(records, sizes, bin_size):
    """Pair U-statistic divided by the uniform-base null probability.

    E[O / choose(N, 2)] = sum_j p_j**2 under independent sampling.
    Actual terminal-bin widths define the null, rather than equal-size bins.
    """
    if bin_size <= 0 or not sizes or any(size <= 0 for size in sizes.values()):
        raise ValueError("Collision metrics require positive bin and chromosome sizes")
    genome_size = sum(sizes.values())
    total_bins = sum(math.ceil(size / bin_size) for size in sizes.values())
    width_squared_sum = 0
    for size in sizes.values():
        full, remainder = divmod(size, bin_size)
        width_squared_sum += full * bin_size**2 + remainder**2
    null_probability = width_squared_sum / genome_size**2
    counts = collections.Counter()
    for rec in records:
        chrom = rec["chrom"]
        midpoint = (rec["start"] + rec["end"]) // 2
        if chrom not in sizes or not 0 <= midpoint < sizes[chrom]:
            raise ValueError(f"Collision midpoint outside reference: {chrom}:{midpoint}")
        counts[chrom, midpoint // bin_size] += 1
    n = sum(counts.values())
    possible_pairs = n * (n - 1) // 2
    observed = sum(v * (v - 1) // 2 for v in counts.values())
    expected = possible_pairs * null_probability
    probability = observed / possible_pairs if possible_pairs else None
    oe = observed / expected if expected else None
    return dict(bin_size=bin_size, fragments=n, genome_bins=total_bins,
                observed_pairs=observed, expected_pairs=expected,
                null_expected_pairs=expected, null_pair_probability=null_probability,
                collision_probability=probability, pair_oe=oe,
                normalized_excess=oe-1 if oe is not None else None)


def collision_metric_fields(stats):
    return dict(collision_bin_size=stats["bin_size"],
                collision_oe=stats["pair_oe"], collision_normalized_excess=stats["normalized_excess"],
                collision_probability=stats["collision_probability"],
                collision_null_pair_probability=stats["null_pair_probability"],
                collision_observed_pairs=stats["observed_pairs"],
                collision_null_expected_pairs=stats["null_expected_pairs"],
                collision_fragments=stats["fragments"],
                collision_normalization="pair_probability / uniform_base_null; terminal-bin-width corrected")


def size_trace_data(values, maximum_bp=2000, bandwidth=.045, points=650, weighting="count"):
    """Count- or length-weighted smoothing on log1p(bp/50); peak = 100.

    Fixed bandwidth keeps smoothing independent of sample depth. Values above
    the visible range are counted but neither clamped nor added to the last bin.
    Weight is 1 per observation or L per fragment of length L. Counts remain
    observation counts in both modes. This is not calibrated DNA fluorescence.
    """
    if maximum_bp <= 0 or bandwidth <= 0 or points < 3:
        raise ValueError("Size trace requires positive display range/bandwidth and >=3 grid points")
    if weighting not in {"count", "length"}:
        raise ValueError("Size trace weighting must be count or length")
    values = list(values)
    if any(not math.isfinite(v) or v < 0 for v in values):
        raise ValueError("Size trace lengths must be finite and non-negative")
    upper = math.log1p(maximum_bp / 50)
    step = upper / (points-1)
    counts = [0] * points
    shown = 0
    shown_weight = 0
    for value in values:
        if value <= maximum_bp:
            weight = value if weighting == "length" else 1
            counts[min(points-1, round(math.log1p(value / 50) / step))] += weight
            shown += 1
            shown_weight += weight
    radius = math.ceil(4*bandwidth/step)
    kernel = [math.exp(-.5*(offset*step/bandwidth)**2) for offset in range(-radius,radius+1)]
    signal = []
    for i in range(points):
        # Reflect at zero; beyond the upper display edge is omitted, not reflected.
        signal.append(sum(weight*counts[abs(i+offset)] for offset,weight in zip(range(-radius,radius+1),kernel)
                          if abs(i+offset) < points))
    peak = max(signal) or 1
    return dict(lengths=[50*math.expm1(i*step) for i in range(points)],
                signal=[100*s/peak for s in signal], n=len(values),
                shown=shown, above_range=sum(v > maximum_bp for v in values),
                zero_count=sum(v == 0 for v in values), bandwidth=bandwidth,
                maximum_bp=maximum_bp, weighting=weighting,
                weight_sum=sum(values) if weighting == "length" else len(values),
                shown_weight=shown_weight,
                above_range_weight=sum(v if weighting == "length" else 1 for v in values if v > maximum_bp))


def svg_size_trace(values, title, maximum_bp=2000, bandwidth=.045, footnote=None, weighting="count"):
    values = list(values)
    data = size_trace_data(values, maximum_bp, bandwidth, weighting=weighting)
    weighting_label = "length-weighted" if weighting == "length" else "count-based"
    color = "#008577" if weighting == "length" else "#176cad"
    width,height,left,top,plot_w,plot_h=940,570,78,82,690,365
    baseline = top+plot_h
    if not values:
        return f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' role='img' data-weighting='{weighting}' aria-label='{html.escape(title,quote=True)}'><rect width='100%' height='100%' fill='white'/><text x='40' y='65'>{html.escape(title)}: no observations</text></svg>"
    elements=[f"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 {width} {height}' role='img' data-weighting='{weighting}' aria-label='{html.escape(title,quote=True)}' style='font-family:Arial,Helvetica,sans-serif'>",
              "<rect width='100%' height='100%' fill='white'/>",
              f"<text x='{left}' y='28' font-size='19' font-weight='600'>{html.escape(title)}</text>",
              f"<text x='{left}' y='49' font-size='13' fill='#576574'>ONT-derived · {weighting_label} · n={data['n']:,} · zero bp={data['zero_count']:,} · &gt;{maximum_bp:,.0f} bp={data['above_range']:,} (outside display)</text>"]
    for value in [0,25,50,75,100]:
        y=baseline-value/100*plot_h
        elements += [f"<line x1='{left}' x2='{left+plot_w}' y1='{y}' y2='{y}' stroke='#e1e7eb'/>",
                     f"<text x='{left-12}' y='{y+4}' text-anchor='end' font-size='13' fill='#485460'>{value}</text>"]
    upper=math.log1p(maximum_bp/50)
    ticks=[0,50,100,150,200,300,500,700,1000,1500,2000,3000,5000,10000]
    for value in ticks:
        if value > maximum_bp: continue
        x=left+math.log1p(value/50)/upper*plot_w
        elements += [f"<line x1='{x:.2f}' x2='{x:.2f}' y1='{top}' y2='{baseline}' stroke='#edf0f2'/>",
                     f"<text x='{x:.2f}' y='{baseline+21}' text-anchor='middle' font-size='13'>{value:,}</text>"]
    points=' '.join(f"{left+i*plot_w/(len(data['signal'])-1):.2f},{baseline-s/100*plot_h:.2f}" for i,s in enumerate(data['signal']))
    elements += [f"<polygon points='{left},{baseline} {points} {left+plot_w},{baseline}' fill='{color}' opacity='.12'/>",
                 f"<polyline points='{points}' fill='none' stroke='{color}' stroke-width='2.0' stroke-linejoin='round'/>",
                 f"<line x1='{left}' x2='{left+plot_w}' y1='{baseline}' y2='{baseline}' stroke='#586069'/>",
                 f"<text x='23' y='{top+plot_h/2}' transform='rotate(-90 23 {top+plot_h/2})' text-anchor='middle' font-size='14'>Relative signal (a.u.; peak = 100)</text>",
                 f"<text x='{left+plot_w/2}' y='{baseline+45}' text-anchor='middle' font-size='14'>Fragment size (bp; nonlinear size axis)</text>"]
    # The right-hand lane uses the same signal and nonlinear size mapping.
    # Larger fragments remain near the top, as in an electrophoresis lane.
    lane_x, lane_w = left+plot_w+50, 32
    elements += ["<g id='virtual-size-lane'>",
                 f"<title>Virtual size lane: larger fragments at top; darkness follows the {weighting_label} ONT trace.</title>",
                 f"<text x='{lane_x+lane_w/2}' y='{top-15}' text-anchor='middle' font-size='13' fill='#485460'>Virtual lane</text>"]
    count = len(data['signal'])
    step = plot_h/(count-1)
    for i, signal in enumerate(data['signal']):
        shade = round(248-205*signal/100)
        center = baseline-i*step
        y0, y1 = max(top, center-step/2), min(baseline, center+step/2)
        elements.append(f"<rect x='{lane_x}' y='{y0:.3f}' width='{lane_w}' height='{y1-y0:.3f}' fill='rgb({shade},{shade},{shade})'/>")
    elements.append(f"<rect x='{lane_x}' y='{top}' width='{lane_w}' height='{plot_h}' fill='none' stroke='#aab4bc' stroke-width='.8'/>")
    lane_ticks = sorted({v for v in [0,100,300,1000,maximum_bp] if v <= maximum_bp})
    for value in lane_ticks:
        y = baseline-math.log1p(value/50)/upper*plot_h
        elements += [f"<line x1='{lane_x+lane_w}' x2='{lane_x+lane_w+5}' y1='{y:.2f}' y2='{y:.2f}' stroke='#586069'/>",
                     f"<text x='{lane_x+lane_w+9}' y='{y+4:.2f}' font-size='10' fill='#485460'>{value:,.0f}</text>"]
    elements.append(f"<text x='{lane_x+lane_w/2}' y='{baseline+23}' text-anchor='middle' font-size='13' fill='#485460'>bp</text></g>")
    note=footnote or ("Length-weighted full-library model; relative mass proxy. Not measured TapeStation fluorescence."
                     if weighting == "length" else
                     "Count-based ONT signal; fixed smoothing; virtual lane. Not measured TapeStation fluorescence.")
    if len(values)<20: note="Low observation count: trace shape is uncertain. " + note
    elements.append(f"<desc>{html.escape(note)}</desc>")
    note_lines = textwrap.wrap(note, width=125)
    for i, line in enumerate(note_lines):
        elements.append(f"<text x='{left}' y='{height-10-13*(len(note_lines)-1-i)}' font-size='10' fill='#576574'>{html.escape(line)}</text>")
    elements.append("</svg>")
    return ''.join(elements)
