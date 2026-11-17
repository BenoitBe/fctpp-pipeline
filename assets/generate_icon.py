"""Generate the Olympus icon (200x240 card) and mark (48x48), same family as Deimos/Pathfinder.

    python assets/generate_icon.py          # writes olympus_icon.svg and olympus_mark.svg at the repo root

Concept: Olympus Mons on the Martian horizon. Its right flank is a melting curve
(soluble fraction vs temperature); a second, golden curve shifted to the right
shows a stabilised protein (ΔTm). Twelve dots under the volcano are the
temperature points, coloured from cold to hot; the spectrum sticks at the top are
the family signature, tinted along the same thermal gradient.
"""

from __future__ import annotations

import math
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

COLD, WARM, HOT = (79, 143, 214), (232, 192, 106), (224, 85, 47)


def thermal(t: float) -> str:
    """Cold blue -> warm gold -> hot red, t in [0, 1]."""
    a, b, u = (COLD, WARM, t / 0.5) if t < 0.5 else (WARM, HOT, (t - 0.5) / 0.5)
    r, g, bl = (round(a[i] + (b[i] - a[i]) * u) for i in range(3))
    return f"#{r:02x}{g:02x}{bl:02x}"


def sigmoid_y(x: float, tm: float, top: float, base: float, s: float) -> float:
    """Melting curve in SVG coordinates (y grows downwards): top at low T, base at high T."""
    return top + (base - top) / (1 + math.exp(-(x - tm) / s))


def volcano_path(cx: float, top: float, base: float, tm_off: float, s: float,
                 x_end: float, step: float = 1.0) -> str:
    """Symmetric shield whose flanks are melting curves; scarp at the foot."""
    right = []
    x = cx
    while x <= x_end:
        right.append((x, sigmoid_y(x, cx + tm_off, top, base, s)))
        x += step
    left = [(2 * cx - px, py) for px, py in reversed(right)]
    pts = left + right[1:]
    d = f"M {pts[0][0]:.1f} {base + 0.6:.1f} L " + " L ".join(f"{px:.1f} {py:.2f}" for px, py in pts)
    d += f" L {pts[-1][0]:.1f} {base + 0.6:.1f} Z"
    return d


def curve_path(x0: float, x1: float, tm: float, top: float, base: float, s: float,
               step: float = 1.0) -> str:
    pts = []
    x = x0
    while x <= x1 + 1e-9:
        pts.append((x, sigmoid_y(x, tm, top, base, s)))
        x += step
    return "M " + " L ".join(f"{px:.1f} {py:.2f}" for px, py in pts)


def spectrum(x0: float, x1: float, y_base: float, h_max: float, n: int, seed: int,
             width: float = 0.8, opacity: float = 0.85) -> str:
    rnd = random.Random(seed)
    out = []
    for i in range(n):
        x = x0 + (x1 - x0) * (i + rnd.uniform(-0.3, 0.3)) / (n - 1)
        h = h_max * (rnd.random() ** 2.2) + 1.2
        if rnd.random() < 0.12:
            h = h_max * rnd.uniform(0.75, 1.0)          # a few base peaks
        col = thermal((x - x0) / (x1 - x0))
        out.append(f'<line x1="{x:.1f}" y1="{y_base:.1f}" x2="{x:.1f}" y2="{y_base - h:.1f}" '
                   f'stroke="{col}" stroke-width="{width}" stroke-opacity="{opacity}"/>')
    return "\n      ".join(out)


def stars(n: int, seed: int, w: float, h: float, r_max: float = 0.8) -> str:
    rnd = random.Random(seed)
    return "".join(
        f'<circle cx="{rnd.uniform(4, w - 4):.1f}" cy="{rnd.uniform(4, h):.1f}" '
        f'r="{rnd.uniform(0.3, r_max):.2f}" opacity="{rnd.uniform(0.35, 0.9):.2f}"/>'
        for _ in range(n))


def icon() -> str:
    # a broad shield: wide summit plateau, gentle flanks (Olympus Mons is ~600 km wide)
    cx, top, base, s = 88.0, 120.0, 178.0, 8.0
    tm_ctrl, tm_shift = 50.0, 67.0            # offsets from the summit (control, stabilised)
    body = volcano_path(cx, top, base, tm_ctrl, s, x_end=210)
    ctrl = curve_path(cx, 204, cx + tm_ctrl, top, base, s)
    shifted = curve_path(cx + 8, 204, cx + tm_shift, top, base, s)
    y_half = (top + base) / 2
    dots = []
    for i in range(12):
        x = 34 + i * (132 / 11)
        dots.append(f'<circle cx="{x:.1f}" cy="189" r="1.55" fill="{thermal(i / 11)}"/>')
    pts_shift = []
    for i in range(12):
        x = 100 + i * 8.8
        pts_shift.append(f'<circle cx="{x:.1f}" cy="{sigmoid_y(x, cx + tm_shift, top, base, s):.2f}" '
                         f'r="1.05" fill="#f6dca0"/>')
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 200 240" width="200" height="240" role="img" aria-label="Olympus">
  <defs>
    <radialGradient id="sky" cx="50%" cy="42%" r="80%">
      <stop offset="0%" stop-color="#0e1524"/>
      <stop offset="100%" stop-color="#05070f"/>
    </radialGradient>
    <linearGradient id="ground" x1="0" y1="0" x2="0" y2="1">
      <stop offset="0%" stop-color="#b0603a"/>
      <stop offset="55%" stop-color="#8a4529"/>
      <stop offset="100%" stop-color="#4e2717"/>
    </linearGradient>
    <linearGradient id="shield" x1="0" y1="0" x2="1" y2="0.35">
      <stop offset="0%" stop-color="#6e3820"/>
      <stop offset="45%" stop-color="#b8693f"/>
      <stop offset="62%" stop-color="#9a522f"/>
      <stop offset="100%" stop-color="#4a2414"/>
    </linearGradient>
    <radialGradient id="heat" cx="50%" cy="50%" r="50%">
      <stop offset="0%" stop-color="#ff9a4a" stop-opacity="0.55"/>
      <stop offset="45%" stop-color="#e0552f" stop-opacity="0.18"/>
      <stop offset="100%" stop-color="#e0552f" stop-opacity="0"/>
    </radialGradient>
    <linearGradient id="haze" x1="0" y1="0" x2="0" y2="1">
      <stop offset="0%" stop-color="#b4786e" stop-opacity="0"/>
      <stop offset="100%" stop-color="#b4786e" stop-opacity="0.40"/>
    </linearGradient>
    <linearGradient id="gold" x1="0" y1="0" x2="1" y2="0">
      <stop offset="0%" stop-color="#f6dca0"/>
      <stop offset="100%" stop-color="#e8a24a"/>
    </linearGradient>
    <filter id="glow" x="-20%" y="-50%" width="140%" height="200%">
      <feGaussianBlur stdDeviation="1.1" result="b"/>
      <feMerge><feMergeNode in="b"/><feMergeNode in="SourceGraphic"/></feMerge>
    </filter>
    <clipPath id="card"><rect x="2" y="2" width="196" height="236" rx="16"/></clipPath>
  </defs>

  <g clip-path="url(#card)">
    <rect x="0" y="0" width="200" height="240" fill="url(#sky)"/>

    <!-- stars -->
    <g fill="#9fb6d4">{stars(26, 3, 200, 112)}</g>

    <!-- spectrum sticks (family signature), tinted cold -> hot -->
    <g>
      {spectrum(40, 160, 86, 32, 44, 11)}
    </g>

    <!-- volcanic heat over the caldera -->
    <ellipse cx="{cx}" cy="{top + 2}" rx="42" ry="16" fill="url(#heat)"/>

    <!-- plain -->
    <rect x="-5" y="{base:.1f}" width="210" height="70" fill="url(#ground)"/>
    <rect x="0" y="{base - 14:.0f}" width="200" height="15" fill="url(#haze)"/>

    <!-- Olympus Mons: its flanks are melting curves -->
    <path d="{body}" fill="url(#shield)"/>
    <ellipse cx="{cx - 2}" cy="{top + 0.7}" rx="11" ry="1.9" fill="#2a140b" opacity="0.85"/>
    <ellipse cx="{cx + 4}" cy="{top + 1.0}" rx="5" ry="1.0" fill="#160a05" opacity="0.9"/>

    <!-- control melting curve (the flank) and stabilised curve (ΔTm) -->
    <path d="{ctrl}" fill="none" stroke="#9fb6d4" stroke-width="1.1" stroke-opacity="0.9"/>
    <path d="{shifted}" fill="none" stroke="url(#gold)" stroke-width="1.6" filter="url(#glow)"/>
    <g>{"".join(pts_shift)}</g>
    <g stroke="#f6dca0" stroke-width="0.7" fill="#f6dca0">
      <line x1="{cx + tm_ctrl + 1.5:.1f}" y1="{y_half:.1f}" x2="{cx + tm_shift - 1.5:.1f}" y2="{y_half:.1f}"/>
      <path d="M {cx + tm_shift - 1.2:.1f} {y_half:.1f} l -2.4 -1.3 l 0 2.6 Z"/>
    </g>

    <!-- the temperature points -->
    <g>{"".join(dots)}</g>

    <text x="100" y="206" text-anchor="middle" fill="#8b949e" font-family="Georgia, 'Times New Roman', serif"
          font-size="4.6" letter-spacing="1.6">· THERMAL · PROTEOME · PROFILING ·</text>
    <text x="100" y="230" text-anchor="middle" fill="#c49c6e" font-family="Georgia, 'Times New Roman', serif"
          font-size="19" letter-spacing="3.4">OLYMPUS</text>
  </g>
  <rect x="2" y="2" width="196" height="236" rx="16" fill="none" stroke="#3c6082" stroke-opacity=".7"/>
</svg>
'''


def mark() -> str:
    cx, top, base, s = 21.0, 24.5, 37.0, 2.4
    body = volcano_path(cx, top, base, 12.0, s, x_end=50, step=0.5)
    shifted = curve_path(cx + 3, 49, cx + 17.0, top, base, s, step=0.5)
    return f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 48 48" width="48" height="48" role="img" aria-label="Olympus">
  <defs>
    <radialGradient id="mSky" cx="50%" cy="40%" r="75%">
      <stop offset="0%" stop-color="#121a2b"/>
      <stop offset="100%" stop-color="#070a12"/>
    </radialGradient>
    <linearGradient id="mGround" x1="0" y1="0" x2="0" y2="1">
      <stop offset="0%" stop-color="#b0603a"/>
      <stop offset="100%" stop-color="#5e2f1c"/>
    </linearGradient>
    <linearGradient id="mShield" x1="0" y1="0" x2="1" y2="0.35">
      <stop offset="0%" stop-color="#6e3820"/>
      <stop offset="50%" stop-color="#b8693f"/>
      <stop offset="100%" stop-color="#4a2414"/>
    </linearGradient>
    <radialGradient id="mHeat" cx="50%" cy="50%" r="50%">
      <stop offset="0%" stop-color="#ff9a4a" stop-opacity="0.6"/>
      <stop offset="100%" stop-color="#e0552f" stop-opacity="0"/>
    </radialGradient>
    <clipPath id="mCard"><rect x="1" y="1" width="46" height="46" rx="10"/></clipPath>
  </defs>
  <g clip-path="url(#mCard)">
    <rect width="48" height="48" fill="url(#mSky)"/>
    <g fill="#9fb6d4">{stars(6, 5, 48, 18, 0.6)}</g>
    <ellipse cx="{cx}" cy="{top + 1}" rx="10" ry="5" fill="url(#mHeat)"/>
    <rect x="-2" y="{base}" width="52" height="14" fill="url(#mGround)"/>
    <path d="{body}" fill="url(#mShield)"/>
    <ellipse cx="{cx}" cy="{top + 0.3}" rx="3" ry="0.7" fill="#2a140b"/>
    <path d="{shifted}" fill="none" stroke="#f0c070" stroke-width="1.3"/>
    <g>{"".join(f'<circle cx="{8 + i * 32 / 11:.1f}" cy="42.5" r="0.75" fill="{thermal(i / 11)}"/>' for i in range(12))}</g>
  </g>
  <rect x="1" y="1" width="46" height="46" rx="10" fill="none" stroke="#3c6082" stroke-opacity=".7"/>
</svg>
'''


if __name__ == "__main__":
    (ROOT / "olympus_icon.svg").write_text(icon(), encoding="utf-8")
    (ROOT / "olympus_mark.svg").write_text(mark(), encoding="utf-8")
    print("olympus_icon.svg, olympus_mark.svg written in", ROOT)
