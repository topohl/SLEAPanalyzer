"""Draw figure stills: one video frame per row of stills.csv (validation/exp9_boris/scripts/27_figure_stills.R), the nose as
tracked (filled dot) and as the analysis uses it (ring), the two ears the head is measured from (small open squares), and
where a cup is given its point (cross) and the contact distance (thin circle), cropped to a square around them. No text is
drawn: labels belong to the figure.

usage:
  python -I tools/figure_stills.py <stills.csv> <out_dir> [--size 900] [--margin-cm 5] [--min-cm 20]
                                   [--tracked "#FF8C1A"] [--analysed "#29B6F6"] [--reference "#FFFFFF"] [--contact-cm 6]
                                   [--inset-cm 8] [--inset-above-cm 30]
A still wider than inset-above-cm gets an inset of the head (inset-cm square around the analysed nose) at the top right.
Colours are hex RGB; pass the figure's palette so the stills and the figure keys agree.
"""
import argparse, csv, os
import numpy as np
import cv2


def bgr(hex_rgb):
    h = hex_rgb.lstrip("#")
    return (int(h[4:6], 16), int(h[2:4], 16), int(h[0:2], 16))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stills"); ap.add_argument("out_dir"); ap.add_argument("--size", type=int, default=900)
    ap.add_argument("--margin-cm", type=float, default=5); ap.add_argument("--min-cm", type=float, default=20)
    ap.add_argument("--tracked", default="#FF8C1A"); ap.add_argument("--analysed", default="#29B6F6")
    ap.add_argument("--reference", default="#FFFFFF"); ap.add_argument("--contact-cm", type=float, default=6)
    ap.add_argument("--inset-cm", type=float, default=8); ap.add_argument("--inset-above-cm", type=float, default=30)
    a = ap.parse_args()
    tracked_c, analysed_c, ref_c = bgr(a.tracked), bgr(a.analysed), bgr(a.reference)
    os.makedirs(a.out_dir, exist_ok=True)
    for row in csv.DictReader(open(a.stills, encoding="utf-8")):
        cap = cv2.VideoCapture(row["video"]); cap.set(cv2.CAP_PROP_POS_FRAMES, int(row["frame"]) - 1)
        ok, img = cap.read(); cap.release()
        if not ok:
            raise SystemExit(f"cannot read frame {row['frame']} of {row['video']}")
        k = float(row["px_per_cm"]); h, w = img.shape[:2]
        pts = [(float(row["tracked_x"]), float(row["tracked_y"])), (float(row["used_x"]), float(row["used_y"]))]
        cup = None
        if row["cup_x"] not in ("", "NA"):
            cup = (float(row["cup_x"]), float(row["cup_y"])); pts.append(cup)
        pts = np.array(pts)
        lo, hi = pts.min(0) - a.margin_cm * k, pts.max(0) + a.margin_cm * k
        side = min(max(max(hi - lo), a.min_cm * k), w, h); c = (lo + hi) / 2
        # a square that touches the frame edge is shifted inward, not cut
        x0 = int(min(max(0, c[0] - side / 2), w - side)); y0 = int(min(max(0, c[1] - side / 2), h - side))
        x1, y1 = int(x0 + side), int(y0 + side)
        clean = img.copy()

        def draw(im, side):   # marker sizes follow the view, so they print at the same size on every still; the dot goes on top
            lw = max(2, int(0.007 * side))
            if cup is not None:
                cv2.circle(im, (int(cup[0]), int(cup[1])), int(a.contact_cm * k), ref_c, max(1, lw // 2), cv2.LINE_AA)
                cv2.drawMarker(im, (int(cup[0]), int(cup[1])), ref_c, cv2.MARKER_CROSS, int(0.04 * side), max(1, lw // 2), cv2.LINE_AA)
            for e in ("lear", "rear"):
                if row.get(e + "_x") not in (None, "", "NA"):
                    ex, ey = float(row[e + "_x"]), float(row[e + "_y"]); q = max(3, int(0.008 * side))
                    cv2.rectangle(im, (int(ex) - q, int(ey) - q), (int(ex) + q, int(ey) + q), ref_c, max(1, lw // 2), cv2.LINE_AA)
            (tx, ty), (ux, uy) = pts[0], pts[1]
            cv2.circle(im, (int(ux), int(uy)), int(0.024 * side), analysed_c, lw, cv2.LINE_AA)
            cv2.circle(im, (int(tx), int(ty)), int(0.012 * side), tracked_c, -1, cv2.LINE_AA)
            return im

        crop = draw(img, side)[y0:y1, x0:x1]
        crop = cv2.resize(crop, (a.size, int(a.size * crop.shape[0] / crop.shape[1])), interpolation=cv2.INTER_AREA)
        if side > a.inset_above_cm * k:
            # a wide view leaves the head small: an inset of the head (inset_cm square around the analysed nose), top right
            ux, uy = pts[1]; r = a.inset_cm * k / 2
            ix0 = int(min(max(0, ux - r), w - 2 * r)); iy0 = int(min(max(0, uy - r), h - 2 * r))
            n = int(a.size * 0.42); sub = draw(clean, 2 * r * 2.2)[iy0:iy0 + int(2 * r), ix0:ix0 + int(2 * r)]
            inset = cv2.resize(sub, (n, n), interpolation=cv2.INTER_CUBIC)
            b = max(2, a.size // 300); crop[b:b + n, a.size - n - b:a.size - b] = inset
            cv2.rectangle(crop, (a.size - n - b, b), (a.size - b - 1, b + n - 1), ref_c, b, cv2.LINE_AA)
        out = os.path.join(a.out_dir, row["still"] + ".png")
        cv2.imwrite(out, crop)
        print(out, crop.shape, "crop", x0, y0, x1, y1)


if __name__ == "__main__":
    main()
