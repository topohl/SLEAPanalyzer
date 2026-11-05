"""Draw figure stills: one video frame per row of stills.csv (validation/exp9_boris/scripts/27_figure_stills.R), the nose as
tracked (filled dot) and as the analysis uses it (ring), and where a cup is given its point (cross) and the contact distance
(thin circle), cropped to a square around them. No text is drawn: labels belong to the figure.

usage:
  python -I tools/figure_stills.py <stills.csv> <out_dir> [--size 900] [--margin-cm 5] [--min-cm 20]
                                   [--tracked "#FF8C1A"] [--analysed "#29B6F6"] [--reference "#FFFFFF"] [--contact-cm 6]
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
        lw = max(2, int(0.007 * side))
        if cup is not None:
            cv2.circle(img, (int(cup[0]), int(cup[1])), int(a.contact_cm * k), ref_c, max(1, lw // 2), cv2.LINE_AA)
            cv2.drawMarker(img, (int(cup[0]), int(cup[1])), ref_c, cv2.MARKER_CROSS, int(0.04 * side), max(1, lw // 2), cv2.LINE_AA)
        tx, ty = pts[0]; ux, uy = pts[1]
        # marker sizes follow the crop, so they print at the same size on every still; the dot goes on top of the ring
        cv2.circle(img, (int(ux), int(uy)), int(0.024 * side), analysed_c, lw, cv2.LINE_AA)
        cv2.circle(img, (int(tx), int(ty)), int(0.012 * side), tracked_c, -1, cv2.LINE_AA)
        crop = img[y0:y1, x0:x1]
        crop = cv2.resize(crop, (a.size, int(a.size * crop.shape[0] / crop.shape[1])), interpolation=cv2.INTER_AREA)
        out = os.path.join(a.out_dir, row["still"] + ".png")
        cv2.imwrite(out, crop)
        print(out, crop.shape, "crop", x0, y0, x1, y1)


if __name__ == "__main__":
    main()
