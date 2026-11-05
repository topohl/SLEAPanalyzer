"""Draw figure stills: one video frame per row of stills.csv (validation/exp9_boris/scripts/27_figure_stills.R), the nose as
tracked (filled yellow dot) and as the analysis uses it (cyan ring), cropped to a square around both (and the cup, if given).
No text is drawn: labels belong to the figure.

usage:
  python -I tools/figure_stills.py <stills.csv> <out_dir> [--size 900] [--margin-cm 5] [--min-cm 20]
"""
import argparse, csv, os
import numpy as np
import cv2

TRACKED = (78, 202, 242)   # BGR of #F2CA4E
USED = (247, 195, 79)      # BGR of #4FC3F7


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stills"); ap.add_argument("out_dir"); ap.add_argument("--size", type=int, default=900)
    ap.add_argument("--margin-cm", type=float, default=5); ap.add_argument("--min-cm", type=float, default=20)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    for row in csv.DictReader(open(a.stills, encoding="utf-8")):
        cap = cv2.VideoCapture(row["video"]); cap.set(cv2.CAP_PROP_POS_FRAMES, int(row["frame"]) - 1)
        ok, img = cap.read(); cap.release()
        if not ok:
            raise SystemExit(f"cannot read frame {row['frame']} of {row['video']}")
        k = float(row["px_per_cm"]); h, w = img.shape[:2]
        pts = [(float(row["tracked_x"]), float(row["tracked_y"])), (float(row["used_x"]), float(row["used_y"]))]
        if row["cup_x"] not in ("", "NA"):
            pts.append((float(row["cup_x"]), float(row["cup_y"])))
        pts = np.array(pts)
        lo, hi = pts.min(0) - a.margin_cm * k, pts.max(0) + a.margin_cm * k
        side = min(max(max(hi - lo), a.min_cm * k), w, h); c = (lo + hi) / 2
        # a square that touches the frame edge is shifted inward, not cut
        x0 = int(min(max(0, c[0] - side / 2), w - side)); y0 = int(min(max(0, c[1] - side / 2), h - side))
        x1, y1 = int(x0 + side), int(y0 + side)
        tx, ty = pts[0]; ux, uy = pts[1]
        # marker sizes follow the crop, so they print at the same size on every still; the dot goes on top of the ring
        cv2.circle(img, (int(ux), int(uy)), int(0.024 * side), USED, max(2, int(0.007 * side)), cv2.LINE_AA)
        cv2.circle(img, (int(tx), int(ty)), int(0.012 * side), TRACKED, -1, cv2.LINE_AA)
        crop = img[y0:y1, x0:x1]
        crop = cv2.resize(crop, (a.size, int(a.size * crop.shape[0] / crop.shape[1])), interpolation=cv2.INTER_AREA)
        out = os.path.join(a.out_dir, row["still"] + ".png")
        cv2.imwrite(out, crop)
        print(out, crop.shape, "crop", x0, y0, x1, y1)


if __name__ == "__main__":
    main()
