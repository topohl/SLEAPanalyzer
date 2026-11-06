"""objects: do the tracked points of fixed objects sit on the objects in the video, the same way in every group?

The points of fixed objects (manifest 'objects', e.g. objL objR, or socl socr for cups) place the zones an analysis
draws around the objects, at their median tracked position. Per recording the median of 15 frames spread over the
clip (the animal moves and drops out) is compared with the floor inside the arena outline, shrunk by 4% of the
arena's size so walls and corners do not count. A pixel differs from the floor when its Lab colour is far from the
floor's median colour, or when its lightness varies within 0.4 cm (edges and texture: grey objects on a grey floor,
wire cups). Each limit is the floor's median + 6 robust SDs, at least 8 (colour) and 3 (texture) Lab units. The
differing pixels, joined (closing 0.5 cm), with holes filled and at least 2 cm2, are the regions; a hole may be closed
by the outline (an object the outline cuts, e.g. a cup under a lid reaching past it) if under 10% of the floor. A
region includes what the camera sees of the object (its top, sides and shadow) and a rim of about 0.2 cm from the
texture window. Each point's median position is matched to the region it lies on, or the nearest one. An analysis
that holds the objects at their median over the frames in which the points are apart (e.g. static_objects) places a
point that jumps elsewhere than this median; object_jumps marks those recordings. An object the camera sees only in
part, the rest under something the floor's colour (a cup under a lid), gives a region of its visible part only, which
moves with the light while the point stays: leave such objects out of 'objects'.

Columns: group and object_set from the manifest, n_regions, regions_share (the share of the floor they cover), and per
object landmark <o>, distances in cm when the recording has a scale, else in px (the column name says):
  <o>_x, <o>_y             median tracked position, video pixels
  <o>_on_region            the median position lies on a region
  <o>_region               that region's number (1 = largest)
  <o>_to_region_cm         distance to the region's edge, 0 when on it
  <o>_dx_cm, <o>_dy_cm     offset of the median position from the region's centroid (image axes: x right, y down)
  <o>_region_area_cm2      the region's area
  <o>_near_median          share of the frames the point is tracked in with the point within 2 cm of its median
  <o>_group_shift_cm       (after all recordings) distance of the group's median offset from the other groups'
and per pair: <a>_<b>_sep_cm (median positions) and <a>_<b>_region_sep_cm (their regions' centroids).

Flags:
  off_object        a point is not on a region and more than 1 cm from the nearest one, or no region was found
  object_jumps      a point is within 2 cm of its median in under 80% of the frames: it sits on something else part
                    of the time
  objects_merged    two points match the same region: one sits on the other's object (or the two objects touch)
  missing_landmark  an object landmark is not in the tracking, or never tracked
  no_outline        the arena outline landmarks are not tracked: there is no floor to compare with, so the points
                    are not checked against the video
  uneven_floor      the regions cover more than a quarter of the floor (lines on it, bedding, uneven light): they are
                    not objects, so the points are not checked against them
  group_offset      the median offset of a landmark in this recording's group (manifest 'group') differs from the
                    median of the other groups' medians by more than --object-offset-cm: the points sit on the
                    objects differently in this group, or the objects look different in it (other objects, light,
                    reflections joining a region); the stills tell which. Groups are compared within the manifest's
                    optional 'object_set' (recordings with the same objects, e.g. one test and trial; without it every
                    recording with the landmark). Needs at least 2 groups; with exactly 2 both are flagged. Recordings
                    flagged off_object, objects_merged or missing_landmark do not count toward the medians, and only
                    recordings with that landmark are flagged.

Without a scale every cm limit is converted with --object-px-per-cm. Stills (--stills): the background with the
regions (red), tracked points (dots) and their median (cross), and a zoom on each object, for flagged recordings.
"""

from __future__ import annotations

import os
import re

import cv2
import numpy as np

import core

N_BACKGROUND = 15       # frames in the median background
WIDTH = 720             # working width; results are mapped back to video pixels
MARGIN = 0.04           # outline shrink, share of the arena's size (square root of its area)
ROBUST_SD = 6           # limits: floor median + this many robust SDs
MIN_COLOUR = 8.0        # Lab units
MIN_TEXTURE = 3.0       # SD of Lab lightness
TEXTURE_CM = 0.4
OPEN_CM, CLOSE_CM = 0.15, 0.5
MIN_AREA_CM2 = 2.0
HOLE_SHARE = 0.1        # enclosed floor patches below this share of the floor are filled
MAX_SHARE = 0.25        # regions covering more of the floor than this are not objects
NEAR_CM, MIN_NEAR = 2.0, 0.8
OFF_CM = 1.0
OFFSET_CM = 1.0         # default of --object-offset-cm
PX_PER_CM = 16.0        # default of --object-px-per-cm
COLOURS = [(0, 255, 255), (255, 255, 0), (255, 0, 255), (0, 255, 0), (0, 165, 255), (255, 128, 128)]
NOT_COUNTED = {"off_object", "objects_merged", "missing_landmark", "no_outline", "uneven_floor", "error", "unreadable"}


def run(rec, ctx):
    objs = core.names(rec, "objects")
    row = {"group": core.field(rec, "group", ""), "object_set": core.field(rec, "object_set", "")}
    if not objs:
        row.update(n_objects=0, flags=[])
        return row
    opts, track, v = ctx["opts"], ctx["track"], ctx["video"]
    scale = ctx["scale"]
    unit = "cm" if scale else "px"
    cm = scale or 1 / getattr(opts, "object_px_per_cm", PX_PER_CM)   # cm per video pixel, for the limits
    to_unit = (lambda px: px * scale) if scale else (lambda px: px)
    flags = []

    pts = {}
    for o in objs:
        p = core.median_point(track, o)
        if p is None:
            flags.append("missing_landmark")
        else:
            pts[o] = p
    row["cm_per_px"] = core.rnd(scale, 5)
    bg = None
    if ctx["poly"] is None:
        flags.append("no_outline")
    else:
        bg = v.background(N_BACKGROUND, width=min(WIDTH, v.width), gray=False)
        if bg is None:
            raise RuntimeError("no frames decoded")
        fx, fy = bg.shape[1] / v.width, bg.shape[0] / v.height
        regions, mask = find_regions(bg, ctx["poly"], fx, fy, cm)
        share = sum(r["area"] for r in regions) / max(mask.sum() / (fx * fy), 1)
        row.update(n_regions=len(regions), regions_share=core.rnd(share))
        if share > MAX_SHARE:
            flags.append("uneven_floor")
    judge = bg is not None and "uneven_floor" not in flags

    match = {}
    for o in objs:
        if o not in pts:
            continue
        p = pts[o]
        x, y = track[o]
        ok = np.isfinite(x) & np.isfinite(y)
        near = float(np.mean(np.hypot(x[ok] - p[0], y[ok] - p[1]) <= NEAR_CM / cm))
        row[f"{o}_x"], row[f"{o}_y"] = core.rnd(p[0], 1), core.rnd(p[1], 1)
        row[f"{o}_near_median"] = core.rnd(near)
        if near < MIN_NEAR:
            flags.append("object_jumps")
        if not judge:
            continue
        if not regions:
            flags.append("off_object")
            continue
        signed = [cv2.pointPolygonTest(r["contour"], (float(p[0]), float(p[1])), True) for r in regions]
        j = int(np.argmax(signed))
        r, dist = regions[j], max(0.0, -signed[j])
        match[o] = j
        row[f"{o}_on_region"] = signed[j] >= 0
        row[f"{o}_region"] = j + 1
        row[f"{o}_to_region_{unit}"] = core.rnd(to_unit(dist), 2)
        row[f"{o}_dx_{unit}"] = core.rnd(to_unit(p[0] - r["centroid"][0]), 2)
        row[f"{o}_dy_{unit}"] = core.rnd(to_unit(p[1] - r["centroid"][1]), 2)
        row[f"{o}_region_area_{unit}2"] = core.rnd(r["area"] * scale * scale if scale else r["area"], 1)
        if dist > OFF_CM / cm:
            flags.append("off_object")
            del match[o]

    present = [o for o in objs if o in pts]
    for i, a in enumerate(present):
        for b in present[i + 1:]:
            row[f"{a}_{b}_sep_{unit}"] = core.rnd(to_unit(np.hypot(*(pts[a] - pts[b]))), 2)
            if a in match and b in match:
                ca, cb = regions[match[a]]["centroid"], regions[match[b]]["centroid"]
                row[f"{a}_{b}_region_sep_{unit}"] = core.rnd(to_unit(np.hypot(*(ca - cb))), 2)
                if match[a] == match[b]:
                    flags.append("objects_merged")
    row["flags"] = list(dict.fromkeys(flags))
    if getattr(opts, "stills", False) and bg is not None:
        img = draw(bg, mask, regions, fx, fy, track, pts, match, row, unit, cm)
        row["_still"] = cv2.imencode(".jpg", img, [cv2.IMWRITE_JPEG_QUALITY, 92])[1].tobytes()
    return row


def find_regions(bg, poly, fx, fy, cm):
    """Regions of the background that differ from the floor, as contours in video pixels (largest first), and the
    floor mask in working pixels."""
    h, w = bg.shape[:2]
    mask = np.zeros((h, w), np.uint8)
    if poly is not None:
        cv2.fillPoly(mask, [np.round(poly * [fx, fy]).astype(np.int32)], 1)
    else:
        mask[:] = 1
    size = np.sqrt(mask.sum())
    mask = cv2.erode(mask, disk(MARGIN * size))
    inside = mask > 0
    if inside.sum() < 100:
        return [], mask
    cm_work = cm / fx
    lab = cv2.cvtColor(bg.astype(np.float32) / 255, cv2.COLOR_BGR2Lab)
    colour = np.sqrt(((lab - np.median(lab[inside], axis=0)) ** 2).sum(axis=2))
    win = max(3, int(round(TEXTURE_CM / cm_work)) | 1)
    light = lab[..., 0]
    mean = cv2.blur(light, (win, win))
    texture = np.sqrt(np.maximum(cv2.blur(light * light, (win, win)) - mean * mean, 0))
    differs = (colour > limit(colour[inside], MIN_COLOUR)) | (texture > limit(texture[inside], MIN_TEXTURE))
    m = (differs & inside).astype(np.uint8)
    m = cv2.morphologyEx(m, cv2.MORPH_OPEN, disk(OPEN_CM / cm_work))
    m = cv2.morphologyEx(m, cv2.MORPH_CLOSE, disk(CLOSE_CM / cm_work))
    # floor patches enclosed by regions and the outline belong to an object the outline cuts (a cup under a lid that
    # reaches past it); the open floor is far larger than HOLE_SHARE. The outline counts as a band as wide as the
    # closing, so a region ending just short of it still encloses its patch.
    edge = mask - cv2.erode(mask, disk(CLOSE_CM / cm_work))
    patches = ((mask > 0) & (m == 0) & (edge == 0)).astype(np.uint8)
    n, labels, stats, _ = cv2.connectedComponentsWithStats(patches, connectivity=4)
    small = [i for i in range(1, n) if stats[i, cv2.CC_STAT_AREA] < HOLE_SHARE * inside.sum()]
    m = (m | np.isin(labels, small)).astype(np.uint8)
    contours, _ = cv2.findContours(m, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    regions = []
    for c in contours:
        c = c.astype(np.float32) / np.float32([fx, fy])
        area = cv2.contourArea(c)
        if area * cm * cm < MIN_AREA_CM2:
            continue
        mo = cv2.moments(c)
        regions.append({"contour": c, "area": area, "centroid": np.array([mo["m10"] / mo["m00"], mo["m01"] / mo["m00"]])})
    regions.sort(key=lambda r: -r["area"])
    return regions, mask


def limit(values, least):
    med = np.median(values)
    return max(least, med + ROBUST_SD * 1.4826 * np.median(np.abs(values - med)))


def disk(radius):
    r = max(1, int(round(radius)))
    return cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * r + 1, 2 * r + 1))


def finish(rows, opts):
    """Group-level flag: the points sit on the objects differently in one group than in the others."""
    stills = {id(r): r.pop("_still", None) for r in rows}
    offset_cm = getattr(opts, "object_offset_cm", OFFSET_CM)
    px_per_cm = getattr(opts, "object_px_per_cm", PX_PER_CM)
    keys = {}
    for r in rows:
        for k in r:
            m = re.fullmatch(r"(.+)_dx_(cm|px)", k)
            if m:
                keys.setdefault((m.group(1), m.group(2)), None)
    for o, unit in keys:
        dx, dy = f"{o}_dx_{unit}", f"{o}_dy_{unit}"
        by_group = {}                                   # (object set, group) -> offsets
        for r in rows:
            if r.get("group") and r.get(dx) not in (None, "") and not NOT_COUNTED & set(r.get("flags", [])):
                by_group.setdefault((r.get("object_set", ""), r["group"]), []).append((float(r[dx]), float(r[dy])))
        med = {g: np.median(np.array(v), axis=0) for g, v in by_group.items()}
        lim = offset_cm if unit == "cm" else offset_cm * px_per_cm
        shift = {}
        for g in med:
            others = [med[h] for h in med if h != g and h[0] == g[0]]
            if others:
                shift[g] = float(np.hypot(*(med[g] - np.median(np.array(others), axis=0))))
        for r in rows:
            g = (r.get("object_set", ""), r.get("group"))
            if g in shift and dx in r:
                r[f"{o}_group_shift_{unit}"] = core.rnd(shift[g], 2)
                flags = r.setdefault("flags", [])
                if shift[g] > lim and "group_offset" not in flags:
                    flags.append("group_offset")
    if getattr(opts, "stills", False):
        out = os.path.join(opts.stills_dir, "objects")
        os.makedirs(out, exist_ok=True)
        for r in rows:
            jpg = stills.get(id(r))
            if r.get("flags") and jpg is not None:
                img = cv2.imdecode(np.frombuffer(jpg, np.uint8), cv2.IMREAD_COLOR)
                bar = np.zeros((30, img.shape[1], 3), np.uint8)
                text = f"{r['id']}  {r.get('group', '')}  " + " ".join(r["flags"])
                cv2.putText(bar, text, (8, 21), cv2.FONT_HERSHEY_SIMPLEX, 0.55, (255, 255, 255), 1, cv2.LINE_AA)
                name = re.sub(r"[^\w.-]+", "_", str(r["id"])) + ".jpg"
                cv2.imwrite(os.path.join(out, name), np.vstack([bar, img]), [cv2.IMWRITE_JPEG_QUALITY, 88])


def draw(bg, mask, regions, fx, fy, track, pts, match, row, unit, cm):
    """The whole view on the left, a zoom on each object on the right."""
    s = np.float32([fx, fy])
    view = bg.copy()
    edge, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    cv2.drawContours(view, edge, -1, (255, 255, 255), 1)
    for r in regions:
        cv2.drawContours(view, [np.round(r["contour"] * s).astype(np.int32)], -1, (0, 0, 255), 1)
    for k, (o, p) in enumerate(pts.items()):
        col = COLOURS[k % len(COLOURS)]
        x, y = track[o]
        step = max(1, len(x) // 600)
        for xi, yi in zip(x[::step], y[::step]):
            if np.isfinite(xi) and np.isfinite(yi):
                cv2.circle(view, (int(xi * fx), int(yi * fy)), 1, col, -1)
        cv2.drawMarker(view, (int(p[0] * fx), int(p[1] * fy)), col, cv2.MARKER_CROSS, 14, 1)
        cv2.putText(view, o, (int(p[0] * fx) + 8, int(p[1] * fy) - 8), cv2.FONT_HERSHEY_SIMPLEX, 0.45, col, 1, cv2.LINE_AA)
    h = view.shape[0]
    tile = h // max(2, len(pts))
    column = np.zeros((h, tile, 3), np.uint8)
    for k, (o, p) in enumerate(pts.items()):
        col = COLOURS[k % len(COLOURS)]
        half = 5 / cm
        if o in match:
            x0, y0, bw, bh = cv2.boundingRect(regions[match[o]]["contour"])
            half = min(12 / cm, max(half, abs(x0 - p[0]) + 1 / cm, abs(x0 + bw - p[0]) + 1 / cm,
                                    abs(y0 - p[1]) + 1 / cm, abs(y0 + bh - p[1]) + 1 / cm))
        z = tile / (2 * half)                      # tile pixels per video pixel
        origin = p - half
        warp = np.float32([[z / fx, 0, -origin[0] * z], [0, z / fy, -origin[1] * z]])
        t = cv2.warpAffine(bg, warp, (tile, tile), flags=cv2.INTER_CUBIC, borderValue=(0, 0, 0))
        to_tile = lambda q: (int(round((q[0] - origin[0]) * z)), int(round((q[1] - origin[1]) * z)))
        if o in match:
            r = regions[match[o]]
            cv2.drawContours(t, [np.round((r["contour"] - origin.astype(np.float32)) * z).astype(np.int32)], -1, (0, 0, 255), 1)
            cv2.circle(t, to_tile(r["centroid"]), 4, (0, 0, 255), 1)
        x, y = track[o]
        step = max(1, len(x) // 300)
        for xi, yi in zip(x[::step], y[::step]):
            if np.isfinite(xi) and np.isfinite(yi):
                cv2.circle(t, to_tile((xi, yi)), 1, col, -1)
        cv2.drawMarker(t, to_tile(p), col, cv2.MARKER_CROSS, 16, 2)
        info = f"{o}  dx {row.get(f'{o}_dx_{unit}', '')} dy {row.get(f'{o}_dy_{unit}', '')} {unit}"
        cv2.putText(t, info, (6, 18), cv2.FONT_HERSHEY_SIMPLEX, 0.45, (255, 255, 255), 1, cv2.LINE_AA)
        column[k * tile:(k + 1) * tile] = t
    return np.hstack([view, column])
