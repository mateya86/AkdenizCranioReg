"""
AkdenizCranioReg
----------------
Three-phase fiducial-based registration tool for cranial STL models.
Aligns skull to Frankfurt plane and midsagittal plane using anatomical landmarks,
then computes distances from a reference point to a point of interest.

Developed by: Muhammad Ateya  
Department of Neurosurgery, Akdeniz University, Antalya, Turkey

Usage:
    python akdeniz_cranio_reg.py

Requirements:
    numpy, scipy, pyvista, tkinter (usually comes with Python)

License: MIT (see LICENSE file in repository)
"""

import numpy as np
from scipy.spatial.transform import Rotation as R
import pyvista as pv
import tkinter as tk
from tkinter import filedialog, messagebox

# -----------------------------
# Math / transformation helpers
# -----------------------------
def rotation_to_4x4(rot: R) -> np.ndarray:
    M = np.eye(4, dtype=np.float64)
    M[:3, :3] = rot.as_matrix()
    return M

def apply_rotation(mesh: pv.DataSet, rot: R):
    m = mesh.copy()
    M = rotation_to_4x4(rot)
    m.transform(M, inplace=True)
    return m

# -----------------------------
# Registration helpers
# -----------------------------
def align_frankfurt(points):
    try:
        p1, p2, p3, p4 = points["LEAC"], points["REAC"], points["LMORB"], points["RMORB"]
        centroid = np.mean([p1, p2, p3, p4], axis=0)
        centered_points = np.array([p1, p2, p3, p4]) - centroid
        _, _, V = np.linalg.svd(centered_points)
        normal = V[-1]
        nrm = np.linalg.norm(normal)
        if nrm < 1e-10:
            raise ValueError("Frankfurt normal undefined: points are collinear or coincident.")
        normal = normal / nrm
        if normal[2] < 0:
            normal = -normal
        target = np.array([0.0, 0.0, 1.0], dtype=np.float64)
        rot1, _ = R.align_vectors([target], [normal])
        if np.linalg.det(rot1.as_matrix()) < 0:
            rot1 = rot1 * R.from_rotvec(np.array([0, 0, np.pi]))
        pts_rot1 = {k: rot1.apply(v) for k, v in points.items()}
        return pts_rot1, rot1
    except Exception as e:
        print(f"Frankfurt alignment error: {e}")
        raise

def align_sagittal_from_two(points_rot1, key_nasion="Nasion", key_pml="PML"):
    try:
        nasion, pml = points_rot1[key_nasion], points_rot1[key_pml]
        midline = pml - nasion
        mlen = np.linalg.norm(midline)
        if mlen < 1e-10:
            raise ValueError("Midline undefined: Nasion and PML coincide.")
        midline = midline / mlen
        proj = np.array([midline[0], midline[1], 0.0], dtype=np.float64)
        plen = np.linalg.norm(proj)
        if plen < 1e-10:
            raise ValueError("Midline projection onto XY is zero (vector parallel to Z).")
        proj = proj / plen
        target = np.array([0.0, 1.0, 0.0], dtype=np.float64)
        rot2, _ = R.align_vectors([target], [proj])
        if np.linalg.det(rot2.as_matrix()) < 0:
            rot2 = rot2 * R.from_rotvec(np.array([0, 0, np.pi]))
        pts_aligned = {k: rot2.apply(v) for k, v in points_rot1.items()}
        return pts_aligned, rot2
    except Exception as e:
        print(f"Sagittal alignment error: {e}")
        raise

def compose_rotations(rot1: R, rot2: R) -> np.ndarray:
    try:
        M = np.eye(4, dtype=np.float64)
        M[:3, :3] = rot2.as_matrix() @ rot1.as_matrix()
        return M
    except Exception as e:
        print(f"Rotation composition error: {e}")
        raise

# -----------------------------
# Plane creation helper
# -----------------------------
def create_plane(point, normal, size=175):
    normal = np.array(normal, dtype=np.float64)
    nrm = np.linalg.norm(normal)
    if nrm < 1e-10:
        raise ValueError("Invalid plane normal: zero magnitude.")
    normal = normal / nrm
    if abs(normal[2]) > 0.99:
        u = np.cross(normal, [1, 0, 0])
        if np.linalg.norm(u) < 1e-10:
            u = np.cross(normal, [0, 1, 0])
    else:
        u = np.cross(normal, [0, 0, 1])
    u = u / np.linalg.norm(u)
    v = np.cross(normal, u)
    v = v / np.linalg.norm(v)
    half_size = size / 2
    vertices = [
        point + u * half_size + v * half_size,
        point + u * half_size - v * half_size,
        point - u * half_size - v * half_size,
        point - u * half_size + v * half_size
    ]
    faces = [4, 0, 1, 2, 3]
    return pv.PolyData(vertices, faces)

# -----------------------------
# Interactive picking
# -----------------------------
def pick_points(mesh: pv.DataSet, labels, message):
    try:
        plotter = pv.Plotter()
        plotter.set_background("white")
        plotter.add_mesh(mesh, color="lightblue", opacity=0.6, show_edges=True)
        text_actor = plotter.add_text(
            message, font_size=10, font="arial", position="upper_edge", color="black"
        )
        picked_points = []
        colors = ["red", "orange", "cyan", "green", "blue", "magenta"][:len(labels)]
        picked_count = [0]

        def extract_coords(*args):
            for a in args:
                if isinstance(a, (list, tuple, np.ndarray)) and len(a) == 3:
                    return np.array(a, dtype=np.float64)
                if hasattr(a, "points"):
                    try:
                        p = np.array(a.points[0], dtype=np.float64)
                        if p.shape == (3,):
                            return p
                    except Exception:
                        pass
            return None

        def on_pick(*args):
            nonlocal picked_count
            point = extract_coords(*args)
            if point is None:
                return
            if picked_count[0] < len(labels):
                picked_points.append(point)
                idx = picked_count[0]
                print(f"Picked point {labels[idx]}: {point}")
                plotter.add_mesh(pv.Sphere(radius=0.1, center=point), color=colors[idx])
                plotter.add_point_labels([point], [labels[idx]], font_size=10, text_color=colors[idx])
                picked_count[0] += 1
                if picked_count[0] < len(labels):
                    next_label = labels[picked_count[0]]
                    text_actor.SetText(0, f"Next:\nPick {next_label}")
                else:
                    text_actor.SetText(0, "Picking complete.\nClose the window.")
                plotter.render()
                plotter.update()

        plotter.enable_point_picking(
            callback=on_pick,
            use_picker=True,
            show_message=True,
            show_point=True,
            font_size=10
        )
        plotter.show()
        if len(picked_points) != len(labels):
            raise ValueError(f"Expected {len(labels)} points, but picked {len(picked_points)}.")
        return {label: point for label, point in zip(labels, picked_points)}
    except Exception as e:
        print(f"Picking error: {e}")
        raise

# -----------------------------
# Main three-phase workflow
# -----------------------------
def run_three_phase_registration(stl_path: str):
    try:
        print(f"Loading STL: {stl_path}")
        mesh_orig = pv.read(stl_path)
        mesh_orig = mesh_orig.triangulate().compute_normals(auto_orient_normals=True)

        print("Starting Phase 1: Frankfurt picking")
        frankfurt_labels = ["LEAC", "REAC", "LMORB", "RMORB"]
        frankfurt_points = pick_points(
            mesh_orig,
            frankfurt_labels,
            "Phase 1:\nPick LEAC, then REAC,\nthen LMORB, then RMORB\n(click on mesh surface,\nthen close window)"
        )
        pts_rot1, rot1 = align_frankfurt(frankfurt_points)

        mesh_frankfurt = apply_rotation(mesh_orig, rot1)

        print("Starting Phase 2: Sagittal picking (2 points)")
        sagittal_labels = ["Nasion", "PML"]
        sagittal_points = pick_points(
            mesh_frankfurt,
            sagittal_labels,
            "Phase 2:\nPick Nasion, then PML\n(Posterior Midline Fiducial)\n(click on mesh surface,\nthen close window)"
        )
        pts_rot1.update(sagittal_points)
        pts_aligned, rot2 = align_sagittal_from_two(pts_rot1, key_nasion="Nasion", key_pml="PML")
        M_total = compose_rotations(rot1, rot2)

        # Frankfurt RMSE calculation
        frankfurt_pts_rot1 = np.array([rot1.apply(frankfurt_points[p])
                                       for p in ["LEAC", "REAC", "LMORB", "RMORB"]])
        z_vals = frankfurt_pts_rot1[:, 2]
        z_centroid = np.mean(z_vals)
        z_dev = z_vals - z_centroid
        rmse_frankfurt = np.sqrt(np.mean(z_dev ** 2))
        print(f"Frankfurt plane alignment RMSE: {rmse_frankfurt:.3f} mm")

        mesh_after_sagittal = mesh_orig.copy()
        mesh_after_sagittal.transform(M_total, inplace=True)

        print("Starting Phase 3: pick Reference and Point of Interest")
        phase3_labels = ["REF", "POI"]
        phase3_points = pick_points(
            mesh_after_sagittal,
            phase3_labels,
            "Phase 3:\nPick Reference point (e.g., left HOM),\nthen Point of Interest (tested point)\n(click on mesh surface,\nthen close window)"
        )
        pts_aligned.update({
            "REF": phase3_points["REF"],
            "POI": phase3_points["POI"]
        })

        leac_reac_vec = pts_aligned["REAC"] - pts_aligned["LEAC"]
        frankfurt_normal = np.array([0.0, 0.0, 1.0], dtype=np.float64)
        plane_normal = np.cross(leac_reac_vec, frankfurt_normal)
        nrm = np.linalg.norm(plane_normal)
        if nrm < 1e-10:
            raise ValueError("LEAC-REAC plane normal undefined: LEAC and REAC coincide.")
        plane_normal = plane_normal / nrm

        X_abs = abs(np.dot(pts_aligned["REF"] - pts_aligned["POI"], plane_normal))
        leac_reac_avg = (pts_aligned["LEAC"] + pts_aligned["REAC"]) / 2
        d_poi = np.dot(pts_aligned["POI"], plane_normal)
        d_ref = np.dot(pts_aligned["REF"], plane_normal)
        d_avg = np.dot(leac_reac_avg, plane_normal)
        dist_poi = abs(d_poi - d_avg)
        dist_ref = abs(d_ref - d_avg)
        X = -X_abs if dist_poi < dist_ref else X_abs

        line_dir = np.cross(plane_normal, frankfurt_normal)
        line_nrm = np.linalg.norm(line_dir)
        if line_nrm < 1e-10:
            raise ValueError("Line direction undefined: plane_normal parallel to Z.")
        line_dir = line_dir / line_nrm

        poi_to_ref = pts_aligned["POI"] - pts_aligned["REF"]
        poi_proj = pts_aligned["POI"] - np.dot(poi_to_ref, plane_normal) * plane_normal
        vec = poi_proj - pts_aligned["REF"]
        Z = np.linalg.norm(vec - np.dot(vec, line_dir) * line_dir)

        mesh_registered = mesh_orig.copy()
        mesh_registered.transform(M_total, inplace=True)
        rot_swap = R.from_rotvec(np.array([0, 0, np.pi/2]))
        M_swap = rotation_to_4x4(rot_swap)
        mesh_registered.transform(M_swap, inplace=True)
        pts_aligned_swapped = {k: rot_swap.apply(v) for k, v in pts_aligned.items()}

        frankfurt_pts_swapped = np.array([
            pts_aligned_swapped["LEAC"],
            pts_aligned_swapped["REAC"],
            pts_aligned_swapped["LMORB"],
            pts_aligned_swapped["RMORB"],
        ], dtype=np.float64)
        frankfurt_centroid = np.mean(frankfurt_pts_swapped, axis=0)
        frankfurt_plane_normal = np.array([0.0, 0.0, 1.0], dtype=np.float64)

        xy = frankfurt_pts_swapped[:, :2]
        center_xy = np.mean(xy, axis=0)
        radial_spans = np.linalg.norm(xy - center_xy, axis=1)
        diameter_xy = 2.0 * np.max(radial_spans) if radial_spans.size > 0 else 175.0
        frankfurt_size = max(175.0, float(diameter_xy + 20.0))
        frankfurt_plane = create_plane(frankfurt_centroid, frankfurt_plane_normal, size=frankfurt_size)

        nasion_sw = pts_aligned_swapped["Nasion"]
        poi_sw = pts_aligned_swapped["POI"]
        ref_sw = pts_aligned_swapped["REF"]

        nf_vec = poi_sw - nasion_sw
        nf_len = np.linalg.norm(nf_vec)
        if nf_len < 1e-10:
            raise ValueError("Nasion–POI plane undefined: points coincide.")
        nf_vec = nf_vec / nf_len
        nf_plane_normal = np.cross(nf_vec, frankfurt_plane_normal)
        n2 = np.linalg.norm(nf_plane_normal)
        if n2 < 1e-10:
            raise ValueError("Nasion–POI plane normal undefined: direction parallel to +Z.")
        nf_plane_normal = nf_plane_normal / n2
        nf_center = 0.5 * (nasion_sw + poi_sw)
        nf_plane = create_plane(nf_center, nf_plane_normal, size=frankfurt_size)

        t_perp = -np.dot(nf_plane_normal, ref_sw - nf_center)
        ref_to_nf_intersection = ref_sw + t_perp * nf_plane_normal
        line1 = pv.Line(ref_sw, ref_to_nf_intersection)

        ref_proj_on_nf = ref_sw - np.dot(nf_plane_normal, ref_sw - nf_center) * nf_plane_normal
        ref_to_poi_zlevel = np.array([ref_proj_on_nf[0], ref_proj_on_nf[1], poi_sw[2]], dtype=np.float64)
        line2_yellow = pv.Line(ref_proj_on_nf, ref_to_poi_zlevel)
        line3_blue = pv.Line(ref_to_poi_zlevel, poi_sw)

        # ───────────────────────────────────────────────
        # FINAL VISUALIZATION
        # ───────────────────────────────────────────────
        print("Starting final visualization")
        plotter = pv.Plotter()
        plotter.set_background("white")
        plotter.add_mesh(mesh_registered, color="lightblue", show_edges=True, opacity=0.75)
        plotter.add_mesh(frankfurt_plane, color="gold", opacity=0.25, show_edges=True, line_width=2)
        plotter.add_mesh(nf_plane, color="lightgreen", opacity=0.25, show_edges=True, line_width=2)
        plotter.add_mesh(line1, color="black", line_width=4)
        plotter.add_mesh(line2_yellow, color="yellow", line_width=6)
        plotter.add_mesh(line3_blue, color="blue", line_width=6)

        for label, coord in pts_aligned_swapped.items():
            plotter.add_mesh(pv.Sphere(radius=1.0, center=coord), color="red")
            plotter.add_point_labels([coord], [label], font_size=10, text_color="black")

        # ── HUD ──
        plotter.add_text(
            f"Distance X : {X:.2f} mm",
            position="upper_left",
            font_size=14, color="blue", font="arial", shadow=True
        )

        plotter.add_text(
            f"Distance Z : {Z:.2f} mm",
            position="upper_right",
            font_size=14, color=[1.0, 0.80, 0.0], font="arial", shadow=True
        )

        plotter.add_text(
            "Plane legend:\n  [gold]   Frankfurt plane\n  [green]  Midsagittal plane",
            position="lower_left",
            font_size=11, color="black", font="arial"
        )

        plotter.add_text(
            f"Frankfurt RMSE:\n  {rmse_frankfurt:.3f} mm",
            position="lower_right",
            font_size=11, color="gray", font="arial", shadow=True
        )

        plotter.add_axes()
        plotter.show_grid()
        plotter.show()

    except Exception as e:
        print(f"Error: {e}")
        raise

# -----------------------------
# GUI
# -----------------------------
def open_stl_and_register():
    filepath = filedialog.askopenfilename(filetypes=[("STL files", "*.stl")])
    if filepath:
        try:
            run_three_phase_registration(filepath)
        except Exception as e:
            messagebox.showerror("Error", str(e))

def main():
    root = tk.Tk()
    root.title("AkdenizCranioReg")
    frame = tk.Frame(root, padx=20, pady=20)
    frame.pack()
    tk.Button(frame, text="Open STL and Register", command=open_stl_and_register, width=30).pack(pady=10)
    tk.Label(
        frame,
        text="AkdenizCranioReg\nThree-Phase Cranial Registration Tool\n"
             "Akdeniz University – Neurosurgery Department, M. Ateya",
        justify="left"
    ).pack(pady=10)
    root.mainloop()

if __name__ == "__main__":
    main()