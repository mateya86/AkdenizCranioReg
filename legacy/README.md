# AkdenizCranioReg — Legacy

**Three-phase fiducial-based rigid registration tool for cranial STL models**

[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Developed at **Akdeniz University – Department of Neurosurgery**, Antalya, Turkey.

---

> **This is the legacy (single-POI) release of AkdenizCranioReg.**
> It is preserved here for reference and reproducibility.
> The current maintained version — which supports multiple points of interest and full X, Y, Z output — is available at the [root of the repository](../).

---

## What is different in this version?

This release measures a **single point of interest (POI)** per session and outputs only **X (anteroposterior)** and **Z (superoinferior)** distances relative to the reference point — i.e. it defines the lateral projection of brain structures.

| Feature | This version (legacy) | Current version |
|---------|----------------------|-----------------|
| Points of interest | Single (REF + POI) | Multiple (REF + POI 1…N) |
| Output axes | X, Z only | X, Y, Z per POI |
| HUD | Two-line X and Z display | Table of X/Y/Z per point |

---

## Purpose

AkdenizCranioReg is an interactive Python application that performs **landmark-based rigid registration** of cranial surface meshes (STL format) to standardized anatomical reference planes.

It supports reproducible allocation of the lateral projection of brain structures, surgical targets, and anatomical landmarks in neurosurgical planning and research.

---

## Workflow

The registration follows a strict **three-phase procedure**:

1. **Frankfurt horizontal plane** alignment
   using the highest points of the left and right ear canals (LEAC, REAC) and the lowest points on the orbital rims (LMORB, RMORB)

2. **Midsagittal plane** alignment
   using Nasion and PML (Posterior Midline Fiducial)

3. **Distance measurement**
   from a user-defined reference point (REF) to a single point of interest (POI)
   outputs **X (anteroposterior)** and **Z (superoinferior)** distances

---

## Features

- Interactive 3D landmark picking with visual feedback (PyVista)
- Frankfurt plane alignment with RMSE reporting
- Automatic midsagittal plane alignment using Nasion–PML
- Single point of interest measurement per session
- Color-coded measurement lines (black, yellow, blue)
- Clean two-line HUD displaying X and Z distances + RMSE
- Single-file executable script (no complex installation)

**Output:**
- **Distance X** — anteroposterior distance from REF to POI
- **Distance Z** — superoinferior distance from REF to POI
- **Frankfurt alignment RMSE** — goodness of fit of the Frankfurt plane

---

## Landmark Definitions

| Landmark | Description | Side | Note |
|----------|-------------|------|------|
| LEAC | Highest point (roof) of left external auditory canal | Left | Fixed anatomical landmark |
| REAC | Highest point (roof) of right external auditory canal | Right | Fixed anatomical landmark |
| LMORB | Lowest point on left orbital rim (orbitale) | Left | Fixed anatomical landmark |
| RMORB | Lowest point on right orbital rim (orbitale) | Right | Fixed anatomical landmark |
| Nasion | Midline point at frontonasal suture | Midline | Fixed anatomical landmark |
| PML | Posterior Midline Fiducial | Midline | e.g. foramen of Monro |
| REF | Reference point (user-defined) | — | e.g. head of mandible (HOM) |
| POI | Single point of interest | — | e.g. foramen of Monro |

---

## Setup & Installation

```bash
git clone https://github.com/mateya86/AkdenizCranioReg.git
cd AkdenizCranioReg/legacy
```

Create a virtual environment (recommended):

```bash
python -m venv venv

# Linux / macOS
source venv/bin/activate

# Windows
venv\Scripts\activate
```

Install dependencies:

```bash
pip install numpy scipy pyvista
```

---

## Usage

```bash
python akdeniz_cranio_reg_legacy.py
```

1. Click **"Open STL and Register"**
2. **Phase 1 – Frankfurt plane**
   Pick in order: LEAC → REAC → LMORB → RMORB
3. **Phase 2 – Midsagittal alignment**
   Pick: Nasion → PML (Posterior Midline Fiducial)
4. **Phase 3 – Measurement**
   Pick: REF (e.g. head of mandible) → POI (e.g. foramen of Monro)
5. Registered model with planes and measurement lines appears — Distance X, Distance Z, and RMSE shown in HUD

---

## Demo Video

[Click here to watch the demo video](https://github.com/mateya86/AkdenizCranioReg/raw/main/legacy/screenshots/AkdenizCranioReg_Demo.mp4)

*Demonstrates fiducial extraction, landmark picking (Frankfurt + Midsagittal), alignment, and final X/Z distance measurement for a single point of interest.*

---

## Screenshots

### Phase 1 — Frankfurt Plane Alignment
*Picking LEAC, REAC, LMORB, and RMORB*

![Phase 1](https://raw.githubusercontent.com/mateya86/AkdenizCranioReg/main/legacy/screenshots/phase1_picking.png)

### Phase 2 — Midsagittal Alignment
*Picking Nasion and PML*

![Phase 2](https://raw.githubusercontent.com/mateya86/AkdenizCranioReg/main/legacy/screenshots/phase2_picking.png)

### Phase 3 — Single POI Measurement
*Picking Reference point (REF) and Point of Interest (POI)*

![Phase 3](https://raw.githubusercontent.com/mateya86/AkdenizCranioReg/main/legacy/screenshots/phase3_picking.png)

### Final Result
*Registered model with Frankfurt plane (gold), midsagittal plane (green), measurement lines, and HUD*

![Final Result](https://raw.githubusercontent.com/mateya86/AkdenizCranioReg/main/legacy/screenshots/final_result.png)

---

## Segmentation and Registration Quality Validation (RMSE)

The tool calculates the Frankfurt Alignment RMSE using:

```
RMSE = sqrt[ sum (zi - z_mean)^2 / n ]
```

Where:
- `zi` = vertical position of each Frankfurt point after alignment
- `z_mean` = mean vertical height (centroid) of all 4 points
- `n` = number of points (4)

In clinical practice, an **RMSE < 2.0 mm** is generally considered an excellent fit, accounting for natural cranial asymmetry and picking precision. If your RMSE is high, re-pick the Frankfurt landmarks ensuring they are properly seated on the bony surface.

---

## License

MIT License — see the [LICENSE](../LICENSE) file for details.

---

## Contact

**Muhammad Ateya**
Akdeniz University – Department of Neurosurgery, Antalya, Turkey

- GitHub: [mateya86](https://github.com/mateya86)
- Email: [tc.muhammad86@gmail.com](mailto:tc.muhammad86@gmail.com)
