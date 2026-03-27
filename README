# AkdenizCranioReg

**Three-phase fiducial-based rigid registration tool for cranial STL models**

[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Developed at **Akdeniz University – Department of Neurosurgery**, Antalya, Turkey.

---

## Purpose

AkdenizCranioReg is an interactive Python application that performs **landmark-based rigid registration** of cranial surface meshes (STL format) to standardized anatomical reference planes.

It supports reproducible measurement of the lateral projection of brain structures, surgical targets, and anatomical landmarks in neurosurgical planning and research.

---

## Workflow

The registration follows a strict **three-phase procedure**:

1. **Frankfurt horizontal plane** alignment
   using the highest points of the left and right ear canals (LEAC, REAC) and the lowest points on the orbital rims (LMORB, RMORB)

2. **Midsagittal plane** alignment
   using Nasion and PML (Posterior Midline Fiducial)

3. **Distance measurement**
   Enter the number of POIs → pick REF first, then POI 1 … POI N sequentially
   Outputs **X (anteroposterior)**, **Y (left–right)**, and **Z (superoinferior)** distances per POI relative to REF

---

## Features

- Interactive 3D landmark picking with visual feedback (PyVista)
- Frankfurt plane alignment with RMSE reporting
- Automatic midsagittal (YZ) plane alignment using Nasion–PML
- **Multiple points of interest** measured in a single session
- Color-coded measurement lines per POI
- HUD table displaying X, Y, Z per point + axis color legend
- Single-file executable script (no complex installation)

**Output per POI:**
- **X** — anteroposterior distance (anterior = +, posterior = −)
- **Y** — left–right distance (right = +, left = −)
- **Z** — superoinferior distance (superior = +, inferior = −)
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
| PML | Posterior Midline Fiducial | Midline | e.g. posterior fossa midline point |
| REF | Reference point (user-defined) | — | e.g. head of mandible (HOM) |
| POI 1…N | Points of interest | — | Multiple targets measured in one session |

---

## Setup & Installation

```bash
git clone https://github.com/mateya86/AkdenizCranioReg.git
cd AkdenizCranioReg
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
pip install -r requirements.txt
```

**Requirements:** `numpy`, `scipy`, `pyvista`

---

## Usage

```bash
python akdeniz_cranio_reg.py
```

1. Click **"Open STL and Register"**
2. **Phase 1 – Frankfurt plane**
   Pick in order: LEAC → REAC → LMORB → RMORB
3. **Phase 2 – Midsagittal alignment**
   Pick: Nasion → PML (Posterior Midline Fiducial)
4. Enter the **number of POIs** when prompted
5. **Phase 3 – Measurement**
   Pick: REF → Point1 → Point2 … → PointN
6. Registered model with all measurement lines appears — X, Y, Z per point shown in HUD

---

## Demo Video

[Click here to watch the demo video](https://github.com/mateya86/AkdenizCranioReg/raw/main/screenshots/AkdenizCranioReg_Demo.mp4)

*Demonstrates fiducial extraction, landmark picking (Frankfurt + Midsagittal), alignment, and final X/Y/Z distance measurement for multiple points of interest.*

---

## Screenshots

### Phase 1 — Frankfurt Plane Alignment
*Picking LEAC, REAC, LMORB, and RMORB*

![Phase 1](https://raw.githubusercontent.com/mateya86/AkdenizCranioReg/main/screenshots/phase1_picking.png)

### Phase 2 — Midsagittal Alignment
*Picking Nasion and PML*

![Phase 2](https://raw.githubusercontent.com/mateya86/AkdenizCranioReg/main/screenshots/phase2_picking.png)

### Phase 3 — Multi-POI Measurement
*Picking REF followed by multiple Points of Interest*

![Phase 3](https://raw.githubusercontent.com/mateya86/AkdenizCranioReg/main/screenshots/phase3_picking.png)

### Final Result
*Registered model with Frankfurt plane (gold), midsagittal plane (green), color-coded measurement lines per POI, and HUD*

![Final Result](https://raw.githubusercontent.com/mateya86/AkdenizCranioReg/main/screenshots/final_result.png)

---

## Test Data

Two sample cranial STL models are provided in the `test_data/` folder for validation and demonstration purposes.

| File | REF point | POI |
|------|-----------|-----|
| `ref point left HOM, poi=FOM.stl` | Left Head of Mandible (HOM) | Foramen of Monro (FOM) |
| `ref point right HOM, poi=FOM.stl` | Right Head of Mandible (HOM) | Foramen of Monro (FOM) |

These models allow users to reproduce the full three-phase registration workflow and verify measurement outputs before applying the tool to their own data.

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

## Repository Structure

```
AkdenizCranioReg/
├── akdeniz_cranio_reg.py       ← main script (multi-POI)
├── README.md
├── requirements.txt
├── LICENSE
├── screenshots/                ← demo video and phase screenshots
└── test_data/                  ← sample STL models
```

---

## Citation

If you use this tool in your research or clinical work, please cite:

```bibtex
@software{ateya_akdenizcranioreg_2026,
  author       = {Ateya, Muhammad},
  title        = {AkdenizCranioReg: Three-phase cranial STL registration tool},
  year         = {2026},
  publisher    = {Akdeniz University - Department of Neurosurgery},
  url          = {https://github.com/mateya86/AkdenizCranioReg},
  howpublished = {Computer software}
}
```

---

## License

MIT License — see the [LICENSE](LICENSE) file for details.

---

## Contact

**Muhammad Ateya**
Akdeniz University – Department of Neurosurgery, Antalya, Turkey

Feel free to open an issue or contact via GitHub for questions, feature requests, or collaboration.
