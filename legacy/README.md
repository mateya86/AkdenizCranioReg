# AkdenizCranioReg

**Three-phase fiducial-based rigid registration tool for cranial STL models**

[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Developed at **Akdeniz University – Department of Neurosurgery**, Antalya, Turkey.

## Purpose

AkdenizCranioReg is an interactive Python application that performs **landmark-based rigid registration** of cranial surface meshes (STL format) to standardized anatomical reference planes.

It supports reproducible measurement of the lateral projection of brain structures, surgical targets and anatomical landmarks in neurosurgical planning and research.

## Workflow

The registration follows a strict **three-phase procedure**:

1. **Frankfurt horizontal plane** alignment  
   using the highest points of the left and right ear canals (LEAC, REAC) and the lowest points on the orbital rims (LMORB, RMORB)

2. **Midsagittal plane** alignment  
   using nasion and PML (Posterior Midline Fiducial)

3. **Distance measurement**  
   from a user-defined reference point (REF) to a point of interest (POI)  
   → outputs perpendicular **X (anteroposterior)** and **Z (superoinferior)** distances

## Features

- Interactive 3D landmark picking with visual feedback (PyVista)
- Frankfurt plane alignment with RMSE reporting
- Automatic midsagittal (YZ) alignment
- Anatomical coordinate system:  
  X = A/P (blue), Y = L/R (black), Z = S/I (yellow)
- Measurement visualization with colored lines
- Clean HUD showing X and Z distances + RMSE
- Single-file executable script (no complex installation)

## Setup & Usage

```bash
git clone https://github.com/mateya86/AkdenizCranioReg.git
cd AkdenizCranioReg

# Create virtual environment (recommended)
python -m venv venv

# Linux / macOS
source venv/bin/activate

# Windows
venv\Scripts\activate

pip install -r requirements.txt

# Run the tool
python akdeniz_cranio_reg.py

1. Click "Open STL and Register"

2. Follow the guided picking phases:  
   Phase 1 – Frankfurt plane  
   Pick in order: LEAC → REAC → LMORB → RMORB  
   Phase 2 – Midsagittal alignment  
   Pick: Nasion → PML (Posterior Midline Fiducial)  
   Phase 3 – Measurement  
   Pick: REF (reference point, e.g. head of mandible / HOM) → POI (point of interest, e.g. foramen of Monro or any structure you want to measure relative to REF)

3. After picking, the registered model with planes and measurement lines appears.

## Landmark Definitions

| Landmark | Description                                          | Side    | Note / typical example                     |
|----------|------------------------------------------------------|---------|--------------------------------------------|
| LEAC     | Highest point (roof) of left external auditory canal | Left    | Fixed anatomical landmark                  |
| REAC     | Highest point (roof) of right external auditory canal| Right   | Fixed anatomical landmark                  |
| LMORB    | Lowest point on left orbital rim (orbitale)          | Left    | Fixed anatomical landmark                  |
| RMORB    | Lowest point on right orbital rim (orbitale)         | Right   | Fixed anatomical landmark                  |
| Nasion   | Midline point at frontonasal suture                  | Midline | Fixed anatomical landmark                  |
| PML      | Posterior Midline Fiducial                           | Midline | Fixed anatomical landmark (e.g. foramen of Monro) |
| REF      | Reference point (user-defined)                       | —       | e.g. head of mandible (HOM)                |
| POI      | Point of interest (structure being measured)         | —       | e.g. foramen of Monro                      |

## Coordinate System & Outputs

After registration:

- **X (blue)** → Anteroposterior (positive = anterior)
- **Y (black)** → Left–Right (positive = right)
- **Z (yellow)** → Superoinferior (positive = superior)

**Displayed in HUD:**
- Distance X (A/P)
- Distance Z (S/I)
- Frankfurt alignment RMSE

## Demo Video

Watch a quick 1-minute demonstration of the full three-phase workflow:

<video width="800" controls>
  <source src="screenshots/AkdenizCranioReg_Demo.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>

*Demonstrates fiducial extraction, landmark picking (Frankfurt + Midsagittal), alignment, and final X/Z distance measurement with visual feedback.*

## Screenshots

Phase 1 — Picking LEAC, REAC, LMORB, and RMORB for Frankfurt horizontal plane alignment

Phase 2 — Picking Nasion and PML (Posterior Midline Fiducial)

Phase 3 — Picking Reference point (REF) and Point of Interest (POI)

Final registered model showing Frankfurt plane (gold), midsagittal plane (green), measurement lines (blue, black, yellow), and HUD with X/Z distances and RMSE

## Data Validation (RMSE)

The tool calculates the Frankfurt Alignment RMSE using the following logic:

RMSE = √[ Σ (zᵢ - z̄)² / n ]

Where:

√ = Square root of the entire result.

Σ = The sum of all values.

zᵢ = The vertical position of each individual Frankfurt point.

z̄ (z-bar) = The average vertical height (centroid) of all 4 points.

n = The number of points (in Phase 1, n = 4).

In clinical practice, an RMSE < 2.0 mm is generally considered an excellent fit, accounting for natural cranial asymmetry and picking precision. If your RMSE is high, we recommend re-picking the Frankfurt landmarks (LEAC, REAC, LMORB, RMORB) to ensure they are properly seated on the bony surface.


## Citation

@software{ateya_akdenizcranioreg_2026,
  author       = {Ateya, Muhammad},
  title        = {AkdenizCranioReg: Three-phase cranial STL registration tool},
  year         = {2026},
  publisher    = {Akdeniz University – Department of Neurosurgery},
  url          = {https://github.com/mateya86/AkdenizCranioReg},
  howpublished = {Computer software}
}

## License

MIT License — see the LICENSE file.

##Contact

Muhammad Ateya
Akdeniz University – Department of Neurosurgery
Antalya, Turkey
Feel free to open an issue or contact me for questions, feature requests or collaboration.
