# SAXS In-Situ Flow-Through Data Reduction & Publication Repository

[![Beamline](https://img.shields.io/badge/Beamline-BESSY%20II%20%C2%B5Spot-blue)](https://www.helmholtz-berlin.de/)
[![Energy](https://img.shields.io/badge/Beam%20Energy-18.00%20keV-orange)]()
[![Detector](https://img.shields.io/badge/Detector-Dectris%20Eiger%209M-green)]()
[![Uncertainty](https://img.shields.io/badge/Uncertainty-7--Component%20Analytical-purple)]()
[![Self--Contained](https://img.shields.io/badge/Data%20Integrity-100%25%20Physical%20(No%20Symlinks)-teal)]()

Publication dataset and reduction pipeline for Small-Angle X-ray Scattering (SAXS) measurements collected at the **BESSY II µSpot (myspot)** beamline (Helmholtz-Zentrum Berlin, March 2026 beamtime).

This repository contains raw detector frames, beamline metadata, geometry calibration files, reduced 1D and 2D profiles, figures, and an HTML report. All files are stored as physical copies without symlinks.

---

## Table of Contents
1. [Experimental & Beamline Parameters](#experimental--beamline-parameters)
2. [Experiment Scans & Flow-Through Beamtime](#experiment-scans--flow-through-beamtime)
3. [Mathematical Data Reduction & Error Propagation](#mathematical-data-reduction--error-propagation)
4. [Directory Structure & Data Dictionary](#directory-structure--data-dictionary)
5. [Reproducibility & Execution Guide](#reproducibility--execution-guide)
6. [Interactive Publication Report](#interactive-publication-report)
7. [Data Governance & Provenance](#data-governance--provenance)

---

## Experimental & Beamline Parameters

| Parameter | Value | Details / Notes |
| :--- | :--- | :--- |
| **Synchrotron Light Source** | **BESSY II** (Berlin, Germany) | Helmholtz-Zentrum Berlin für Materialien und Energie |
| **Beamline** | **µSpot (myspot)** | Microfocus SAXS/WAXS/XRF endstation |
| **X-ray Photon Energy** | **18.00 keV** | Monochromatic beam ($\mathrm{Si}(111)$ double crystal monochromator) |
| **X-ray Wavelength ($\lambda$)** | **$0.0688916\ \mathrm{nm}$** ($0.6889\ \mathrm{Å}$) | Calibrated via `calib.poni` |
| **Detector System** | **Dectris Eiger 9M** | Hybrid photon counting (HPC) active pixel detector |
| **Detector Dimensions** | $3110 \times 3269$ pixels | Active area $\approx 233.25 \times 245.18\ \mathrm{mm}^2$ |
| **Pixel Pitch** | $75.0 \times 75.0\ \mu\mathrm{m}^2$ | Isotropic square pixels |
| **Sample-to-Detector Distance** | **$0.37328\ \mathrm{m}$ ($37.33\ \mathrm{cm}$)** | `calib.poni` parameter `Distance: 0.373278` |
| **Direct Beam Center $(X_0, Y_0)$** | **$(1638.14, 1412.06)\ \mathrm{px}$** | `poni1 = 0.105904 m`, `poni2 = 0.122860 m` |
| **Detector Tilt Angles** | $\mathrm{Rot}_1 = 0.00113\ \mathrm{rad}$, $\mathrm{Rot}_2 = -0.00643\ \mathrm{rad}$ | Planar orthogonal alignment |
| **Absolute Standard Standardizer** | **Glassy Carbon Standard** | $K_{std} = 385.2002 \pm 1.3760\ \mathrm{cm}^{-1}/(\mathrm{cts}/\mathrm{mon})$ |
| **Calibration Reference Standard**| **Quartz Capillary** | Reference run scans 57–58 (`quartz57_58.tiff`) |

---

## Experiment Scans & Flow-Through Beamtime

The in-situ flow-through reaction was monitored continuously over **16.15 hours**, comprising **3,134 distinct acquisition frames** across consecutive SPEC scans:

| Series Identifier | SPEC Scan Range | Primary Focus & Description |
| :--- | :--- | :--- |
| **Exp2_1** | **Scan #184 to #274** | **In-situ kinetic flow-through series (Part 1, 0.0 h to 11.2 h)** |
| **Exp2_2** | **Scan #275 to #292** | **In-situ kinetic flow-through series (Part 2, 11.2 h to 16.15 h)** |
| **Exp2_1 Scan #184 (Frame 1)** | Scan #184, frame 1 | **Primary publication frame** ($t = 0.0\ \mathrm{h}$) |
| **CapillaryTop Scan #63** | Scan #63 (all 11 frames) | **Buffer baseline scan** ($15\ \mathrm{s}$, $I_0 = 973,266$) |

---

## Mathematical Data Reduction & Error Propagation

Data reduction transforms raw 2D pixel count arrays into absolute macroscopic differential scattering cross-sections $\frac{d\Sigma}{d\Omega}(q)\ [\mathrm{cm}^{-1}]$ with analytical variance propagation matching Spiger `integrator20.py`:

### 1. Detector Masking (Gaps, Borders, Dead and Hot Pixels)
To prevent pane border lines and detector defect artifacts, an expanded mask $M(x, y)$ ($787,976$ pixels, $7.75\%$ of sensor) was constructed by extending the Spiger baseline mask:
1. **Beamstop & Baseline Geometry**: Beamstop arm and shadow from `Start_ide.tiff`.
2. **Detector Sentinels**: Overflow pixels ($I_{\text{raw}} \ge 4,294,967,000$).
3. **Dilated Module Pane Gaps**: The 5 horizontal and 2 vertical inter-module gaps of the Dectris Eiger 9M dilated by 2 iterations (`iterations=2`). This covers the pane boundaries and double-sized silicon border pixels (rows 512–552, 1063–1103, 1614–1654, 2165–2205, 2716–2756, and cols 1028–1041, 2068–2081), removing edge effects and pane boundary lines.
4. **Outer Boundary Guard**: 2-pixel margin around the detector perimeter.
5. **Dead Pixels**: 30 pixels with count $\le 0$ or counts $< 25\%$ of local $3 \times 3$ median.
6. **Hot Pixels**: 68 pixels exceeding $7\sigma$ above local $3 \times 3$ median (counts up to $2.3 \times 10^7$).
$$M(x, y) = M_{\text{beamstop}}(x, y) \cup \mathcal{D}_2[M_{\text{gaps}}(x, y)] \cup M_{\text{dead}}(x, y) \cup M_{\text{hot}}(x, y) \cup M_{\text{sat}}(x, y)$$
The resulting mask is saved in [`calibration/Start_rigorous_mask.tiff`](file:///run/media/tomek/data/Mrc_26/saxs_publication_project/calibration/Start_rigorous_mask.tiff).


### 2. Monitor flux normalization and transmission correction (pulse pileup physics)
Incident photon flux fluctuations are scaled using the Keithley ion chamber counter $I_0$. Transmission factors $T$ and relative sample transmission $T_{\text{rel}}$ are computed from transmitted beam diode readings ($\text{roi3}$):
$$T_{\text{sample}} = \frac{\text{roi3}_{\text{sample}} \cdot t_{\text{sec}}}{\text{ltime} \cdot I_{0, \text{sample}}}, \quad T_{\text{rel}} = \frac{T_{\text{sample}}}{T_{\text{bkg}}}$$
$$\sigma_{T_{\text{rel}}} = T_{\text{rel}} \sqrt{\left(\frac{\sigma_{T, \text{sam}}}{T_{\text{sam}}}\right)^2 + \left(\frac{\sigma_{T, \text{bkg}}}{T_{\text{bkg}}}\right)^2}, \quad \frac{\sigma_T}{T} = \sqrt{\frac{1}{\text{roi3}} + 0.02^2}$$

> **Detector physics: ROI3 livetime and pulse pileup correction**  
> $\text{roi3}$ is the integrated intensity of the transmitted beam measured by an energy-dispersive XRF detector / multi-channel analyzer (MCA). At high photon count rates, incoming pulses overlap within the electronic shaping time (pulse pileup), rendering the hardware temporarily dead ($t_{\text{dead}} = t_{\text{sec}} - t_{\text{live}}$). In our beamtime, Scan #184 had $26.68\%$ dead time ($t_{\text{live}} = 7.33\text{ s}$ out of $10\text{ s}$), whereas Background #63 had $37.50\%$ dead time ($t_{\text{live}} = 9.37\text{ s}$ out of $15\text{ s}$) due to higher transmitted flux through the empty capillary. Dividing by $t_{\text{live}}$ ($\text{Rate}_{\text{ROI3}} = \text{roi3} / t_{\text{live}}$) recovers the true incoming photon rate independent of whether the system counted for longer or shorter effective time due to pileup.


### 3. Quantitative reduction: 1D subtraction preceded by integration
Azimuthal radial integration is executed directly on the masked raw 2D detector patterns before subtraction, preserving pure Poisson counting statistics ($\sigma = \sqrt{I}$) and avoiding negative-count artifacts in radial binning:
$$I_{\text{sam}}(q) = \text{integrate1d}(I_{\text{sample,raw}}), \quad I_{\text{bkg}}(q) = \text{integrate1d}(I_{\text{bkg,raw}})$$
Profiles are subsequently normalized by monitor and relative transmission in 1D, and subtracted:
$$I_{\text{sam,norm}}(q) = \frac{I_{\text{sam}}(q)}{I_{0, \text{sample}} \cdot T_{\text{rel}}}, \quad I_{\text{bkg,norm}}(q) = \frac{I_{\text{bkg}}(q)}{I_{0, \text{bkg}}}$$
$$I_{\text{sub}}(q) = K_{\text{std}} \left[ I_{\text{sam,norm}}(q) - I_{\text{bkg,norm}}(q) \right]$$
where $K_{\text{std}} = 385.2002\ \mathrm{cm}^{-1}$ is the glassy carbon standard factor for absolute intensity calibration. (2D subtracted frames $I_{\text{sub}}(x, y)$ are generated strictly for qualitative visual verification of parasitic scatter removal).

#### Empty beamline background protocol
The empty beamline background (air and instrument parasitic scatter without a capillary, Scan #47 `AirGlassy_000047`) is **inherently accounted for within the capillary buffer background** (Scan #63 `CapillaryTop`). The capillary sits directly in the beamline path and captures both borosilicate glass scatter and transmitted empty beamline air scatter: $I_{\text{cap}}(q) = I_{\text{glass}}(q) + T_{\text{cap}} I_{\text{air}}(q)$. Therefore, subtracting the capillary background inherently removes the empty beamline background. Subtracting the empty beamline background a second time would cause unphysical negative intensities at intermediate $q$ (double-subtraction error).


### 4. Propagated uncertainty analysis
#### 1D analytical variance propagation (Spiger `integrator20.py`):
$$\sigma_{\text{sub}}^2(q) = K_{\text{std}}^2 \left[ \sigma_{\text{sam,norm}}^2(q) + \sigma_{\text{bkg,norm}}^2(q) \right] + \sigma_{K_{\text{std}}}^2 \left[ I_{\text{sam,norm}}(q) - I_{\text{bkg,norm}}(q) \right]^2$$
$$\sigma_{\text{sam,norm}}^2(q) = \underbrace{\left( \frac{\sigma_{\text{sam}}(q)}{I_{0, \text{sam}} T_{\text{rel}}} \right)^2}_{\text{Sample Poisson counting}} + \underbrace{\left( \frac{I_{\text{sam}}(q) \sigma_{T_{\text{rel}}}}{I_{0, \text{sam}} T_{\text{rel}}^2} \right)^2}_{\text{Relative transmission error}} + \underbrace{\left( \frac{I_{\text{sam}}(q) \sigma_{I_{0, \text{sam}}}}{T_{\text{rel}} I_{0, \text{sam}}^2} \right)^2}_{\text{Sample monitor fluctuation}}$$
$$\sigma_{\text{bkg,norm}}^2(q) = \underbrace{\left( \frac{\sigma_{\text{bkg}}(q)}{I_{0, \text{bkg}}} \right)^2}_{\text{Background Poisson counting}} + \underbrace{\left( \frac{I_{\text{bkg}}(q) \sigma_{I_{0, \text{bkg}}}}{I_{0, \text{bkg}}^2} \right)^2}_{\text{Background monitor fluctuation}}$$

#### 2D spatial uncertainty and SNR mapping (diagnostic):
The pixel variance $\sigma^2(x, y)$ propagates seven independent components across the detector face for spatial diagnostic inspection and pixel signal-to-noise mapping ($\text{SNR}(x, y) = |I_{\text{sub}}(x, y)| / \sigma(x, y)$).


### 5. Azimuthal Integration, Sector Slicing & Cake Transformation
Using `pyFAI` CSR Cython algorithms with solid-angle and polarization corrections:
- **1D Profile**: Azimuthal integration over $0^\circ \le \phi \le 360^\circ$ into $q = \frac{4\pi}{\lambda}\sin\theta$ ($2,150$ radial bins from $0.089\ \mathrm{nm}^{-1}$ to $42.08\ \mathrm{nm}^{-1}$).
- **Azimuthal Sectors**: Meridian ($45^\circ \le \phi \le 135^\circ$) and Equatorial ($-45^\circ \le \phi \le 45^\circ$), integrated independently for sample and background before subtraction.
- **2D Caking**: Transformation to $I(q, \chi)$ ($400$ radial $\times 360$ azimuthal bins), integrated independently and subtracted.

---

## Directory Structure & Data Dictionary

```
saxs_publication_project/
├── README.md                                 # This scientific and technical manual
├── METADATA.json                             # Machine-readable schema (Zenodo / OSF compliant)
├── requirements.txt                          # Python dependencies with pinned versions
├── saxs_reduction_pipeline_demo.py           # Standalone executable reduction pipeline
├── saxs_data_reduction_report.html           # Interactive, publication-ready HTML report
│
├── metadata/                                 # Experimental logs, motor positions, and macros
│   ├── spec_files/
│   │   └── ML_18keV_17th_March_2026.dat      # Master SPEC file (scans #1 to #300+, counters, motor positions)
│   └── macros/
│       ├── Tomasz_March_2026_FlowTop.mac     # SPEC acquisition macro for the in-situ flow-through series
│       ├── Tomasz_March_2026_FlowTop_restart.mac # Restart macro for continuous acquisition
│       └── Tomasz_March_2026_Quartz.mac      # Alignment and quartz capillary calibration macro
│
├── calibration/                              # PyFAI detector geometry, Spiger project state, and masks
│   ├── calib.poni                            # pyFAI PONI file (18 keV, 37.33 cm distance, beam center)
│   ├── Start.spig                            # Spiger Python pickle project file (contains Glassy Carbon factor)
│   ├── Start_ide.tiff                        # Composite 2D detector mask (3110 x 3269 pixels)
│   ├── Start_rigorous_mask.tiff              # Expanded detector mask (pane dilation, dead/hot pixels)
│   └── Start_calibration.hdf5                # Calibration reference curves and integration bins
│
├── raw_data/                                 # Representative raw Eiger 9M detector files (HDF5 format)
│   ├── sample_scan184/                       # Flow-through sample scan (Scan #184, 55 frames, ~517 MB)
│   │   ├── LUKU2026030601PF_000184_data_000001.h5
│   │   └── ... (frames 000001 to 000056)
│   └── background_scan063/                   # Capillary buffer background scan (Scan #63, 11 frames, ~98 MB)
│       ├── CapillaryTop_000063_data_000001.h5
│       └── ... (frames 000001 to 000012)
│
├── reduced_data/                             # Reduced 1D and 2D profiles with uncertainties
│   ├── Exp2_1/                               # In-situ FlowThrough Part 1 (Scans #184 to #274)
│   │   ├── Diffraction/HDF5/
│   │   │   ├── MultiRed_All/MultiRed_All.hdf5 # Full azimuthal integration (0°-360°)
│   │   │   ├── MultiRed_Equatorial/          # Equatorial sector (-45° to +45°)
│   │   │   └── MultiRed_Meridian/            # Meridian sector (45° to 135°)
│   │   └── Transmission/Txt/Transmission.txt # Transmission table and ion chamber counter logs
│   └── Exp2_2/                               # In-situ FlowThrough Part 2 (Scans #275 to #292)
│       ├── Diffraction/HDF5/                 # MultiRed_All, MultiRed_Equatorial, MultiRed_Meridian
│       └── Transmission/Txt/Transmission.txt
│
└── publication_figures/                      # Figures in PDF and PNG format
    ├── Figure_SAXS_Reduction_Pipeline.pdf/.png # Multi-panel overview figure (panels a-f)
    ├── step1_raw_2d_detector.pdf/.png        # Raw 2D photon counts for sample and buffer
    ├── step2_mask_application.pdf/.png       # Mask and active pixel area
    ├── step3_monitor_normalization.pdf/.png  # Flux-normalized 2D patterns
    ├── step4_background_subtraction_2d.pdf/.png # Transmission-corrected net subtracted pattern
    ├── step4b_uncertainty_and_snr_2d.pdf/.png # 2D analytical uncertainty map and pixel SNR
    ├── step5_azimuthal_cake_plot.pdf/.png    # 2D polar cake transformation I(q, chi)
    └── step6_1d_integrated_profiles.pdf/.png # Calibrated 1D profile with +/-1sigma uncertainty band
```

### HDF5 Reduced Data Schema (`MultiRed_All.hdf5`)
Each measurement frame is stored as an individual HDF5 group `Iq_XXXXXX`:
- **`profile` dataset**: Shape `(2150, 4)` float64 array:
  - Column 0: Scattering vector $q\ [\mathrm{nm}^{-1}]$
  - Column 1: Calibrated absolute differential scattering cross-section $\frac{d\Sigma}{d\Omega}(q)\ [\mathrm{cm}^{-1}]$
  - Column 2: Propagated $1\sigma_I(q)$ uncertainty $[\mathrm{cm}^{-1}]$
  - Column 3: Instrumental resolution / $q$-smearing $\Delta q\ [\mathrm{nm}^{-1}]$
- **Group Attributes**:
  - `SPEC Scan - Data`: Integer scan number in SPEC log
  - `Epoch - Data`: Unix timestamp of the acquisition
  - `Transmission - Data`: Measured transmission factor and uncertainty
  - `Monitor - Data`: Incident beam ion chamber counts $I_0$
  - `Seconds - Data`: Exposure time in seconds

---

## Reproducibility & Execution Guide

### 1. Environment Setup
Create and activate a Python virtual environment, then install dependencies:
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### 2. Single-Scan Data Reduction & Figure Generation
Run the reduction pipeline on the primary frame (Scan #184 vs Background #63):
```bash
python3 saxs_reduction_pipeline_demo.py --scan 184 --frame 1 --bkg 63
```
*Generates Step 1 to Step 6 figures, 2D uncertainty and SNR maps, and the multi-panel overview figure in `publication_figures/`.*

### 3. Execute Full Pipeline Suite
```bash
python3 saxs_reduction_pipeline_demo.py --all
```

---

## Interactive Report

The report is provided in [`saxs_data_reduction_report.html`](file:///run/media/tomek/data/Mrc_26/saxs_publication_project/saxs_data_reduction_report.html). It includes:
- Mathematical derivations and processing steps.
- Figure modal viewers with zoom.
- Uncertainty budget breakdown table.
- Dark/Light mode toggle.
- Links to data files and reproduction commands.

To view the report, open `saxs_data_reduction_report.html` in a web browser:
```bash
xdg-open saxs_data_reduction_report.html
```

---

## Data Governance & Provenance

- **Storage**: All files within this directory are physical data copies without symlinks (`find . -type l` returns `0`).
- **Archive Format**: Compliant with scientific open-access repositories (Zenodo, Materials Commons, OSF). See [`METADATA.json`](file:///run/media/tomek/data/Mrc_26/saxs_publication_project/METADATA.json) for machine-readable JSON-LD metadata.
- **License**: Creative Commons Attribution 4.0 International (CC BY 4.0) for data and documentation; MIT License for analysis scripts.
