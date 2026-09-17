#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SAXS Data Reduction Pipeline Demonstration & Publication Figures
================================================================
Focus: Experiment 2_1, Measurements Scan #184 to #292
Instrument: BESSY II µSpot (myspot) beamline (18 keV, Eiger 9M detector)
Software baseline: Spiger 2.0 (Start.spig)

Includes Comprehensive Uncertainty Analysis & Error Propagation:
  1. Pixel-level 7-component analytical error propagation (integrator20.py)
  2. 2D spatial uncertainty and signal-to-noise ratio (SNR) mapping
  3. 1D radial profile uncertainty bands, error bars, and fractional error delta(I)/I
  4. In-situ 16.15-hour kinetic evolution across Scan #184 to #292

Usage:
  python3 saxs_reduction_pipeline_demo.py --scan 184 --frame 1
  python3 saxs_reduction_pipeline_demo.py --series 184 292
  python3 saxs_reduction_pipeline_demo.py --all

Author: Antigravity Pair-Programming for Tomasz
Date: March 2026 Beamtime / September 2026 Analysis
"""

import os
import sys
import argparse
import pickle
import numpy as np
import h5py
import hdf5plugin  # Required for Eiger bitshuffle/lz4 decompression
import fabio
import pyFAI
from pyFAI.integrator.azimuthal import AzimuthalIntegrator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm, Normalize
import matplotlib.gridspec as gridspec

# -----------------------------------------------------------------------------
# Configuration & Paths (Dynamic & Self-Contained)
# -----------------------------------------------------------------------------
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PUB_PROJECT_DIR = '/run/media/tomek/data/Mrc_26/saxs_publication_project'

if os.path.exists(os.path.join(CURRENT_DIR, 'calibration', 'calib.poni')):
    PROJECT_DIR = CURRENT_DIR
elif os.path.exists(os.path.join(PUB_PROJECT_DIR, 'calibration', 'calib.poni')):
    PROJECT_DIR = PUB_PROJECT_DIR
else:
    PROJECT_DIR = CURRENT_DIR

# Primary local paths
CALIB_DIR = os.path.join(PROJECT_DIR, 'calibration')
PONI_FILE = os.path.join(CALIB_DIR, 'calib.poni')
SPIG_FILE = os.path.join(CALIB_DIR, 'Start.spig')
MASK_FILE = os.path.join(CALIB_DIR, 'Start_ide.tiff')
RIGOROUS_MASK_FILE = os.path.join(CALIB_DIR, 'Start_rigorous_mask.tiff')

SPEC_FILE = os.path.join(PROJECT_DIR, 'metadata', 'spec_files', 'ML_18keV_17th_March_2026.dat')

REDUCED_EXP2_1 = os.path.join(PROJECT_DIR, 'reduced_data', 'Exp2_1', 
                              'Diffraction', 'HDF5', 'MultiRed_All', 'MultiRed_All.hdf5')
REDUCED_EXP2_2 = os.path.join(PROJECT_DIR, 'reduced_data', 'Exp2_2', 
                              'Diffraction', 'HDF5', 'MultiRed_All', 'MultiRed_All.hdf5')

OUTPUT_DIR = os.path.join(CURRENT_DIR, 'publication_figures')
os.makedirs(OUTPUT_DIR, exist_ok=True)

def find_raw_eiger_file(scan_num, frame_num, prefix=None):
    """Locates the raw Eiger .h5 file in local raw_data/ folders or original muspot_data."""
    import glob
    search_dirs = [
        os.path.join(PROJECT_DIR, 'raw_data', f'sample_scan{scan_num:03d}'),
        os.path.join(PROJECT_DIR, 'raw_data', f'background_scan{scan_num:03d}'),
        os.path.join(PROJECT_DIR, 'raw_data'),
        os.path.join(PROJECT_DIR, 'raw_data', f'sample_scan{scan_num}'),
        os.path.join(PROJECT_DIR, 'raw_data', f'background_scan{scan_num}'),
        '/run/media/tomek/data/Mrc_26/muspot_data/eiger_files'
    ]
    patterns = [
        f"*_{scan_num:06d}_data_{frame_num:06d}.h5",
        f"*{scan_num:06d}*{frame_num:06d}*.h5",
        f"*{scan_num}*data*{frame_num:06d}*.h5"
    ]
    if prefix:
        patterns.insert(0, f"{prefix}_{scan_num:06d}_data_{frame_num:06d}.h5")
        
    for sdir in search_dirs:
        if not os.path.exists(sdir):
            continue
        for pat in patterns:
            matches = glob.glob(os.path.join(sdir, pat))
            if matches:
                return matches[0]
            matches = glob.glob(os.path.join(sdir, '**', pat), recursive=True)
            if matches:
                return matches[0]
    raise FileNotFoundError(f"Raw Eiger file not found for scan #{scan_num}, frame #{frame_num}")


# Publication styling
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 14,
    'mathtext.fontset': 'dejavusans',
    'figure.autolayout': False
})


def load_spiger_state():
    """Loads calibration geometry, mask, and glassy carbon factor from Spiger Start.spig."""
    with open(SPIG_FILE, 'rb') as f:
        spig = pickle.load(f)

    ai = pyFAI.load(PONI_FILE)
    mask_ide = fabio.open(MASK_FILE).data
    base_mask = (mask_ide != 0)
    K_std = float(spig['glassy']['standard factor'])
    K_std_err = float(spig['glassy']['standard error'])

    return spig, ai, base_mask, K_std, K_std_err


def build_rigorous_mask(img_sample_raw, img_bkg_raw, base_mask, ai, save_path=None):
    """
    Constructs an expanded detector mask (extending the baseline mask) that covers:
      1. Beamstop arm and shadows from base mask
      2. Detector sentinel values (>= 4,294,967,000)
      3. Inter-module gaps with 2-pixel dilation
      4. Outer detector perimeter guard margin (2 pixels)
      5. Dead and cold pixels (counts <= 0 or counts < 25% of local 3x3 median)
      6. Hot pixels (> 7 sigma above local 3x3 median)
    """
    from scipy.ndimage import binary_dilation, median_filter

    rigorous_mask = base_mask.copy().astype(bool)

    # 1. Detector sentinels
    sentinels = (img_sample_raw >= 4294967000) | (img_bkg_raw >= 4294967000)
    rigorous_mask |= sentinels

    # 2. Module pane gaps + 2-pixel perimeter dilation (covers double pixels and border crosstalk)
    det = pyFAI.detector_factory('Eiger9M')
    pyfai_gaps = det.mask != 0
    dilated_gaps = binary_dilation(pyfai_gaps, iterations=2)
    rigorous_mask |= dilated_gaps

    # 3. Outer detector perimeter margin (2 pixels)
    rigorous_mask[0:2, :] = True
    rigorous_mask[-2:, :] = True
    rigorous_mask[:, 0:2] = True
    rigorous_mask[:, -2:] = True

    # 4. Dead & cold pixels
    s_clean = np.where(sentinels, 0, img_sample_raw)
    b_clean = np.where(sentinels, 0, img_bkg_raw)
    med_s = median_filter(s_clean, size=3)
    med_b = median_filter(b_clean, size=3)

    zero_dead = (img_sample_raw <= 0) | (img_bkg_raw <= 0)
    cold_dead = ((img_sample_raw < 0.25 * med_s) & (med_s > 50)) | ((img_bkg_raw < 0.25 * med_b) & (med_b > 40))
    dead_pixels = (zero_dead | cold_dead) & (~sentinels)
    rigorous_mask |= dead_pixels

    # 5. Hot pixels / isolated spikes (> 7 sigma above local median outside beam center)
    diff_s = s_clean - med_s
    diff_b = b_clean - med_b
    hot_s = (diff_s > np.maximum(7.0 * np.sqrt(np.maximum(med_s, 1.0)), 50.0)) & (s_clean > 1.8 * med_s)
    hot_b = (diff_b > np.maximum(7.0 * np.sqrt(np.maximum(med_b, 1.0)), 50.0)) & (b_clean > 1.8 * med_b)
    hot_pixels = (hot_s | hot_b) & (~sentinels)
    rigorous_mask |= hot_pixels

    print(f"  Detector Mask Generated: {np.sum(rigorous_mask)} pixels masked ({np.sum(rigorous_mask)/rigorous_mask.size*100:.2f}%)")
    print(f"    - Pane gaps & borders: {np.sum(dilated_gaps)} px")
    print(f"    - Dead/cold pixels: {np.sum(dead_pixels)} px")
    print(f"    - Hot pixels/spikes: {np.sum(hot_pixels)} px")

    if save_path:
        import fabio
        fabio.tifimage.TifImage(data=rigorous_mask.astype(np.uint8)).write(save_path)
        print(f"    - Saved mask to: {os.path.basename(save_path)}")

    return rigorous_mask


def compute_2d_uncertainty(d_raw, b_raw, i0_d, i0_b, t_rel, t_rel_err, gc, gc_err, mon_perc=0.02):
    """
    Computes exact pixel-by-pixel analytical variance propagation matching Spiger integrator20.py:
    sigma^2 = Var(P_sam) + Var(T_rel) + Var(I0_sam) + Var(Gc_sam) + Var(P_bkg) + Var(I0_bkg) + Var(Gc_bkg)
    """
    i0_err_d = i0_d * mon_perc
    i0_err_b = i0_b * mon_perc

    d_safe = np.maximum(d_raw, 0)
    b_safe = np.maximum(b_raw, 0)

    # 7 analytical terms from integrator20.py lines 537-538:
    term1 = (gc * np.sqrt(d_safe) / (i0_d * t_rel)) ** 2                # Sample Poisson noise
    term2 = (gc * d_raw / i0_d * t_rel_err / (t_rel ** 2)) ** 2        # Relative Transmission uncertainty
    term3 = (gc * d_raw / t_rel * i0_err_d / (i0_d ** 2)) ** 2        # Sample Monitor uncertainty
    term4 = (gc_err * d_raw / (i0_d * t_rel)) ** 2                    # Standard Factor (Sample)
    term5 = (gc * np.sqrt(b_safe) / i0_b) ** 2                         # Background Poisson noise
    term6 = (gc * b_raw * i0_err_b / (i0_b ** 2)) ** 2                 # Background Monitor uncertainty
    term7 = (gc_err * b_raw / i0_b) ** 2                              # Standard Factor (Background)

    sigma_2d = np.sqrt(term1 + term2 + term3 + term4 + term5 + term6 + term7)
    
    terms = {
        'poisson_sample': term1,
        'trans_rel': term2,
        'monitor_sample': term3,
        'std_sample': term4,
        'poisson_bkg': term5,
        'monitor_bkg': term6,
        'std_bkg': term7
    }
    return sigma_2d, terms


def run_single_scan_reduction(scan_num=184, frame_num=1, bkg_scan=63, bkg_frame=1):
    """
    Executes the step-by-step data reduction on a single raw frame using Spiger math,
    and produces all individual step figures and the master multi-panel figure.
    """
    print(f"\n--- Running Single-Scan Reduction for Scan #{scan_num} (frame {frame_num}) vs Bkg #{bkg_scan} ---")
    spig, ai, base_mask, K_std, K_std_err = load_spiger_state()
    spec_data = spig['spec']['data']
    motors = spec_data['motors']

    # Locate SPEC rows
    row_sample = None
    for r in spec_data['data']:
        if r[motors['run']] == str(scan_num) and r[motors['scan']] == (frame_num - 1):
            row_sample = r
            break
    if row_sample is None:
        for r in spec_data['data']:
            if r[motors['run']] == str(scan_num):
                row_sample = r
                break
    if row_sample is None:
        raise ValueError(f"Could not find scan #{scan_num} in SPEC data!")

    row_bkg = None
    for r in spec_data['data']:
        if r[motors['run']] == str(bkg_scan) and r[motors['scan']] == (bkg_frame - 1):
            row_bkg = r
            break
    if row_bkg is None:
        for r in spec_data['data']:
            if r[motors['run']] == str(bkg_scan):
                row_bkg = r
                break
    if row_bkg is None:
        raise ValueError(f"Could not find background scan #{bkg_scan} in SPEC data!")

    # Extract sample metadata
    sec_d = float(row_sample[motors['Seconds']])
    i0_d = float(row_sample[motors['IoniCh']])
    roi3_d = float(row_sample[motors['roi3']])
    ltime_d = float(row_sample[motors['ltime']])
    trans_d = (roi3_d / ltime_d * sec_d) / i0_d
    trans_err_d = trans_d * np.sqrt(1.0 / roi3_d + 0.02**2)

    # Extract background metadata
    sec_b = float(row_bkg[motors['Seconds']])
    i0_b = float(row_bkg[motors['IoniCh']])
    roi3_b = float(row_bkg[motors['roi3']])
    ltime_b = float(row_bkg[motors['ltime']])
    trans_b = (roi3_b / ltime_b * sec_b) / i0_b
    trans_err_b = trans_b * np.sqrt(1.0 / roi3_b + 0.02**2)

    # Relative transmission & propagated uncertainty
    t_rel = trans_d / trans_b
    t_rel_err = t_rel * np.sqrt((trans_err_d / trans_d)**2 + (trans_err_b / trans_b)**2)

    print(f"  Sample #{scan_num}: Sec={sec_d:.1f}s, IoniCh={i0_d:.0f}, Trans={trans_d:.5f} ± {trans_err_d:.5f}")
    print(f"  Bkg #{bkg_scan}: Sec={sec_b:.1f}s, IoniCh={i0_b:.0f}, Trans={trans_b:.5f} ± {trans_err_b:.5f}")
    print(f"  Relative Transmission T_rel = {t_rel:.5f} ± {t_rel_err:.5f}")

    # File paths
    f_sample = find_raw_eiger_file(scan_num, frame_num, prefix="LUKU2026030601PF")
    f_bkg = find_raw_eiger_file(bkg_scan, bkg_frame, prefix="CapillaryTop")

    print(f"  Loading Sample raw image: {os.path.basename(f_sample)}")
    img_sample_raw = fabio.open(f_sample).data.astype(np.float64)
    print(f"  Loading Bkg raw image: {os.path.basename(f_bkg)}")
    img_bkg_raw = fabio.open(f_bkg).data.astype(np.float64)

    # Build detector mask (dilated pane gaps, dead pixels, hot pixels)
    rigorous_mask_path = os.path.join(CALIB_DIR, 'Start_rigorous_mask.tiff')
    mask_full = build_rigorous_mask(
        img_sample_raw, img_bkg_raw, base_mask, ai,
        save_path=rigorous_mask_path
    )

    # Spiger mathematical reduction:
    img_sample_norm = img_sample_raw / i0_d
    img_bkg_norm = img_bkg_raw / i0_b

    # Subtraction with transmission scaling & absolute standard factor:
    # I_sub = K_std * [ I_sample / (i0_d * t_rel) - I_bkg / i0_b ]
    img_sub_2d = K_std * (img_sample_raw / (i0_d * t_rel) - img_bkg_raw / i0_b)
    img_sub_2d_masked = img_sub_2d.copy()
    img_sub_2d_masked[mask_full] = np.nan

    # Compute 2D analytical uncertainty map
    print("  Calculating analytical 2D error propagation map...")
    sigma_2d, err_terms = compute_2d_uncertainty(
        img_sample_raw, img_bkg_raw, i0_d, i0_b, t_rel, t_rel_err, K_std, K_std_err
    )
    sigma_2d_masked = sigma_2d.copy()
    sigma_2d_masked[mask_full] = np.nan

    # Signal-to-Noise Ratio (SNR) map: |I_sub| / sigma
    snr_2d = np.abs(img_sub_2d) / np.maximum(sigma_2d, 1e-9)
    snr_2d_masked = snr_2d.copy()
    snr_2d_masked[mask_full] = np.nan

    # Azimuthal Radial Integrations (with exact Poisson error propagation)
    print("  Integrating 1D profiles and 2D Cake...")
    int_all = ai.integrate1d(img_sub_2d, npt=2150, correctSolidAngle=True, unit='q_nm^-1',
                             mask=mask_full, error_model="poisson", method=('full', 'CSR', 'cython'))
    q_all, i_all = int_all.radial, int_all.intensity

    # Benchmark comparison against Spiger saved profile (Iq_000166)
    with h5py.File(REDUCED_EXP2_1, 'r') as h5_ref:
        ref_profile = h5_ref['Iq_000166/profile'][:]
        q_ref = ref_profile[:, 0]
        i_ref = ref_profile[:, 1]
        sigma_ref = ref_profile[:, 2]     # Propagated 1D uncertainty
        dq_ref = ref_profile[:, 3]        # Resolution / q-smearing

    int_meridian = ai.integrate1d(img_sub_2d, npt=2150, correctSolidAngle=True, unit='q_nm^-1',
                                  mask=mask_full, azimuth_range=(45, 135), method=('full', 'CSR', 'cython'))
    q_meridian, i_meridian = int_meridian.radial, int_meridian.intensity

    int_equatorial = ai.integrate1d(img_sub_2d, npt=2150, correctSolidAngle=True, unit='q_nm^-1',
                                    mask=mask_full, azimuth_range=(-45, 45), method=('full', 'CSR', 'cython'))
    q_equatorial, i_equatorial = int_equatorial.radial, int_equatorial.intensity

    # 2D Cake transformation
    cake_res = ai.integrate2d(img_sub_2d, npt_rad=400, npt_azim=360, correctSolidAngle=True,
                              unit='q_nm^-1', mask=mask_full, method=('full', 'CSR', 'cython'))
    q_cake, chi_cake, img_cake = cake_res.radial, cake_res.azimuthal, cake_res.intensity

    # Direct beam center
    xcen = ai.poni2 / ai.pixel2
    ycen = ai.poni1 / ai.pixel1

    # -------------------------------------------------------------------------
    # Generate Figures
    # -------------------------------------------------------------------------
    print("  Plotting individual figures...")

    # Step 1: Raw 2D Frame
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    axes[0].set_facecolor('#111827')
    axes[1].set_facecolor('#111827')
    cmap_inf0 = plt.cm.inferno.copy()
    cmap_inf0.set_bad('#111827')
    raw_disp_sample = np.where(img_sample_raw >= 4294967000, np.nan, img_sample_raw)
    im0 = axes[0].imshow(raw_disp_sample, origin='lower', cmap=cmap_inf0, norm=LogNorm(vmin=1, vmax=5e4))
    axes[0].plot(xcen, ycen, 'r+', markersize=14, markeredgewidth=2, label='Beam Center')
    axes[0].set_title(rf"$\mathbf{{(a)}}$ Sample Raw 2D Frame (#{scan_num})", pad=10)
    axes[0].set_xlabel("Detector Pixel X")
    axes[0].set_ylabel("Detector Pixel Y")
    axes[0].legend(loc='upper right', framealpha=0.85)
    cb0 = fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)
    cb0.set_label(r"Raw Photon Counts [$\mathrm{log}_{10}$]")

    cmap_inf1 = plt.cm.inferno.copy()
    cmap_inf1.set_bad('#111827')
    raw_disp_bkg = np.where(img_bkg_raw >= 4294967000, np.nan, img_bkg_raw)
    im1 = axes[1].imshow(raw_disp_bkg, origin='lower', cmap=cmap_inf1, norm=LogNorm(vmin=1, vmax=5e4))
    axes[1].plot(xcen, ycen, 'r+', markersize=14, markeredgewidth=2, label='Beam Center')
    axes[1].set_title(rf"$\mathbf{{(b)}}$ Background Raw 2D Frame (#{bkg_scan})", pad=10)
    axes[1].set_xlabel("Detector Pixel X")
    axes[1].set_ylabel("Detector Pixel Y")
    axes[1].legend(loc='upper right', framealpha=0.85)
    cb1 = fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)
    cb1.set_label(r"Raw Photon Counts [$\mathrm{log}_{10}$]")
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "step1_raw_2d_detector.png"), dpi=300)
    fig.savefig(os.path.join(OUTPUT_DIR, "step1_raw_2d_detector.pdf"))
    plt.close(fig)

    # Step 2: Mask Application
    from matplotlib.colors import ListedColormap
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    axes[0].set_facecolor('#111827')
    mask_binary_cmap = ListedColormap(['#111827', '#00e5ff'])
    axes[0].imshow(mask_full, origin='lower', cmap=mask_binary_cmap, interpolation='nearest')
    axes[0].plot(xcen, ycen, 'r+', markersize=14, markeredgewidth=2, label='Beam Center')
    axes[0].set_title(r"$\mathbf{(a)}$ Detector Mask Matrix (Cyan = Masked)", pad=10)
    axes[0].set_xlabel("Detector Pixel X")
    axes[0].set_ylabel("Detector Pixel Y")
    axes[0].legend(loc='upper right', framealpha=0.85)

    axes[1].set_facecolor('#111827')
    cmap_inf_s2 = plt.cm.inferno.copy()
    cmap_inf_s2.set_bad('#111827')
    im_mask_ov = axes[1].imshow(raw_disp_sample, origin='lower', cmap=cmap_inf_s2, norm=LogNorm(vmin=1, vmax=5e4))

    mask_cyan_cmap = ListedColormap(['#00e5ff'])  # Pure high-contrast cyan
    mask_cyan_cmap.set_bad(alpha=0.0)
    mask_overlay = np.where(mask_full, 1.0, np.nan)
    axes[1].imshow(mask_overlay, origin='lower', cmap=mask_cyan_cmap, alpha=1.0, interpolation='nearest')
    axes[1].plot(xcen, ycen, 'r+', markersize=14, markeredgewidth=2, label='Beam Center')
    axes[1].set_title(r"$\mathbf{(b)}$ Active Data Area (Mask Shaded in Cyan)", pad=10)
    axes[1].set_xlabel("Detector Pixel X")
    axes[1].set_ylabel("Detector Pixel Y")
    axes[1].legend(loc='upper right', framealpha=0.85)
    cb_m = fig.colorbar(im_mask_ov, ax=axes[1], fraction=0.046, pad=0.04)
    cb_m.set_label(r"Photon Counts [$\mathrm{log}_{10}$]")
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "step2_mask_application.png"), dpi=300)
    fig.savefig(os.path.join(OUTPUT_DIR, "step2_mask_application.pdf"))
    plt.close(fig)

    # Step 3: Monitor Normalization
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    axes[0].set_facecolor('#111827')
    axes[1].set_facecolor('#111827')
    cmap_inf3_0 = plt.cm.inferno.copy()
    cmap_inf3_0.set_bad('#111827')
    norm_sample_disp = np.ma.masked_where(mask_full, img_sample_norm)
    im_ns = axes[0].imshow(norm_sample_disp, origin='lower', cmap=cmap_inf3_0, norm=LogNorm(vmin=1e-6, vmax=1e-1))
    axes[0].set_title(rf"$\mathbf{{(a)}}$ Sample Normalized by Monitor ($I_0 = {i0_d:.0f}$)", pad=10)
    axes[0].set_xlabel("Detector Pixel X")
    axes[0].set_ylabel("Detector Pixel Y")
    cb_ns = fig.colorbar(im_ns, ax=axes[0], fraction=0.046, pad=0.04)
    cb_ns.set_label(r"Normalized Intensity $I / I_0$ [counts/monitor]")

    cmap_inf3_1 = plt.cm.inferno.copy()
    cmap_inf3_1.set_bad('#111827')
    norm_bkg_disp = np.ma.masked_where(mask_full, img_bkg_norm)
    im_nb = axes[1].imshow(norm_bkg_disp, origin='lower', cmap=cmap_inf3_1, norm=LogNorm(vmin=1e-6, vmax=1e-1))
    axes[1].set_title(rf"$\mathbf{{(b)}}$ Background Normalized by Monitor ($I_{{0,bkg}} = {i0_b:.0f}$)", pad=10)
    axes[1].set_xlabel("Detector Pixel X")
    axes[1].set_ylabel("Detector Pixel Y")
    cb_nb = fig.colorbar(im_nb, ax=axes[1], fraction=0.046, pad=0.04)
    cb_nb.set_label(r"Normalized Intensity $I / I_0$ [counts/monitor]")
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "step3_monitor_normalization.png"), dpi=300)
    fig.savefig(os.path.join(OUTPUT_DIR, "step3_monitor_normalization.pdf"))
    plt.close(fig)

    # Step 4: Transmission & Subtraction
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    axes[0].set_facecolor('#111827')
    axes[1].set_facecolor('#111827')
    axes[2].set_facecolor('#111827')
    cmap_inf4_0 = plt.cm.inferno.copy()
    cmap_inf4_0.set_bad('#111827')
    sample_trans_disp = np.ma.masked_where(mask_full, img_sample_raw / (i0_d * t_rel))
    im_st = axes[0].imshow(sample_trans_disp, origin='lower', cmap=cmap_inf4_0, norm=LogNorm(vmin=1e-6, vmax=1e-1))
    axes[0].set_title(rf"$\mathbf{{(a)}}$ Sample Trans.-Corrected ($T_{{rel}} = {t_rel:.4f}$)", pad=10)
    axes[0].set_xlabel("Detector Pixel X")
    axes[0].set_ylabel("Detector Pixel Y")
    cb_st = fig.colorbar(im_st, ax=axes[0], fraction=0.046, pad=0.04)
    cb_st.set_label(r"$I_{sample} / (I_0 \cdot T_{rel})$")

    im_bt = axes[1].imshow(norm_bkg_disp, origin='lower', cmap=cmap_inf4_0, norm=LogNorm(vmin=1e-6, vmax=1e-1))
    axes[1].set_title(r"$\mathbf{(b)}$ Background Buffer ($I_{bkg} / I_{0,bkg}$)", pad=10)
    axes[1].set_xlabel("Detector Pixel X")
    axes[1].set_ylabel("Detector Pixel Y")
    cb_bt = fig.colorbar(im_bt, ax=axes[1], fraction=0.046, pad=0.04)
    cb_bt.set_label(r"$I_{bkg} / I_{0,bkg}$")

    cmap_cw = plt.cm.coolwarm.copy()
    cmap_cw.set_bad('#111827')
    im_sub = axes[2].imshow(img_sub_2d_masked, origin='lower', cmap=cmap_cw, 
                           norm=SymLogNorm(linthresh=0.01, linscale=1.0, vmin=-0.2, vmax=5.0))
    axes[2].set_title(r"$\mathbf{(c)}$ Net Sample 2D Pattern ($I_{sub}$)", pad=10)
    axes[2].set_xlabel("Detector Pixel X")
    axes[2].set_ylabel("Detector Pixel Y")
    cb_sub = fig.colorbar(im_sub, ax=axes[2], fraction=0.046, pad=0.04)
    cb_sub.set_label(r"Absolute Cross-Section $[\mathrm{cm}^{-1}]$")
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "step4_background_subtraction_2d.png"), dpi=300)
    fig.savefig(os.path.join(OUTPUT_DIR, "step4_background_subtraction_2d.pdf"))
    plt.close(fig)

    # Step 4b: Comprehensive 2D Spatial Uncertainty & SNR Mapping
    print("  Generating Fig 4b: 2D Spatial Uncertainty & SNR Map...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    axes[0].set_facecolor('#111827')
    axes[1].set_facecolor('#111827')
    cmap_inf4b_0 = plt.cm.inferno.copy()
    cmap_inf4b_0.set_bad('#111827')
    im_sig = axes[0].imshow(sigma_2d_masked, origin='lower', cmap=cmap_inf4b_0, norm=LogNorm(vmin=5e-3, vmax=0.5))
    axes[0].plot(xcen, ycen, 'r+', markersize=14, markeredgewidth=2)
    axes[0].set_title(r"$\mathbf{(a)}$ Propagated 2D Uncertainty $\sigma(x, y)\ [\mathrm{cm}^{-1}]$", pad=10)
    axes[0].set_xlabel("Detector Pixel X")
    axes[0].set_ylabel("Detector Pixel Y")
    cb_sig = fig.colorbar(im_sig, ax=axes[0], fraction=0.046, pad=0.04)
    cb_sig.set_label(r"Pixel Uncertainty $\sigma(x, y)\ [\mathrm{cm}^{-1}]$")

    cmap_inf4b_1 = plt.cm.inferno.copy()
    cmap_inf4b_1.set_bad('#111827')
    im_snr = axes[1].imshow(snr_2d_masked, origin='lower', cmap=cmap_inf4b_1, norm=LogNorm(vmin=0.1, vmax=50))
    axes[1].plot(xcen, ycen, 'r+', markersize=14, markeredgewidth=2)
    axes[1].set_title(r"$\mathbf{(b)}$ Pixel Signal-to-Noise Ratio $\mathrm{SNR}(x, y) = |I_{sub}| / \sigma$", pad=10)
    axes[1].set_xlabel("Detector Pixel X")
    axes[1].set_ylabel("Detector Pixel Y")
    cb_snr = fig.colorbar(im_snr, ax=axes[1], fraction=0.046, pad=0.04)
    cb_snr.set_label(r"$\mathrm{SNR}(x, y)$")
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "step4b_uncertainty_and_snr_2d.png"), dpi=300)
    fig.savefig(os.path.join(OUTPUT_DIR, "step4b_uncertainty_and_snr_2d.pdf"))
    plt.close(fig)

    # Step 5: Azimuthal Cake Transformation
    fig, ax = plt.subplots(figsize=(10, 6))
    cake_disp = np.ma.masked_invalid(img_cake)
    im_cake = ax.imshow(cake_disp, origin='lower', aspect='auto', cmap='inferno',
                        extent=[q_cake.min(), q_cake.max(), chi_cake.min(), chi_cake.max()],
                        norm=LogNorm(vmin=1e-2, vmax=10))
    ax.set_title(r"$\mathbf{2D\ Azimuthal\ Cake\ Transformation}\ I(q, \chi)$", pad=12)
    ax.set_xlabel(r"Scattering Vector $q\ [\mathrm{nm}^{-1}]$")
    ax.set_ylabel(r"Azimuthal Angle $\chi\ [^\circ]$")
    ax.axhline(45, color='cyan', linestyle='--', linewidth=1.5, label=r'Meridian Sector ($45^\circ - 135^\circ$)')
    ax.axhline(135, color='cyan', linestyle='--', linewidth=1.5)
    ax.axhline(-45, color='lime', linestyle=':', linewidth=1.5, label=r'Equatorial Sector ($-45^\circ - 45^\circ$)')
    ax.axhline(45, color='lime', linestyle=':', linewidth=1.5)
    ax.legend(loc='upper right', framealpha=0.85)
    cb_c = fig.colorbar(im_cake, ax=ax, fraction=0.03, pad=0.02)
    cb_c.set_label(r"Intensity $I(q, \chi)\ [\mathrm{cm}^{-1}]$")
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "step5_azimuthal_cake_plot.png"), dpi=300)
    fig.savefig(os.path.join(OUTPUT_DIR, "step5_azimuthal_cake_plot.pdf"))
    plt.close(fig)

    # Step 6: 1D Absolute Intensity I(q) with Uncertainty Bands & Relative Error
    print("  Generating Fig 6: 1D Profile with Propagated Uncertainty Bands...")
    fig, (ax_main, ax_diff, ax_err) = plt.subplots(3, 1, figsize=(9.5, 9.5), sharex=True, 
                                                   gridspec_kw={'height_ratios': [3.0, 1.0, 1.2]})

    # Main plot: I(q) with 1-sigma uncertainty ribbon
    ax_main.plot(q_ref, i_ref, color='#0072B2', linewidth=1.8, label=r'Complete Radial Average ($0^\circ - 360^\circ$)')
    ax_main.fill_between(q_ref, i_ref - sigma_ref, i_ref + sigma_ref, color='#0072B2', alpha=0.25, 
                         label=r'Propagated Uncertainty ($\pm 1\sigma_I$)')
    ax_main.plot(q_meridian, i_meridian, color='#E69F00', linewidth=1.3, linestyle='--', label=r'Meridian Sector ($45^\circ - 135^\circ$)')
    ax_main.plot(q_equatorial, i_equatorial, color='#009E73', linewidth=1.3, linestyle=':', label=r'Equatorial Sector ($-45^\circ - 45^\circ$)')

    ax_main.set_yscale('log')
    ax_main.set_xscale('log')
    ax_main.set_xlim(0.09, 42.0)
    ax_main.set_ylim(5e-4, 5e1)
    ax_main.set_ylabel(r"Cross-Section $\frac{d\Sigma}{d\Omega}(q)\ [\mathrm{cm}^{-1}]$")
    ax_main.set_title(rf"$\mathbf{{Absolute\ 1D\ SAXS\ Profile\ with\ Propagated\ Uncertainty}}\ I(q) \pm \sigma(q)$", pad=10)
    ax_main.grid(True, which="both", ls="-", alpha=0.2)
    ax_main.legend(loc='lower left', framealpha=0.9)

    # Residual against PyFAI integration
    diff_val = i_all - i_ref
    ax_diff.plot(q_all, diff_val, color='black', linewidth=1.0)
    ax_diff.fill_between(q_all, -sigma_ref, sigma_ref, color='gray', alpha=0.15, label=r'$\pm 1\sigma$ Band')
    ax_diff.axhline(0, color='red', linestyle='--', linewidth=0.8)
    ax_diff.set_xscale('log')
    ax_diff.set_xlim(0.09, 42.0)
    ax_diff.set_ylim(-0.02, 0.02)
    ax_diff.set_ylabel(r"Residual $[\mathrm{cm}^{-1}]$")
    ax_diff.grid(True, which="both", ls="-", alpha=0.2)
    ax_diff.legend(loc='upper right', framealpha=0.8, fontsize=8)

    # Relative Fractional Uncertainty: sigma(q) / |I(q)| in percent
    rel_err_perc = (sigma_ref / np.maximum(np.abs(i_ref), 1e-6)) * 100.0
    ax_err.plot(q_ref, rel_err_perc, color='#CC79A7', linewidth=1.3)
    ax_err.axhline(5.0, color='darkgreen', linestyle=':', linewidth=1.0, label=r'5% Relative Error Threshold')
    ax_err.set_xscale('log')
    ax_err.set_yscale('log')
    ax_err.set_xlim(0.09, 42.0)
    ax_err.set_ylim(0.1, 100.0)
    ax_err.set_xlabel(r"Scattering Vector $q = \frac{4\pi}{\lambda}\sin\theta\ [\mathrm{nm}^{-1}]$")
    ax_err.set_ylabel(r"Relative Error $\frac{\sigma_I(q)}{|I(q)|}\ [\%]$")
    ax_err.grid(True, which="both", ls="-", alpha=0.2)
    ax_err.legend(loc='upper left', framealpha=0.8, fontsize=8)

    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "step6_1d_integrated_profiles.png"), dpi=300)
    fig.savefig(os.path.join(OUTPUT_DIR, "step6_1d_integrated_profiles.pdf"))
    plt.close(fig)

    # Master Multi-Panel Publication Figure (incorporating uncertainty ribbon in panel F)
    fig_pub = plt.figure(figsize=(16, 10.5))
    gs = gridspec.GridSpec(2, 3, figure=fig_pub, hspace=0.26, wspace=0.28)
    axA = fig_pub.add_subplot(gs[0, 0])
    axB = fig_pub.add_subplot(gs[0, 1])
    axC = fig_pub.add_subplot(gs[0, 2])
    axD = fig_pub.add_subplot(gs[1, 0])
    axE = fig_pub.add_subplot(gs[1, 1])
    axF = fig_pub.add_subplot(gs[1, 2])

    axA.set_facecolor('#111827')
    cmap_A = plt.cm.inferno.copy()
    cmap_A.set_bad('#111827')
    imA = axA.imshow(raw_disp_sample, origin='lower', cmap=cmap_A, norm=LogNorm(vmin=1, vmax=5e4))
    mask_cyan_pub = ListedColormap(['#00e5ff'])
    mask_cyan_pub.set_bad(alpha=0.0)
    axA.imshow(np.where(mask_full, 1.0, np.nan), origin='lower', cmap=mask_cyan_pub, alpha=1.0, interpolation='nearest')
    axA.plot(xcen, ycen, 'r+', markersize=12, markeredgewidth=2)
    axA.set_title(r"$\mathbf{(a)}$ Raw Detector Image & Mask (Cyan)", pad=8)
    axA.set_xlabel("Pixel X")
    axA.set_ylabel("Pixel Y")
    cbA = fig_pub.colorbar(imA, ax=axA, fraction=0.046, pad=0.04)
    cbA.set_label(r"Photon Counts $[\mathrm{log}_{10}]$", fontsize=9)

    axB.set_facecolor('#111827')
    cmap_B = plt.cm.inferno.copy()
    cmap_B.set_bad('#111827')
    imB = axB.imshow(norm_sample_disp, origin='lower', cmap=cmap_B, norm=LogNorm(vmin=1e-6, vmax=1e-1))
    axB.set_title(r"$\mathbf{(b)}$ Beam Monitor Normalized ($I / I_0$)", pad=8)
    axB.set_xlabel("Pixel X")
    axB.set_ylabel("Pixel Y")
    cbB = fig_pub.colorbar(imB, ax=axB, fraction=0.046, pad=0.04)
    cbB.set_label(r"Normalized Counts [$I_0^{-1}$]", fontsize=9)

    axC.set_facecolor('#111827')
    cmap_C = plt.cm.inferno.copy()
    cmap_C.set_bad('#111827')
    imC = axC.imshow(norm_bkg_disp, origin='lower', cmap=cmap_C, norm=LogNorm(vmin=1e-6, vmax=1e-1))
    axC.set_title(r"$\mathbf{(c)}$ Capillary Background ($I_{bkg} / I_{0,bkg}$)", pad=8)
    axC.set_xlabel("Pixel X")
    axC.set_ylabel("Pixel Y")
    cbC = fig_pub.colorbar(imC, ax=axC, fraction=0.046, pad=0.04)
    cbC.set_label(r"Background Counts [$I_{0,bkg}^{-1}$]", fontsize=9)

    axD.set_facecolor('#111827')
    cmap_D = plt.cm.coolwarm.copy()
    cmap_D.set_bad('#111827')
    imD = axD.imshow(img_sub_2d_masked, origin='lower', cmap=cmap_D,
                    norm=SymLogNorm(linthresh=0.01, linscale=1.0, vmin=-0.2, vmax=5.0))
    axD.set_title(r"$\mathbf{(d)}$ Net Sample 2D Pattern ($I_{sub}$)", pad=8)
    axD.set_xlabel("Pixel X")
    axD.set_ylabel("Pixel Y")
    cbD = fig_pub.colorbar(imD, ax=axD, fraction=0.046, pad=0.04)
    cbD.set_label(r"Calibrated Intensity $[\mathrm{cm}^{-1}]$", fontsize=9)

    imE = axE.imshow(cake_disp, origin='lower', aspect='auto', cmap='inferno',
                     extent=[q_cake.min(), q_cake.max(), chi_cake.min(), chi_cake.max()],
                     norm=LogNorm(vmin=1e-2, vmax=10))
    axE.set_title(r"$\mathbf{(e)}$ Azimuthal Transformation $I(q, \chi)$", pad=8)
    axE.set_xlabel(r"$q\ [\mathrm{nm}^{-1}]$")
    axE.set_ylabel(r"$\chi\ [^\circ]$")
    cbE = fig_pub.colorbar(imE, ax=axE, fraction=0.046, pad=0.04)
    cbE.set_label(r"$I(q, \chi)\ [\mathrm{cm}^{-1}]$", fontsize=9)

    # Panel F with shaded uncertainty ribbon
    axF.plot(q_ref, i_ref, color='#1f77b4', linewidth=1.8, label=r'All ($0^\circ - 360^\circ$)')
    axF.fill_between(q_ref, i_ref - sigma_ref, i_ref + sigma_ref, color='#1f77b4', alpha=0.25, label=r'$\pm 1\sigma_I$ Error Band')
    axF.plot(q_meridian, i_meridian, color='#ff7f0e', linewidth=1.3, linestyle='--', label=r'Meridian ($45^\circ - 135^\circ$)')
    axF.plot(q_equatorial, i_equatorial, color='#2ca02c', linewidth=1.3, linestyle=':', label=r'Equatorial ($-45^\circ - 45^\circ$)')
    axF.set_yscale('log')
    axF.set_xscale('log')
    axF.set_xlim(0.09, 42.0)
    axF.set_ylim(1e-3, 4e1)
    axF.set_title(r"$\mathbf{(f)}$ Calibrated 1D Profile with Uncertainty", pad=8)
    axF.set_xlabel(r"Scattering Vector $q\ [\mathrm{nm}^{-1}]$")
    axF.set_ylabel(r"Cross-Section $\frac{d\Sigma}{d\Omega}(q)\ [\mathrm{cm}^{-1}]$")
    axF.grid(True, which="both", ls="-", alpha=0.2)
    axF.legend(loc='lower left', framealpha=0.9, fontsize=8.5)

    fig_pub.savefig(os.path.join(OUTPUT_DIR, "Figure_SAXS_Reduction_Pipeline.png"), dpi=300)
    fig_pub.savefig(os.path.join(OUTPUT_DIR, "Figure_SAXS_Reduction_Pipeline.pdf"))
    plt.close(fig_pub)
    print(f"  Step figures, 2D uncertainty maps, and master figure generated successfully in {OUTPUT_DIR}")


def run_time_series_analysis(start_scan=184, end_scan=292):
    """
    Constructs the 16-hour in-situ kinetic evolution across scan #184 to #292,
    including transmission T(t), 2D contour map I(q, t), and stacked waterfall profiles.
    """
    print(f"\n--- Running In-Situ Kinetic Time Series Analysis (Scans #{start_scan} to #{end_scan}) ---")
    time_series = []

    # Read from Exp2_1 HDF5
    if os.path.exists(REDUCED_EXP2_1):
        with h5py.File(REDUCED_EXP2_1, 'r') as h1:
            for gname in h1.keys():
                s = h1[gname].attrs.get('SPEC Scan - Data')
                if s and start_scan <= int(s) <= min(274, end_scan):
                    epoch = float(h1[gname].attrs.get('Epoch - Data', 0))
                    trans = float(h1[gname].attrs.get('Transmission - Data', '0').split()[0])
                    time_series.append((int(s), gname, epoch, trans, REDUCED_EXP2_1))

    # Read from Exp2_2 HDF5 if range extends beyond 274
    if end_scan > 274 and os.path.exists(REDUCED_EXP2_2):
        with h5py.File(REDUCED_EXP2_2, 'r') as h2:
            for gname in h2.keys():
                s = h2[gname].attrs.get('SPEC Scan - Data')
                if s and 274 < int(s) <= end_scan:
                    epoch = float(h2[gname].attrs.get('Epoch - Data', 0))
                    trans = float(h2[gname].attrs.get('Transmission - Data', '0').split()[0])
                    time_series.append((int(s), gname, epoch, trans, REDUCED_EXP2_2))

    if not time_series:
        print("  No reduced measurements found in the specified range.")
        return

    time_series.sort(key=lambda x: x[2])
    n_points = len(time_series)
    epochs = np.array([x[2] for x in time_series])
    t0 = epochs[0]
    time_hours = (epochs - t0) / 3600.0
    transmissions = np.array([x[3] for x in time_series])

    print(f"  Total time points: {n_points} over {time_hours[-1]:.2f} hours")

    # Downsample time points for matrix rendering (e.g. 200 time steps) to keep file compact
    n_samples = min(300, n_points)
    sample_indices = np.linspace(0, n_points - 1, n_samples).astype(int)
    t_sampled = time_hours[sample_indices]
    trans_sampled = transmissions[sample_indices]

    # Load 1D curves
    curves = []
    q_grid = None
    for idx in sample_indices:
        item = time_series[idx]
        file_path, gname = item[4], item[1]
        with h5py.File(file_path, 'r') as h:
            prof = h[f"{gname}/profile"][:]
            if q_grid is None:
                q_grid = prof[:, 0]
            curves.append(prof[:, 1])

    I_matrix = np.array(curves)  # shape: (n_samples, n_q)

    # Plot 4-Panel In-Situ Kinetic Publication Figure
    fig = plt.figure(figsize=(15, 10))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.28, wspace=0.25)
    ax_trans = fig.add_subplot(gs[0, 0])
    ax_contour = fig.add_subplot(gs[0, 1])
    ax_waterfall = fig.add_subplot(gs[1, 0])
    ax_invariant = fig.add_subplot(gs[1, 1])

    # 1. Transmission vs Time
    ax_trans.plot(time_hours, transmissions, color='#0072B2', alpha=0.4, linewidth=0.8)
    ax_trans.plot(t_sampled, trans_sampled, color='#0072B2', linewidth=2.0, label=r'Relative Transmission $T(t)$')
    ax_trans.set_title(r"$\mathbf{(a)}$ Sample Transmission Evolution", pad=8)
    ax_trans.set_xlabel("Elapsed Time [hours]")
    ax_trans.set_ylabel(r"Transmission Factor $T_{rel}$")
    ax_trans.set_ylim(0.60, 0.72)
    ax_trans.grid(True, which='both', ls='-', alpha=0.2)
    ax_trans.legend(loc='upper right', framealpha=0.85)

    # 2. 2D Kinetic Map I(q, t)
    I_disp = np.clip(I_matrix, 1e-4, 50.0)
    im_k = ax_contour.imshow(I_disp, origin='lower', aspect='auto', cmap='inferno',
                             extent=[q_grid.min(), q_grid.max(), t_sampled.min(), t_sampled.max()],
                             norm=LogNorm(vmin=1e-3, vmax=10))
    ax_contour.set_title(r"$\mathbf{(b)}$ Time-Resolved SAXS Map $I(q, t)$", pad=8)
    ax_contour.set_xlabel(r"Scattering Vector $q\ [\mathrm{nm}^{-1}]$")
    ax_contour.set_ylabel("Elapsed Time [hours]")
    ax_contour.set_xlim(0.1, 20.0)
    cb_k = fig.colorbar(im_k, ax=ax_contour, fraction=0.046, pad=0.04)
    cb_k.set_label(r"$I(q, t)\ [\mathrm{cm}^{-1}]$")

    # 3. Stacked Waterfall 1D Profiles at Representative Time Intervals
    time_slices = np.linspace(0, len(t_sampled) - 1, 6).astype(int)
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(time_slices)))
    for i, s_idx in enumerate(time_slices):
        t_val = t_sampled[s_idx]
        ax_waterfall.plot(q_grid, curves[s_idx], color=colors[i], linewidth=1.5,
                          label=rf"$t = {t_val:.1f}\ \mathrm{{h}}$ (#{time_series[sample_indices[s_idx]][0]})")
    ax_waterfall.set_yscale('log')
    ax_waterfall.set_xscale('log')
    ax_waterfall.set_xlim(0.1, 35.0)
    ax_waterfall.set_ylim(1e-3, 30.0)
    ax_waterfall.set_title(r"$\mathbf{(c)}$ Representative SAXS Curves Over Time", pad=8)
    ax_waterfall.set_xlabel(r"Scattering Vector $q\ [\mathrm{nm}^{-1}]$")
    ax_waterfall.set_ylabel(r"Intensity $\frac{d\Sigma}{d\Omega}(q)\ [\mathrm{cm}^{-1}]$")
    ax_waterfall.grid(True, which='both', ls='-', alpha=0.2)
    ax_waterfall.legend(loc='lower left', framealpha=0.9, fontsize=9)

    # 4. Integrated SAXS Scattering Invariant Q(t)
    q_mask = (q_grid >= 0.2) & (q_grid <= 5.0)
    q_sub = q_grid[q_mask]
    trap_func = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
    invariant = np.array([trap_func(c[q_mask] * (q_sub**2), q_sub) for c in curves])
    ax_invariant.plot(t_sampled, invariant, color='#D55E00', linewidth=2.0, marker='o', markersize=3, label=r'Integrated Invariant $Q(t) = \int q^2 I(q) dq$')
    ax_invariant.set_title(r"$\mathbf{(d)}$ Structural Scattering Invariant vs Time", pad=8)
    ax_invariant.set_xlabel("Elapsed Time [hours]")
    ax_invariant.set_ylabel(r"Scattering Invariant $Q\ [\mathrm{cm}^{-1}\cdot\mathrm{nm}^{-3}]$")
    ax_invariant.grid(True, which='both', ls='-', alpha=0.2)
    ax_invariant.legend(loc='upper right', framealpha=0.85)

    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, "step7_time_series_evolution_184_to_292.png"), dpi=300)
    fig.savefig(os.path.join(OUTPUT_DIR, "step7_time_series_evolution_184_to_292.pdf"))
    plt.close(fig)
    print(f"  Time series evolution figure generated: step7_time_series_evolution_184_to_292.png / .pdf")


def main():
    parser = argparse.ArgumentParser(description="Publication SAXS Reduction Pipeline for Exp 2_1")
    parser.add_argument("--scan", type=int, default=184, help="Target measurement scan number (default: 184)")
    parser.add_argument("--frame", type=int, default=1, help="Target frame number within scan (default: 1)")
    parser.add_argument("--bkg", type=int, default=63, help="Background scan number (default: 63)")
    parser.add_argument("--series", nargs=2, type=int, metavar=('START', 'END'), default=None,
                        help="Generate kinetic time series across scan range, e.g. --series 184 292")
    parser.add_argument("--all", action="store_true", help="Execute both single-scan pipeline and full time-series")
    args = parser.parse_args()

    print("=" * 70)
    print("SAXS DATA REDUCTION PIPELINE & PUBLICATION SUITE (WITH UNCERTAINTIES)")
    print("=" * 70)

    if args.all or (args.series is None):
        run_single_scan_reduction(scan_num=args.scan, frame_num=args.frame, bkg_scan=args.bkg)

    if args.series is not None:
        start_s = args.series[0]
        end_s = args.series[1]
        run_time_series_analysis(start_scan=start_s, end_scan=end_s)

    print("\n" + "=" * 70)
    print("ALL PROCESSING COMPLETE. Figures available in:")
    print(f"  {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == '__main__':
    main()
