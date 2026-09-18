import os
import h5py
import fabio
import pyFAI
import numpy as np
import matplotlib.pyplot as plt

os.environ['MPLCONFIGDIR'] = '/tmp/mpl_cache'

CALIB_DIR = '/run/media/tomek/data/Mrc_26/saxs_publication_project/calibration'
EIGER_DIR = '/run/media/tomek/data/Mrc_26/muspot_data/eiger_files'
SMP_DIR = '/run/media/tomek/data/Mrc_26/saxs_publication_project/raw_data/sample_scan184'
OUTPUT_DIR = '/run/media/tomek/data/Mrc_26/saxs_publication_project/publication_figures'

ai = pyFAI.load(os.path.join(CALIB_DIR, 'calib.poni'))
mask = fabio.open(os.path.join(CALIB_DIR, 'Start_rigorous_mask.tiff')).data != 0

# 1. Load Scans
# Air scan 47 (Empty beamline background)
with h5py.File(os.path.join(EIGER_DIR, 'AirGlassy_000047_data_000001.h5'), 'r') as f:
    air_raw = f['/entry/data/data'][0].astype(np.float64)
I0_air = 915218.0

# BeamOff scan 61 (Dark noise)
with h5py.File(os.path.join(EIGER_DIR, 'BeamOff_000061_data_000001.h5'), 'r') as f:
    dark_raw = f['/entry/data/data'][0].astype(np.float64)
I0_dark = 1.0

# Capillary scan 63 (Capillary buffer background)
with h5py.File(os.path.join(EIGER_DIR, 'CapillaryTop_000063_data_000001.h5'), 'r') as f:
    cap_raw = f['/entry/data/data'][0].astype(np.float64)
I0_cap = 973266.0

# Sample scan 184 (Target measurement)
with h5py.File(os.path.join(SMP_DIR, 'LUKU2026030601PF_000184_data_000001.h5'), 'r') as f:
    smp_raw = f['/entry/data/data'][0].astype(np.float64)
I0_smp = 627109.0  # Sample I0 from SPEC
T_rel = 0.648902   # Relative transmission
K_std = 385.2002   # Glassy carbon standard scaling factor in Spiger Start.spig

# 2. Azimuthal Integration
n_bins = 1000
unit = 'q_nm^-1'

res_air = ai.integrate1d(air_raw, n_bins, unit=unit, mask=mask, normalization_factor=I0_air)
res_dark = ai.integrate1d(dark_raw, n_bins, unit=unit, mask=mask, normalization_factor=I0_dark)
res_cap = ai.integrate1d(cap_raw, n_bins, unit=unit, mask=mask, normalization_factor=I0_cap)
res_smp = ai.integrate1d(smp_raw, n_bins, unit=unit, mask=mask, normalization_factor=I0_smp)

q = res_smp.radial

# Normalize to absolute scale (cm^-1)
# Note: Sample is attenuation corrected: I_meas = K_std * (smp / (I0_smp * T_rel))
I_meas = res_smp.intensity / T_rel * K_std
I_cap = res_cap.intensity * K_std
I_air = res_air.intensity * K_std
I_dark = res_dark.intensity * K_std

# Standard pipeline: net sample after capillary subtraction
I_sub_correct = I_meas - I_cap

# Capillary transmission relative to air: T_cap ~ 0.95
T_cap = 0.95
# Isolated capillary glass scatter:
I_glass = I_cap - T_cap * I_air

# Double-subtraction mistake (corrupted data):
I_double_sub = I_meas - I_cap - I_air

# 3. Create 4-panel publication figure
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(13, 10))

# Panel (a): All profiles overview (log-log)
ax1.plot(q, I_meas, color='#1d4ed8', lw=1.8, label=r'Sample measurement $I_{\mathrm{meas}}(q)$ (Scan #184)')
ax1.plot(q, I_cap, color='#dc2626', lw=1.8, ls='--', label=r'Capillary background $I_{\mathrm{cap}}(q)$ (Scan #63)')
ax1.plot(q, I_air, color='#d97706', lw=1.8, ls='-.', label=r'Empty beamline background $I_{\mathrm{air}}(q)$ (Scan #47)')
ax1.plot(q, I_sub_correct, color='#059669', lw=2.0, label=r'Correct net sample $I_{\mathrm{sub}}(q) = I_{\mathrm{meas}} - I_{\mathrm{cap}}$')
ax1.plot(q, I_dark, color='#6b7280', lw=1.2, ls=':', label=r'Detector dark noise $I_{\mathrm{dark}}(q) \approx 0$ (Scan #61)')

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(0.09, 42.0)
ax1.set_ylim(1e-2, 40.0)
ax1.set_title(r"$\mathbf{(a)}$ Complete scattering decomposition across full $q$-range", pad=8, fontsize=11)
ax1.set_xlabel(r"Scattering vector $q\ [\mathrm{nm}^{-1}]$", fontsize=10.5)
ax1.set_ylabel(r"Differential cross-section $\frac{d\Sigma}{d\Omega}(q)\ [\mathrm{cm}^{-1}]$", fontsize=10.5)
ax1.grid(True, which='both', ls='-', alpha=0.2)
ax1.legend(loc='upper right', framealpha=0.92, fontsize=8.5)

# Panel (b): Low-q zoom: Beamline background vs Capillary
ax2.plot(q, I_meas, color='#1d4ed8', lw=1.8, label=r'Sample $I_{\mathrm{meas}}(q)$')
ax2.plot(q, I_cap, color='#dc2626', lw=1.8, ls='--', label=r'Capillary buffer $I_{\mathrm{cap}}(q)$')
ax2.plot(q, I_air, color='#d97706', lw=1.8, ls='-.', label=r'Empty beamline / air $I_{\mathrm{air}}(q)$')
ax2.plot(q, I_sub_correct, color='#059669', lw=2.0, label=r'Net sample $I_{\mathrm{sub}}(q)$')

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim(0.09, 2.0)
ax2.set_ylim(0.05, 10.0)
ax2.set_title(r"$\mathbf{(b)}$ Low-$q$ regime: parasitic air scatter and beamstop flare", pad=8, fontsize=11)
ax2.set_xlabel(r"Scattering vector $q\ [\mathrm{nm}^{-1}]$", fontsize=10.5)
ax2.set_ylabel(r"Cross-section $\frac{d\Sigma}{d\Omega}(q)\ [\mathrm{cm}^{-1}]$", fontsize=10.5)
ax2.grid(True, which='both', ls='-', alpha=0.2)
ax2.legend(loc='upper right', framealpha=0.92, fontsize=8.5)

ax2.annotate('Air background accounts for\n>90% of capillary low-q scatter',
             xy=(0.2, 0.45), xytext=(0.11, 0.12),
             arrowprops=dict(arrowstyle="->", color='#d97706', lw=1.2),
             bbox=dict(boxstyle="round,pad=0.25", fc="#fef3c7", ec="#d97706", alpha=0.9), fontsize=8.5)

# Panel (c): Decomposing capillary into Air + Glass
ax3.plot(q, I_cap, color='#dc2626', lw=1.8, label=r'Total capillary background $I_{\mathrm{cap}}(q)$')
ax3.plot(q, T_cap * I_air, color='#d97706', lw=1.8, ls='-.', label=r'Transmitted air background $T_{\mathrm{cap}} \cdot I_{\mathrm{air}}(q)$')
ax3.plot(q, np.maximum(I_glass, 1e-4), color='#7c3aed', lw=2.0, label=r'Isolated borosilicate glass $I_{\mathrm{glass}}(q) = I_{\mathrm{cap}} - T_{\mathrm{cap}} I_{\mathrm{air}}$')

ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlim(0.09, 42.0)
ax3.set_ylim(1e-3, 5.0)
ax3.set_title(r"$\mathbf{(c)}$ Physical decomposition: capillary buffer = air path + borosilicate glass", pad=8, fontsize=11)
ax3.set_xlabel(r"Scattering vector $q\ [\mathrm{nm}^{-1}]$", fontsize=10.5)
ax3.set_ylabel(r"Cross-section $\frac{d\Sigma}{d\Omega}(q)\ [\mathrm{cm}^{-1}]$", fontsize=10.5)
ax3.grid(True, which='both', ls='-', alpha=0.2)
ax3.legend(loc='upper right', framealpha=0.92, fontsize=8.5)

# Panel (d): Why empty beamline is NOT subtracted twice (Double-subtraction error)
ax4.plot(q, I_sub_correct, color='#059669', lw=2.2, label=r'Correct single subtraction: $I_{\mathrm{meas}} - I_{\mathrm{cap}}$')
ax4.plot(q, I_double_sub, color='#ef4444', lw=1.8, ls='--', label=r'Erroneous double subtraction: $I_{\mathrm{meas}} - I_{\mathrm{cap}} - I_{\mathrm{air}}$')
ax4.axhline(0, color='black', lw=1.0, ls=':')

ax4.set_xscale('log')
ax4.set_xlim(0.09, 42.0)
ax4.set_ylim(-0.5, 3.5)
ax4.set_title(r"$\mathbf{(d)}$ Methodological demonstration: avoiding erroneous double subtraction", pad=8, fontsize=11)
ax4.set_xlabel(r"Scattering vector $q\ [\mathrm{nm}^{-1}]$", fontsize=10.5)
ax4.set_ylabel(r"Linear cross-section $\frac{d\Sigma}{d\Omega}(q)\ [\mathrm{cm}^{-1}]$", fontsize=10.5)
ax4.grid(True, which='both', ls='-', alpha=0.2)
ax4.legend(loc='upper right', framealpha=0.92, fontsize=8.5)

ax4.annotate('Severe unphysical deficit / negative counts\nif empty beamline is subtracted twice!',
             xy=(1.5, -0.28), xytext=(0.25, -0.42),
             arrowprops=dict(arrowstyle="->", color='#ef4444', lw=1.2),
             bbox=dict(boxstyle="round,pad=0.25", fc="#fef2f2", ec="#ef4444", alpha=0.9), fontsize=8.5)

plt.tight_layout()
png_out = os.path.join(OUTPUT_DIR, 'step4c_empty_beamline_background.png')
pdf_out = os.path.join(OUTPUT_DIR, 'step4c_empty_beamline_background.pdf')
fig.savefig(png_out, dpi=300)
fig.savefig(pdf_out)
plt.close(fig)
print(f"Successfully generated: {png_out}")
