Understanding Metastability in SAR ADCs: Part II: Asynchronous 
==============================================================
MATLAB code used for models described in article: A. Yu, D. Bankman, K. Zheng, B. Murmann, "Understanding Metastability in SAR ADCs: Part II: Asynchronous," *IEEE Solid-State Circuits Magazine* 11, 3 (2019) https://ieeexplore.ieee.org/document/8811772

All code written using MATLAB R2016a.

Folder Structure and File Descriptions
---------------------------------------
```
figs/                                       - svg files for some figures and cover art
src/
 ├─ results/                                - folder to store generated results
 ├─ asar_meta_pmf_ideal.m                   - output pmf for ideal noiseless async SAR
 ├─ asar_meta_pmf_noise_combined.m          - output pmf with analytic thermal noise
 ├─ asar_meta_pmf_noise_combined_parallel.m - same as above but parallelized compute with parfor
 ├─ figure_settings.m                       - common figure/plotting settings
 ├─ plot_helper_functions.m                 - plot analytic helper functions (DAC, residuals, output codes,...)
 ├─ plot_pmeta_async_vs_sync.m              - compare Pmeta of async vs. sync
 ├─ plot_pmf_ideal_reset.m                  - ideal ASAR PMF vs. reset code
 ├─ plot_pmf_ideal_tfix.m                   - ideal ASAR PMF vs. TFIX and TSAR
 ├─ plot_pmf_noise_vs_N.m                   - main script for sweeping/plotting error PMFs with thermal noise
 ├─ plot_t_avg_easy.m                       - compare average vs easy regeneration
 ├─ plot_timing_bands.m                     - generate timing band outlines for Fig. 4b
 ├─ t_avg.m                                 - average regeneration time function
 ├─ t_easy.m                                - "easy" regeneration time function
 ├─ test_asar_meta_noise_combined.m         - test combined noise metastability PMF and PDF at single point
 └─ test_pdf_noise_meta.m                   - analyzing noise metastability PDF behavior
```