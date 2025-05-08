# Prototype Filter Design for WOLA Filter Bank

This module provides MATLAB code for designing **distortion-minimizing** and **norm-minimizing**  synthesis prototype filter for **Weighted Overlap-Add (WOLA)** filter banks.

### Reference
Please cite the following paper if you use this code in your work:
> M. Sharma and M. Moonen, "Prototype filter design for weighted overlap-add filter bank based sub-band adaptive filtering applications," 2023 31st European Signal Processing Conference (EUSIPCO), Helsinki, Finland, 2023, pp. 366-370.

**Download:**  [IEEE Xplore](https://ieeexplore.ieee.org/document/10289725) | [Open Access Version (EURASIP)](https://eurasip.org/Proceedings/Eusipco/Eusipco2023/pdfs/0000366.pdf)  

### Key Files:
- `design_wola_synthesis_window.m` – Core function to generate the synthesis window given an analysis window and design criterion.
- `demo_design_wola_synthesis_window.m` – Example script demonstrating usage and visual comparison of window designs.