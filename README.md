# GeneralizedWOLA-SystemIdentification
**GeneralizedWOLA-SystemIdentification** is a MATLAB implementation of a generalized Weighted Overlap-Add (WOLA) filter bank framework for **subband system identification**, specifically tailored for **acoustic echo cancellation (AEC)** task.

The project includes support for both:
- Full-complexity **Sliding Windowed DFT (SWDFT)** implementation
- Low-complexity **Per-Tone WOLA (PT-WOLA)** implementation

> **License Notice**  
> This repository contains code associated with unpublished research. Redistribution, modification, or reuse is **not permitted** until the corresponding paper is accepted and a formal open-source license is applied.  
> For collaboration or access inquiries, please contact the author.  
> See `LICENSE_PREPRINT.txt` for current restrictions.

---

## How to Run

Use the main script from the project root folder:

```matlab
run_WOLA_AEC_Simulation.m
```
---

## Overview of Key Files

| File                             | Description                                                   |
|----------------------------------|---------------------------------------------------------------|
| `run_WOLA_AEC_Simulation.m`      | Entry-point script to run AEC simulation                      |
| `Generalized_WOLA_AEC_System.m`  | Main controller class                                         |
| `AECParameters.m`                | Global simulation parameters                                  |
| `PTWOLAConfig.m`                 | Configuration for generalized WOLA filter bank setup          |
| `RIRGenerator.m`                 | Room impulse response generator                               |
| `SignalGenerator.m`              | Builds signals (far-end, near-end, echo, mic)                 |
| `WindowDesigner.m`               | Creates analysis and synthesis windows                        |
| `AECSimulation.m`                | Runs the simulation loop                   |
| `Generalized_WOLA_AEC_Processor.m` | Core adaptive filtering logic                               |
| `ResultsManager.m`               | Plots and manages results                                |
| `design_wola_synthesis_window.m` | Synthesis window design function                              |

