# Spectral Transfer Function Calculator

## Description
`rtf_H-He.f90` computes the hydrogen/helium spectral transfer function, the primary product presented in the accompanying research article on Cosmic Dawn and Reionization. The code integrates observational constraints, photoionization cross-sections, and Gunn–Peterson optical depths to produce frequency-dependent transmission tables across a grid of redshifts.

## Contents Referenced in the Article
1. Observational data on the reionization history
2. Optical depth of the diffuse intergalactic gas in hydrogen lines and continuum
3. Optical depth of the diffuse intergalactic gas in the helium lines and continuum
4. Spectral flux from halos of Cosmic Dawn and Reionization epochs
5. Conclusions

## How to Use
1. Ensure the required atomic data files exist under `../AtomicData/HOS.dat` and `../AtomicData/HeIOS.dat` relative to the executable.
2. Compile the code:
   - **macOS**: `gfortran -isysroot $(xcrun --sdk macosx --show-sdk-path) -O3 -o rtf ./rtf_H-He.f90`
   - **Linux**: `gfortran -O3 -o rtf rtf_H-He.f90` or `ifort -O3 -o rtf rtf_H-He.f90`
   - **Windows (MinGW/Intel)**: `gfortran -O3 -o rtf.exe rtf_H-He.f90`
3. Run the executable from the directory containing `../AtomicData/`:
   ```bash
   ./rtf
   ```
4. Inspect the generated spectral transfer function tables named `rtf_H-He_2_z_*.dat`, each corresponding to one of the predefined redshift nodes.

### Sample Session
```bash
gfortran -O3 -o rtf ./rtf_H-He.f90
./rtf
ls rtf_H-He_2_z_*.dat
head -n 5 rtf_H-He_2_z_7.0.dat
```
This produces frequency vs. transmission values that can be ingested into analysis notebooks or plotting scripts.

## Outputs
- `rtf_H-He_2_z_<z>.dat`: Spectral transfer function tables (frequency in 10^12 Hz vs. transmitted flux).
- `gff.dat`: Auxiliary diagnostics written during runtime.

## Authors
- Bohdan Novosyadlyj <bnovos@gmail.com>
- Pavlo Kopach <luspav@gmail.com>

## License
GNU General Public License v3.0

Copyright (c) 2025 Bohdan Novosyadlyj and Pavlo Kopach

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. See `LICENSE` for the full text of the license and your obligations
when redistributing or modifying the code.
