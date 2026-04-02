# Spectral Transfer Function for H and He

## Description

`rtf_H-He.f90` computes the frequency-dependent spectral transfer function **F(nu, z)** for radiation propagating through the intergalactic medium (IGM) during cosmic reionization. The code integrates contributions from:

- **HI** photoionization continuum and Lyman-series lines (Ly-alpha through Ly-39)
- **HeI** photoionization continuum and resonance lines (10 transitions)
- **HeII** photoionization continuum and Lyman-series lines (hydrogen-like, Z=2)
- **Thomson scattering** by free electrons (computed, optionally included)

The neutral hydrogen fraction xHI(z) is modelled with piecewise analytic fits from Bosman, Becker & Inoue (2022), extended to high redshift by either Kageura & Bosman (2025) or Glazer & Bosman (2018), selectable via the `klxHI` parameter.

Cosmological parameters follow Planck 2018 values.

## Contents Referenced in the Article

1. Observational data on the reionization history
2. Optical depth of the diffuse intergalactic gas in hydrogen lines and continuum
3. Optical depth of the diffuse intergalactic gas in the helium lines and continuum
4. Spectral flux from halos of Cosmic Dawn and Reionization epochs
5. Conclusions

## Input Files

The following atomic data files must be present in the working directory:

| File | Contents | Format |
|------|----------|--------|
| `HOS.dat` | HI Lyman-series oscillator strengths (39 lines) | index, wavelength [A], f |
| `HeIOS.dat` | HeI line oscillator strengths (10 lines) | index, wavelength [A], f |

## How to Use

1. Place `HOS.dat` and `HeIOS.dat` in the same directory as the executable.
2. Compile:
   - **macOS**: `gfortran -isysroot $(xcrun --sdk macosx --show-sdk-path) -O3 -o rtf ./rtf_H-He.f90`
   - **Linux**: `gfortran -O3 -o rtf rtf_H-He.f90` or `ifort -O3 -o rtf rtf_H-He.f90`
   - **Windows**: `gfortran -O3 -o rtf.exe rtf_H-He.f90`
3. Run:
   ```bash
   ./rtf
   ```

### Sample Session

```bash
gfortran -O3 -o rtf ./rtf_H-He.f90
./rtf
ls rtf_H-He_1_z_*.dat
head -n 5 rtf_H-He_1_z_7.0.dat
```

## Output

For each redshift node, the program writes a two-column table:

- **File pattern**: `rtf_H-He_1_z_<z>.dat`
- **Columns**: observed frequency [THz] | F(nu)

Default redshift nodes: z = 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0, 12.0, 15.0

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `klxHI` | 1 | xHI model: 1 = Kageura+Bosman (2025), 2 = Glazer+Bosman (2018) |
| `n_xf` | 32000 | Frequency grid points (log-spaced, 100–30000 THz) |
| `n_z` | 9 | Number of source redshifts |
| `zo` | 0.0 | Observer redshift |

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
