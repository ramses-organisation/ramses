# Hot Rotating CGM Inflow Patch for RAMSES

This patch implements the hot rotating CGM (circumgalactic medium) inflow
simulation setup based on **Stern et al. (2024) MNRAS 530, 1711**.

## Physics

The simulation models:
- Hot (~10⁶ K) rotating CGM in quasi-hydrostatic equilibrium
- Isothermal gravitational potential (v_c = constant)
- Radiative cooling using RAMSES built-in cooling tables
- Angular momentum conservation leading to spin-up during inflow
- Non-reflecting boundary conditions for subsonic acoustic waves

## Key Features

### Initial Conditions (`condinit.f90`)
- Density profile: n_H ∝ r^(-1.5) (Eq. 10 in Stern+24)
- Temperature: T ~ 2×10⁶ K (virial temperature)
- Rotation: v_φ = v_c × R_c,max / r × sin(θ) for angular momentum conservation
- Radial inflow: v_r = -r / t_cool (subsonic)

### Gravity (`gravana.f90`)
- Isothermal halo: Φ = v_c² × ln(r)
- Acceleration: g = -v_c² / r × r̂
- Use `gravity_type=4` in namelist

### Boundary Conditions (`hydro_boundary.f90`)
- **Non-reflecting BC (type 31-36)**: Allows perturbations about the
  analytical profile to exit the domain without reflection
- Key for highly subsonic flows where acoustic waves travel much faster
  than the inflow velocity
- Implemented as: u_boundary = u_analytical + (u_interior - u_analytical_interior)

## Files

```
patch/hot_cgm_inflow/
├── cgm_commons.f90      # CGM parameters and module
├── condinit.f90         # Initial conditions
├── gravana.f90          # Isothermal gravity
├── boundana.f90         # Analytical boundary profiles
├── hydro_boundary.f90   # Non-reflecting boundary conditions
├── hot_cgm.nml          # Example namelist
└── README.md            # This file
```

## Compilation

```bash
cd ramses/bin
make clean
make PATCH=../patch/hot_cgm_inflow NDIM=3 SOLVER=hydro
```

## Running

```bash
./ramses3d ../patch/hot_cgm_inflow/hot_cgm.nml
```

## Key Parameters (CGM_PARAMS namelist)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cgm_v_c` | 200 km/s | Circular velocity |
| `cgm_Rc_max` | 15 kpc | Circularization radius |
| `cgm_n_H_0` | 8×10⁻⁴ cm⁻³ | Density at r=10 kpc |
| `cgm_Lambda` | 3×10⁻²³ erg cm³/s | Cooling function |
| `cgm_mu` | 0.6 | Mean molecular weight |
| `cgm_Z` | 0.3 Z_sun | Metallicity |

## Boundary Types

Use these in `bound_type` array:
- `31`: Non-reflecting, low X boundary
- `32`: Non-reflecting, high X boundary
- `33`: Non-reflecting, low Y boundary
- `34`: Non-reflecting, high Y boundary
- `35`: Non-reflecting, low Z boundary
- `36`: Non-reflecting, high Z boundary

## Expected Results

- Hot gas inflows while remaining at ~T_vir
- Rotation increases as r decreases (angular momentum conservation)
- Gas cools and joins disc at R ~ R_c,max
- Total rotation before accretion: ~t_cool/t_ff radians (~12 rad for MW)

## References

Stern J., Fielding D., Hafen Z., et al. (2024), MNRAS, 530, 1711
"Accretion onto disc galaxies via hot and rotating CGM inflows"
https://doi.org/10.1093/mnras/stae824
