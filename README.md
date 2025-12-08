# B02 - Aircraft Wingbox Design
The greatest thing that never flew

## Project Description

This repository contains the design calculations for an aircraft wingbox structure. The wingbox is the primary load-bearing structure in an aircraft wing, consisting of spars, ribs, and skin panels that form a closed box section.

## Features

The main calculation module (`wingbox_calculations.py`) provides:

- **Geometric Analysis**: Wing dimensions, taper ratio, area, mean aerodynamic chord, aspect ratio
- **Structural Calculations**: Wingbox height estimation based on airfoil thickness
- **Load Analysis**: Lift force and bending moment calculations
- **Stress Analysis**: Bending stress calculations with margin of safety
- **Weight Estimation**: Structural weight based on material and geometry
- **Material Properties**: Aluminum 2024-T3 specifications

## Usage

### Running the Calculator

Simply execute the main script:

```bash
python3 wingbox_calculations.py
```

### Using as a Module

You can also import the `WingboxDesign` class in your own scripts:

```python
from wingbox_calculations import WingboxDesign

# Create a wingbox design
wingbox = WingboxDesign(
    span=15.0,          # 15m half-span
    chord_root=4.0,     # 4m root chord
    chord_tip=1.5,      # 1.5m tip chord
    thickness=0.002     # 2mm skin thickness
)

# Calculate wing properties
area = wingbox.calculate_wing_area()
taper = wingbox.calculate_taper_ratio()
mac = wingbox.calculate_mean_aerodynamic_chord()

# Perform structural analysis
lift = wingbox.calculate_lift_force(weight=100000, load_factor=2.5)
moment = wingbox.calculate_bending_moment(lift)
stress_results = wingbox.calculate_bending_stress(moment)

# Print comprehensive summary
wingbox.print_summary(aircraft_weight=100000, load_factor=2.5)
```

## Requirements

- Python 3.6 or higher
- No external dependencies (uses only standard library)

## Calculations

### Geometry
- **Taper Ratio**: λ = c_tip / c_root
- **Wing Area**: S = b × (c_root + c_tip) / 2
- **Mean Aerodynamic Chord**: MAC = (2/3) × c_root × (1 + λ + λ²) / (1 + λ)
- **Aspect Ratio**: AR = b² / S

### Loads
- **Lift Force**: L = W × n / 2 (for half-wing)
- **Bending Moment**: M = L × (4b / 3π) (elliptical distribution)

### Stress
- **Bending Stress**: σ = M × y / I
- **Margin of Safety**: MS = (σ_yield / σ) - 1

### Weight
- Estimated based on skin area, thickness, and material density

## Example Output

```
============================================================
WINGBOX DESIGN SUMMARY
============================================================

Geometry:
  Span (half-wing):        15.00 m
  Root chord:              4.00 m
  Tip chord:               1.50 m
  Taper ratio:             0.375
  Wing area (half):        41.25 m²
  Mean aerodynamic chord:  2.94 m
  Aspect ratio:            10.91
  Wingbox height:          0.480 m

Load Analysis (Load factor n=2.5):
  Aircraft weight:         100.0 kN
  Lift force (half-wing):  125.0 kN
  Root bending moment:     795.8 kN·m

Stress Analysis:
  Maximum bending stress:  10.4 MPa
  Margin of safety:        32.30
  Status:                  ✓ SAFE

Weight Estimation:
  Wingbox weight:          538.8 kg
  Weight per area:         13.1 kg/m²
============================================================
```

## Project Team

B02 Student Project Team

## License

This is a student project for educational purposes.
