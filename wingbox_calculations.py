"""
Aircraft Wingbox Design Calculations
=====================================

This module provides basic calculations for an aircraft wingbox design,
including structural dimensions, loads, stresses, and material properties.

Author: B02 Team
"""

import math


class WingboxDesign:
    """
    A class to perform basic wingbox design calculations for an aircraft wing.
    
    The wingbox is the primary structural component of an aircraft wing,
    consisting of spars, ribs, and skin panels that form a closed box structure.
    """
    
    def __init__(self, span, chord_root, chord_tip, thickness=0.002):
        """
        Initialize wingbox parameters.
        
        Parameters:
        -----------
        span : float
            Wing span (half-span) in meters
        chord_root : float
            Root chord length in meters
        chord_tip : float
            Tip chord length in meters
        thickness : float
            Skin thickness in meters (default: 2mm)
        """
        # Input validation
        if span <= 0:
            raise ValueError("Span must be positive")
        if chord_root <= 0:
            raise ValueError("Root chord must be positive")
        if chord_tip < 0:
            raise ValueError("Tip chord must be non-negative")
        if thickness <= 0:
            raise ValueError("Thickness must be positive")
        if chord_tip > chord_root:
            raise ValueError("Tip chord should not exceed root chord for typical designs")
            
        self.span = span
        self.chord_root = chord_root
        self.chord_tip = chord_tip
        self.thickness = thickness
        
        # Material properties (Aluminum 2024-T3)
        self.E = 73.1e9  # Young's modulus (Pa)
        self.rho = 2780  # Density (kg/m^3)
        self.sigma_yield = 345e6  # Yield strength (Pa)
        
    def calculate_taper_ratio(self):
        """
        Calculate the wing taper ratio.
        
        Returns:
        --------
        float : Taper ratio (chord_tip / chord_root)
        """
        return self.chord_tip / self.chord_root
    
    def calculate_wing_area(self):
        """
        Calculate the wing planform area (for half-wing).
        
        Returns:
        --------
        float : Wing area in m^2
        """
        return self.span * (self.chord_root + self.chord_tip) / 2
    
    def calculate_mean_aerodynamic_chord(self):
        """
        Calculate the mean aerodynamic chord (MAC).
        
        Returns:
        --------
        float : Mean aerodynamic chord in meters
        """
        taper = self.calculate_taper_ratio()
        mac = (2/3) * self.chord_root * (1 + taper + taper**2) / (1 + taper)
        return mac
    
    def calculate_aspect_ratio(self):
        """
        Calculate the wing aspect ratio.
        
        Returns:
        --------
        float : Aspect ratio (dimensionless)
        """
        area = self.calculate_wing_area()
        # For half-span, total span is 2 * span
        return (2 * self.span)**2 / (2 * area)
    
    def calculate_lift_force(self, weight, load_factor=1.0):
        """
        Calculate the lift force on the wing.
        
        Parameters:
        -----------
        weight : float
            Aircraft weight in Newtons
        load_factor : float
            Load factor (n) for maneuver (default: 1.0 for level flight)
            
        Returns:
        --------
        float : Lift force in Newtons
        """
        return weight * load_factor / 2  # Divided by 2 for half-wing
    
    def calculate_bending_moment(self, lift_force):
        """
        Calculate the maximum bending moment at the wing root.
        Assumes elliptical lift distribution.
        
        Parameters:
        -----------
        lift_force : float
            Total lift force on half-wing in Newtons
            
        Returns:
        --------
        float : Root bending moment in N·m
        """
        # For elliptical distribution, moment arm is approximately 4b/(3π)
        moment_arm = (4 * self.span) / (3 * math.pi)
        return lift_force * moment_arm
    
    def calculate_wingbox_height(self):
        """
        Estimate wingbox height based on airfoil thickness.
        Typical wingbox extends from 15% to 65% chord.
            
        Returns:
        --------
        float : Wingbox height in meters (at root)
        """
        # Assuming 12% thick airfoil (typical for transport aircraft)
        thickness_ratio = 0.12
        return self.chord_root * thickness_ratio
    
    def calculate_bending_stress(self, bending_moment):
        """
        Calculate the bending stress in the wingbox.
        Simplified calculation using beam theory.
        
        Parameters:
        -----------
        bending_moment : float
            Bending moment in N·m
            
        Returns:
        --------
        dict : Dictionary containing stress values
            - 'max_stress': Maximum bending stress (Pa)
            - 'margin_of_safety': Margin of safety (dimensionless)
        """
        if bending_moment < 0:
            raise ValueError("Bending moment must be non-negative")
            
        # Simplified wingbox cross-section
        h = self.calculate_wingbox_height()
        b = 0.5 * self.chord_root  # Width approximation
        
        # Moment of inertia for rectangular box (simplified)
        I = (b * h**3) / 12
        
        # Maximum stress at outer fiber
        max_stress = bending_moment * (h/2) / I
        
        # Margin of safety (handle zero stress case)
        if max_stress > 0:
            margin_of_safety = (self.sigma_yield / max_stress) - 1
        else:
            margin_of_safety = float('inf')  # Infinite margin if no stress
        
        return {
            'max_stress': max_stress,
            'margin_of_safety': margin_of_safety
        }
    
    def calculate_weight(self):
        """
        Estimate the structural weight of the wingbox.
        Simplified calculation based on skin area and thickness.
        
        Returns:
        --------
        float : Wingbox weight in kg
        """
        # Surface area (top and bottom skins + front and rear spars)
        area = self.calculate_wing_area()
        h = self.calculate_wingbox_height()
        
        # Approximate total skin area
        total_area = 2 * area + 2 * h * self.span
        
        # Weight = volume * density
        volume = total_area * self.thickness
        weight = volume * self.rho
        
        return weight
    
    def print_summary(self, aircraft_weight=100000, load_factor=2.5):
        """
        Print a comprehensive summary of wingbox calculations.
        
        Parameters:
        -----------
        aircraft_weight : float
            Total aircraft weight in Newtons (default: 100 kN)
        load_factor : float
            Design load factor (default: 2.5 for transport aircraft)
        """
        print("=" * 60)
        print("WINGBOX DESIGN SUMMARY")
        print("=" * 60)
        print(f"\nGeometry:")
        print(f"  Span (half-wing):        {self.span:.2f} m")
        print(f"  Root chord:              {self.chord_root:.2f} m")
        print(f"  Tip chord:               {self.chord_tip:.2f} m")
        print(f"  Taper ratio:             {self.calculate_taper_ratio():.3f}")
        print(f"  Wing area (half):        {self.calculate_wing_area():.2f} m²")
        print(f"  Mean aerodynamic chord:  {self.calculate_mean_aerodynamic_chord():.2f} m")
        print(f"  Aspect ratio:            {self.calculate_aspect_ratio():.2f}")
        print(f"  Wingbox height:          {self.calculate_wingbox_height():.3f} m")
        
        print(f"\nMaterial Properties (Aluminum 2024-T3):")
        print(f"  Young's modulus:         {self.E/1e9:.1f} GPa")
        print(f"  Density:                 {self.rho:.0f} kg/m³")
        print(f"  Yield strength:          {self.sigma_yield/1e6:.0f} MPa")
        print(f"  Skin thickness:          {self.thickness*1000:.1f} mm")
        
        lift = self.calculate_lift_force(aircraft_weight, load_factor)
        moment = self.calculate_bending_moment(lift)
        stress_results = self.calculate_bending_stress(moment)
        weight = self.calculate_weight()
        
        print(f"\nLoad Analysis (Load factor n={load_factor}):")
        print(f"  Aircraft weight:         {aircraft_weight/1000:.1f} kN")
        print(f"  Lift force (half-wing):  {lift/1000:.1f} kN")
        print(f"  Root bending moment:     {moment/1000:.1f} kN·m")
        
        print(f"\nStress Analysis:")
        print(f"  Maximum bending stress:  {stress_results['max_stress']/1e6:.1f} MPa")
        print(f"  Margin of safety:        {stress_results['margin_of_safety']:.2f}")
        
        if stress_results['margin_of_safety'] > 0:
            print(f"  Status:                  ✓ SAFE")
        else:
            print(f"  Status:                  ✗ OVERSTRESSED")
        
        print(f"\nWeight Estimation:")
        print(f"  Wingbox weight:          {weight:.1f} kg")
        print(f"  Weight per area:         {weight/self.calculate_wing_area():.1f} kg/m²")
        
        print("=" * 60)


def main():
    """
    Main function demonstrating basic wingbox calculations.
    """
    print("\nAircraft Wingbox Design Calculator")
    print("B02 Student Project\n")
    
    # Example: Regional transport aircraft wingbox
    # Typical dimensions for a small-to-medium transport aircraft
    wingbox = WingboxDesign(
        span=15.0,          # 15 meter half-span (30m total wingspan)
        chord_root=4.0,     # 4 meter root chord
        chord_tip=1.5,      # 1.5 meter tip chord
        thickness=0.002     # 2mm skin thickness
    )
    
    # Print comprehensive summary
    # Assuming 100 kN aircraft weight and 2.5 load factor
    wingbox.print_summary(aircraft_weight=100000, load_factor=2.5)
    
    print("\n✓ Calculations completed successfully!\n")


if __name__ == "__main__":
    main()
