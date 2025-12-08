"""
Example usage of the wingbox design calculator.

This script demonstrates how to use the WingboxDesign class
to analyze different wing configurations.
"""

from wingbox_calculations import WingboxDesign


def example_small_aircraft():
    """Example: Small general aviation aircraft"""
    print("\n" + "=" * 60)
    print("EXAMPLE 1: Small General Aviation Aircraft")
    print("=" * 60)
    
    wingbox = WingboxDesign(
        span=6.0,           # 6m half-span (12m total)
        chord_root=1.8,     # 1.8m root chord
        chord_tip=0.9,      # 0.9m tip chord
        thickness=0.0015    # 1.5mm skin
    )
    
    wingbox.print_summary(aircraft_weight=15000, load_factor=3.8)


def example_regional_transport():
    """Example: Regional transport aircraft"""
    print("\n" + "=" * 60)
    print("EXAMPLE 2: Regional Transport Aircraft")
    print("=" * 60)
    
    wingbox = WingboxDesign(
        span=15.0,          # 15m half-span (30m total)
        chord_root=4.0,     # 4m root chord
        chord_tip=1.5,      # 1.5m tip chord
        thickness=0.002     # 2mm skin
    )
    
    wingbox.print_summary(aircraft_weight=100000, load_factor=2.5)


def example_large_airliner():
    """Example: Large commercial airliner"""
    print("\n" + "=" * 60)
    print("EXAMPLE 3: Large Commercial Airliner")
    print("=" * 60)
    
    wingbox = WingboxDesign(
        span=30.0,          # 30m half-span (60m total)
        chord_root=8.0,     # 8m root chord
        chord_tip=2.0,      # 2m tip chord
        thickness=0.003     # 3mm skin
    )
    
    wingbox.print_summary(aircraft_weight=600000, load_factor=2.5)


def compare_configurations():
    """Compare different wing configurations"""
    print("\n" + "=" * 60)
    print("CONFIGURATION COMPARISON")
    print("=" * 60)
    
    configs = [
        ("Low aspect ratio", 10.0, 5.0, 3.0),
        ("Medium aspect ratio", 15.0, 4.0, 1.5),
        ("High aspect ratio", 20.0, 3.0, 1.0),
    ]
    
    print(f"\n{'Configuration':<25} {'AR':<8} {'Area (m²)':<12} {'Weight (kg)':<12}")
    print("-" * 60)
    
    for name, span, c_root, c_tip in configs:
        wb = WingboxDesign(span, c_root, c_tip)
        ar = wb.calculate_aspect_ratio()
        area = wb.calculate_wing_area()
        weight = wb.calculate_weight()
        print(f"{name:<25} {ar:<8.2f} {area:<12.2f} {weight:<12.1f}")


if __name__ == "__main__":
    print("\nAircraft Wingbox Design Examples")
    print("B02 Student Project")
    
    example_small_aircraft()
    example_regional_transport()
    example_large_airliner()
    compare_configurations()
    
    print("\n✓ All examples completed successfully!\n")
