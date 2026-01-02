# gearGeneration
This repository will generate gears and other 3d models in stl format using C++

# 🔧 Parametric Gear Generator

A powerful C++ library for generating high-quality 3D gear models with accurate involute tooth profiles. Export directly to STL format for 3D printing, CNC machining, or mechanical design visualization.

<img width="785" height="600" alt="image" src="https://github.com/user-attachments/assets/c43fc6f8-8446-43e4-a74c-24c341c02877" />
<img width="785" height="600" alt="image" src="https://github.com/user-attachments/assets/a42f024d-14c0-4c2c-b2a6-bed44582843e" />
<img width="785" height="600" alt="image" src="https://github.com/user-attachments/assets/e27e09c1-4b27-4323-82f0-003bae9e041b" />




## ✨ Features

- **🦷 True Involute Tooth Profiles** - Mathematically accurate gear teeth based on engineering standards
- **📐 Multiple Pressure Angles** - Support for 14.5°, 20°, 25° and custom angles
- **⚙️ Four Gear Types** - Spur, Helical, Bevel, and Worm gears
- **🎨 Flexible Tooth Styles** - Involute, rounded, or simple rectangular teeth
- **🕳️ Center Hole Support** - Optional shaft holes with configurable diameter
- **📤 STL Export** - Direct export to industry-standard STL format
- **🎯 Fully Parametric** - Control every aspect of gear geometry

## 🎯 Gear Types Supported

### 1. **Spur Gears** 🔄
The most common type - flat gears with straight teeth parallel to the axis.
- **Use case:** Parallel shaft transmission
- **Advantages:** Simple, efficient, cost-effective
- **Applications:** Clocks, washing machines, mechanical toys

### 2. **Helical Gears** 🌀
Gears with angled teeth that form a helix around the gear body.
- **Use case:** High-speed parallel shafts
- **Advantages:** Smoother, quieter operation; higher load capacity
- **Applications:** Car transmissions, industrial machinery

### 3. **Bevel Gears** 📐
Cone-shaped gears for transmitting motion between intersecting shafts.
- **Use case:** 90° angle transmission
- **Advantages:** Changes rotation axis direction
- **Applications:** Hand drills, differential gears, printing presses


## 📊 Pressure Angles Explained

The pressure angle determines how gear teeth interact and affects strength, smoothness, and bearing loads.

| Pressure Angle | Smoothness | Tooth Strength | Bearing Load | Best For |
|----------------|------------|----------------|--------------|----------|
| **14.5°** | ⭐⭐⭐⭐⭐ Very Smooth | ⭐⭐ Low | ⭐⭐ Low | Light-duty, legacy systems |
| **20°** | ⭐⭐⭐⭐ Smooth | ⭐⭐⭐⭐ High | ⭐⭐⭐ Medium | **Standard choice** - Most applications |
| **25°** | ⭐⭐⭐ Less Smooth | ⭐⭐⭐⭐⭐ Very High | ⭐⭐⭐⭐ High | Heavy-duty industrial systems |

**Visual Comparison:**
- **14.5°** → Wider, thicker teeth (gentle curves)
- **20°** → Balanced medium teeth (industry standard)
- **25°** → Narrower, pointed teeth (stronger, more force)

## 🚀 Quick Start

### Prerequisites
- C++ compiler with C++11 support (g++, clang++, MSVC)
- Standard C++ library (no external dependencies!)

### Compilation

```bash
# Linux/macOS
g++ -std=c++11 -O2 main.cpp -o main

# Windows (MSVC)
cl /std:c++11 /O2 main.cpp

# With warnings
g++ -std=c++11 -O2 -Wall -Wextra main.cpp -o main
```

### Basic Usage

```bash
# Run the generator
./main

# Output files will be created in the current directory
# gear_pressure_20.stl, gear_helical.stl, etc.
```

## 📐 Parameter Reference

### GearParams Structure

```cpp
struct GearParams {
    int numTeeth;           // Number of teeth (minimum: 3)
    float innerRadius;      // Root radius (valley depth)
    float pitchRadius;      // Working radius where gears mesh
    float outerRadius;      // Tip radius (tooth height)
    float thickness;        // Gear thickness/width
    float toothWidth;       // Tooth width ratio (0.0-1.0)
    int segments;           // Smoothness (higher = smoother)
    GearType type;          // SPUR, HELICAL, BEVEL, WORM
    float helixAngle;       // Twist angle for helical gears (degrees)
    float pressureAngle;    // 14.5°, 20°, or 25° typical
    bool useInvolute;       // true = engineering profile
    bool roundedTeeth;      // true = simplified rounded teeth
    float roundingRadius;   // Rounding amount for rounded teeth
    float holeRadius;       // Center hole size (0 = solid)
};
```

### Parameter Guidelines

#### **numTeeth**
- **Minimum:** 3 teeth
- **Typical:** 12-80 teeth
- **Rule:** Fewer teeth = larger/chunkier, more teeth = smaller/finer

#### **Radius Values** (innerRadius < pitchRadius < outerRadius)
- **innerRadius:** Root circle - typically 70-80% of pitch radius
- **pitchRadius:** The "working" diameter - main dimension
- **outerRadius:** Tip circle - typically 110-115% of pitch radius

#### **toothWidth**
- **Range:** 0.3 to 0.5
- **Standard:** 0.45 (45% tooth, 55% gap)
- **Higher:** Stronger teeth but may cause interference

#### **segments**
- **Low (2-3):** Fast generation, blocky appearance
- **Medium (4-6):** Good balance (recommended)
- **High (8+):** Very smooth but slower, larger files

#### **pressureAngle**
- **14.5°:** Smooth, quiet, weak (older standard)
- **20°:** Standard modern gears ⭐ **RECOMMENDED**
- **25°:** Strong, handles high loads

## 💻 Code Examples

### Example 1: Standard 20° Spur Gear
```cpp
GearParams standardGear = {
    20,      // 20 teeth
    35.0f,   // Inner radius
    45.0f,   // Pitch radius
    50.0f,   // Outer radius
    10.0f,   // Thickness
    0.45f,   // Tooth width
    4,       // Segments
    SPUR,    // Type
    0,       // No helix angle
    20.0f,   // 20° pressure angle
    true,    // Use involute
    false,   // Not rounded
    2.0f,    // Rounding (unused)
    8.0f     // Center hole
};

std::ofstream file("my_gear.stl");
file << "solid my_gear\n";
generateGear(file, standardGear);
file << "endsolid my_gear\n";
file.close();
```

### Example 2: Helical Gear for Smooth Operation
```cpp
GearParams helicalGear = {
    24,      // 24 teeth
    38.0f,   // Inner radius
    48.0f,   // Pitch radius
    53.0f,   // Outer radius
    20.0f,   // Thicker for helix
    0.45f,   // Tooth width
    4,       // Segments
    HELICAL, // Helical type
    30.0f,   // 30° helix angle
    20.0f,   // 20° pressure angle
    true,    // Use involute
    false,   // Not rounded
    2.0f,    // Rounding (unused)
    12.0f    // Larger center hole
};
```

### Example 3: High-Strength 25° Gear
```cpp
GearParams strongGear = {
    18,      // 18 teeth
    40.0f,   // Inner radius
    50.0f,   // Pitch radius
    55.0f,   // Outer radius
    15.0f,   // Thicker body
    0.45f,   // Tooth width
    5,       // Smoother
    SPUR,    // Type
    0,       // No helix
    25.0f,   // 25° for strength!
    true,    // Use involute
    false,   // Not rounded
    2.0f,    // Rounding (unused)
    10.0f    // Center hole
};
```

### Example 4: Worm Gear for Speed Reduction
```cpp
GearParams wormGear = {
    3,       // 3 thread starts
    8.0f,    // Inner radius
    12.0f,   // Pitch radius
    15.0f,   // Outer radius
    80.0f,   // Long body
    0.3f,    // Thread width
    8,       // High smoothness
    WORM,    // Worm type
    0,       // (unused for worm)
    20.0f,   // Pressure angle
    false,   // Special worm profile
    false,   // Not rounded
    2.0f,    // Rounding (unused)
    5.0f     // Center hole
};
```

## 🎨 Customization Tips

### Creating Matching Gear Pairs

For two gears to mesh properly:

1. **Same pressure angle** - Both must use identical pressure angles
2. **Same pitch radius ratio** - Gear ratio = teeth ratio
3. **Complementary helix angles** - One positive, one negative (helical)

```cpp
// Gear 1: 20 teeth, 40mm pitch radius
GearParams gear1 = {20, 30.0f, 40.0f, 45.0f, ...};

// Gear 2: 40 teeth, 80mm pitch radius (2:1 ratio)
GearParams gear2 = {40, 70.0f, 80.0f, 85.0f, ...};
```

### Optimizing for 3D Printing

```cpp
GearParams printable = {
    numTeeth: 16,         // Not too many
    segments: 4,          // Medium quality
    thickness: 8.0f,      // At least 3-5mm
    holeRadius: 4.0f,     // For 8mm shaft with tolerance
    useInvolute: true,    // Better meshing
    pressureAngle: 20.0f  // Standard
};
```

### Creating Display/Decorative Gears

```cpp
GearParams decorative = {
    numTeeth: 12,
    segments: 6,          // Higher for visual quality
    roundedTeeth: true,   // Smoother look
    holeRadius: 0,        // Solid for mounting
    thickness: 5.0f       // Thinner for wall mounting
};
```

## 📁 Output Format

The library generates **ASCII STL files** with the following structure:

```
solid gear_name
    facet normal nx ny nz
        outer loop
            vertex x1 y1 z1
            vertex x2 y2 z2
            vertex x3 y3 z3
        endloop
    endfacet
    ...
endsolid gear_name
```

**File sizes:**
- Simple gear (20 teeth, 4 segments): ~500 KB
- Complex gear (40 teeth, 8 segments): ~2-3 MB
- Worm gear (long body): ~1-2 MB

## 🔧 Advanced Features

### Involute Curve Mathematics

The library implements true involute curves using the involute function:

```
inv(α) = tan(α) - α
```

Where `α` is the pressure angle. This ensures:
- ✅ Constant velocity ratio during meshing
- ✅ Smooth power transmission
- ✅ Industry-standard tooth profiles

### Tooth Narrowing Algorithm

Teeth automatically narrow from root to tip based on:

```cpp
toothAngleAtRadius = toothAngleAtPitch × (pitchRadius / currentRadius)
```

This creates the characteristic curved involute profile where higher pressure angles produce more pronounced narrowing.

## 🛠️ Troubleshooting

### Issue: Gears look blocky
**Solution:** Increase `segments` parameter (try 5-8)

### Issue: File size too large
**Solution:** Decrease `segments` parameter or reduce `numTeeth`

### Issue: Teeth look wrong/split
**Solution:** Ensure `innerRadius < pitchRadius < outerRadius` and `useInvolute = true`

### Issue: Gears won't mesh
**Solution:** 
- Use same pressure angle for both gears
- Check pitch radius ratio matches teeth ratio
- Ensure proper center distance

### Issue: Compilation errors
**Solution:**
```bash
# Add M_PI definition if needed
g++ -D_USE_MATH_DEFINES main.cpp
```

## 📚 Technical Background

### Involute Gear Theory

Involute gears are the industry standard because:
1. **Constant velocity ratio** - No speed fluctuation
2. **Center distance tolerance** - Works with slight variations
3. **Easy manufacturing** - Can be cut with standard tools
4. **Efficient power transmission** - Minimal friction losses

### Pressure Angle Effects

| Aspect | 14.5° | 20° | 25° |
|--------|-------|-----|-----|
| Tooth thickness | Thickest | Medium | Thinnest |
| Contact ratio | Highest | Medium | Lowest |
| Radial force | Lowest | Medium | Highest |
| Bending strength | Lowest | Medium | Highest |

## 🎓 Learning Resources

- [Gear Theory Basics](https://www.engineeringtoolbox.com/gears-d_1856.html)
- [Involute Curve Mathematics](https://en.wikipedia.org/wiki/Involute_gear)
- [AGMA Standards](https://www.agma.org/) (American Gear Manufacturers Association)
- [3D Printing Gears Guide](https://www.prusaprinters.org/prints/category/mechanical)

## 🤝 Contributing

Contributions are welcome! Areas for improvement:

- [ ] Binary STL output for smaller files
- [ ] Gear pair meshing validation
- [ ] Internal/ring gear support
- [ ] Rack gear generation
- [ ] Epicyclic/planetary gear systems
- [ ] Automatic tolerance calculation
- [ ] GUI parameter interface
- [ ] Batch generation tools

## 📄 License

This project is open source. Feel free to use, modify, and distribute for both personal and commercial projects.

## 🙏 Acknowledgments

Based on classical involute gear theory and mechanical engineering principles. Inspired by the need for accessible, high-quality parametric gear generation tools.

## 📞 Support

- **Issues:** Open an issue on GitHub
- **Questions:** Check the examples and parameter reference
- **Feature Requests:** Submit via GitHub issues

---

**Made with ⚙️ by engineers, for engineers**

*Generate, print, and build amazing mechanical systems!*
