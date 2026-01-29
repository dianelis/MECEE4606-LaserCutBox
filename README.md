# Parametric Acrylic Box Generator

A self-contained Python CLI tool to generate SVG laser cutter files for **Acrylic (Finger Joint)** boxes.

## Features
- **Task 1: Small Rectangle**: Strict check for <= 2x2 inch rectangle generation.
- **Task 2: Acrylic Box (Finger Joints)**: Interactive generation of a rigid box with:
    - 5 Interlocking Parts (Base, Front, Back, Left, Right).
    - **Finger Joints**: Standard rectangular tabs for acrylic assembly.
    - **Features**:
        - **Divider**: Removable insert with automatic vertical supporting slots on Front/Back walls.
        - **Engraving**:
            - Base Text: Defaults to "Digital Manufacturing".
            - Columbia Logo: Optional on Front Wall.
- **Automatic Formatting**:
    - **Red**: Cut lines (0.01mm stroke).
    - **Blue**: Engrave lines (0.01mm stroke).
- **Software Description**: Prints calculations and formulas after generation.

## Requirements
- Python 3.6+
- No external libraries required.

## Usage

Run the script from the terminal:

```bash
python3 box_generator.py
```

### Modes
1.  **Task 1: Small Rectangle**: Quick generation of a simple rectangle (max 50.8mm).
2.  **Task 2: Acrylic Box**: The main parametric mode. You will be prompted for:
    -   Stock size (mm)
    -   Material thickness ($t$)
    -   Box Dimensions ($S, H$) [Inner Size]
    -   Features (Divider, Text, Logo, Fractal)

## Assembly (Acrylic)
1.  **Cut** the Red lines.
2.  **Engrave** the Blue lines.
3.  **Assemble** using acrylic cement:
    -   The Base has tabs on all sides.
    -   Front/Back walls mate with the Base and overlap the Side walls.
    -   Use the finger joints for alignment.
