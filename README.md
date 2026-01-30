# Parametric Acrylic Box Generator

A Python tool to generate **fabrication-ready SVG files** for laser-cut acrylic boxes.

## 🚀 Features

### Core Design
- **True Rectangular Dimensions**: Customize **Width**, **Depth**, and **Height** independently (e.g., Shoe Box, Tower, Cube).
- **Finger Joints**: Dynamic tab generation for rigid, self-aligning assembly.
- **Glue-Only Assembly**: sleek, screw-free design for a clean aesthetic.

### Components
- **Removable Lid**: Friction-fit lid (step design) that lifts off easily.
- **Divider Insert**: Optional internal divider that slides into slots.
- **Engravings**:
  - **Sierpinski Triangle Fractal**: Generative art on side walls.
  - **Text Labels**: Customizable text on Base and Side walls.
  - **Columbia Logo**: Optional insertion.

## 📂 File Structure

```text
.
├── box_generator.py          # Main CLI Application
├── README.md                 # Documentation
├── output/                   # Generated .svg files go here
│   └── box_acrylic_parts.svg
├── generate_scripts/         # Helper scripts for preset sizes
│   ├── generate_15cm.py      # 15cm Cube
│   ├── generate_20x10.py     # 20cm x 10cm Rectangular
│   ├── generate_15x10.py     # 15cm x 10cm Rectangular
│   └── generate_10x20.py     # 10cm x 20cm Tall Tower
└── logo/                     # Asset directory
```

## 🛠 Usage

### 1. Interactive Mode
Run the main script and follow the prompts:

```bash
python3 box_generator.py
```

**Prompts:**
1.  **Stock Size**: e.g., `600` (mm)
2.  **Thickness**: e.g., `3.0` (mm)
3.  **Dimensions**:
    -   `Width (Front/Back) W`: Length of the box front.
    -   `Depth (Left/Right) D`: Length of the box side.
    -   `Height H`: Vertical height.
4.  **Features**: Toggle Divider, Lid, Engravings, etc.

### 2. Preset Scripts
Quickly generate standard sizes using the scripts in `generate_scripts/`:

```bash
# Generate a 20cm x 10cm x 10cm box
python3 generate_scripts/generate_20x10.py

# Generate a 15cm Cube
python3 generate_scripts/generate_15cm.py
```

## 🧩 Assembly Instructions

**Materials**: 3mm Acrylic Sheet, Acrylic Cement (e.g., Weld-On 4).

1.  **Cut & Engrave**: Use `output/box_acrylic_parts.svg`.
    -   **Red Lines**: CUT (0.01mm stroke)
    -   **Blue Lines**: ENGRAVE (0.01mm stroke)
2.  **Base Assembly**:
    -   Apply cement to the finger joints of the **Base**.
    -   Attach **Front**, **Back**, **Left**, and **Right** walls.
    -   *Tip: Use the Divider insert to square the walls while drying.*
3.  **Lid Assembly**:
    -   The Lid consists of a top panel and 4 "lip strips".
    -   Glue the lip strips to the **underside** of the top panel, inset by 3mm.
    -   This creates a step that fits inside the box opening.
4.  **Finish**: Slide in the divider and place the lid on top.
