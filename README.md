# Parametric Open Bin Cardboard Box Generator

A self-contained Python CLI tool to generate SVG laser cutter files for square open bin boxes.

## Features
- **Parametric Generation**: Customize dimensions, thickness, and tab sizes.
- **Thickness Compensation**: Automatically adjusts wall widths to ensure inner dimensions match requirements ("Simple Mode").
- **Glue Tabs**: Generates tabs on North/South walls for easy assembly.
- **Fractal Engraving**: Option to engrave a Sierpinski triangle pattern.
- **Auxiliary Parts**: Can generate a matching lid and divider insert.
- **Rectangle Generator**: Includes a specific mode for generating cut-ready rectangles.

## Requirements
- Python 3.6+
- No external libraries required (uses standard library).

## Usage

Run the script from the terminal:

```bash
python3 box_generator.py
```

### Modes
1.  **Generate Rectangle SVG**: Quick generation of a simple rectangle (Task 1).
2.  **Generate Box SVGs**: The main parametric mode. You will be prompted for:
    -   Stock size (mm)
    -   Material thickness ($t$)
    -   Box Side ($S$) and Height ($H$)
    -   Tab dimensions
    -   Optional logic: Text engraving, Logo, Fractal, Lid, Divider.
3.  **Generate Mailer Box SVG**: (New) Generates a detailed FEFCO 0427 "Rollover" Mailer Box. Features double side walls (rollover) with locking tabs that fit into base slots, and a hinged lid with dust flaps.
4.  **Generate Tray SVG**: (New) Generates a lidless version of the Rollover Mailer (FEFCO 0422 style). Symmetric design with double side walls and corner flaps.
5.  **Run Example Suite**: Generates test files including the complex mailer and tray.

## Output Details
-   **File Format**: SVG (Scalable Vector Graphics), mm units.
-   **Colors**:
    -   **Red** (`rgb(255,0,0)`): Cut lines.
    -   **Blue** (`rgb(0,0,255)`): Engrave/Score lines.
-   **Stroke Width**: 0.2mm.

## Assembly
1.  **Cut** the Red lines.
2.  **Score/Engrave** the Blue lines.
3.  **Fold** along the engraved lines at the base and tabs.
4.  **Glue** the tabs on the North/South walls to the inner face of the East/West walls.
