#!/usr/bin/env python3
"""
Parametric Open Bin Cardboard Box Generator
Author: Antigravity (Google Deepmind)

This program generates SVG laser cutter files for a square open bin box.
It runs as a CLI and produces clean, compliant SVG files.
"""

import math
import os

# --- Constants ---

STROKE_CUT = "rgb(255,0,0)"
STROKE_ENGRAVE = "rgb(0,0,255)"
STROKE_WIDTH = "0.2"

# --- SVG Writer Class ---

class SVGGenerator:
    """A simple self-contained SVG generator."""
    def __init__(self, filename, width, height, unit="mm"):
        self.filename = os.path.join("output", filename)
        # Ensure output directory exists
        os.makedirs("output", exist_ok=True)
        self.width = width
        self.height = height
        self.unit = unit
        self.elements = []
        # Pre-create standard groups
        self.groups = {
            "CUT": [],
            "ENGRAVE": []
        }

    def add_line(self, x1, y1, x2, y2, mode="CUT"):
        """Adds a line to the specified group (CUT or ENGRAVE)."""
        line_str = f'<line x1="{x1:.4f}" y1="{y1:.4f}" x2="{x2:.4f}" y2="{y2:.4f}" />'
        if mode in self.groups:
            self.groups[mode].append(line_str)
        else:
            print(f"Warning: Unknown mode {mode}, defaulting to CUT")
            self.groups["CUT"].append(line_str)

    def add_rect(self, x, y, w, h, mode="CUT"):
        """Adds a rectangle (as lines) for compatibility."""
        self.add_line(x, y, x + w, y, mode)
        self.add_line(x + w, y, x + w, y + h, mode)
        self.add_line(x + w, y + h, x, y + h, mode)
        self.add_line(x, y + h, x, y, mode)

    def add_polyline(self, points, mode="ENGRAVE"):
        """Adds a polyline from a list of (x,y) tuples."""
        if not points:
            return
        points_str = " ".join([f"{p[0]:.4f},{p[1]:.4f}" for p in points])
        # Polyline element
        poly_str = f'<polyline points="{points_str}" />'
        if mode in self.groups:
            self.groups[mode].append(poly_str)

    def add_text(self, x, y, text, font_size=5, mode="ENGRAVE"):
        """Adds text centered at x,y."""
        # Note: SVG text centering uses text-anchor="middle"
        # We will follow spec: stroke blue, fill none.
        t_str = (f'<text x="{x:.4f}" y="{y:.4f}" font-family="sans-serif" font-size="{font_size}" '
                 f'text-anchor="middle" dominant-baseline="middle">{text}</text>')
        if mode in self.groups:
            self.groups[mode].append(t_str)

    def save(self):
        """Writes the SVG file."""
        with open(self.filename, 'w') as f:
            f.write('<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n')
            f.write(f'<svg width="{self.width}{self.unit}" height="{self.height}{self.unit}" '
                    f'viewBox="0 0 {self.width} {self.height}" '
                    f'xmlns="http://www.w3.org/2000/svg">\n')
            
            # Write Groups
            # CUT
            f.write(f'  <g id="CUT" stroke="{STROKE_CUT}" stroke_width="{STROKE_WIDTH}" fill="none">\n')
            for el in self.groups["CUT"]:
                f.write(f'    {el}\n')
            f.write('  </g>\n')

            # ENGRAVE
            f.write(f'  <g id="ENGRAVE" stroke="{STROKE_ENGRAVE}" stroke_width="{STROKE_WIDTH}" fill="none">\n')
            for el in self.groups["ENGRAVE"]:
                f.write(f'    {el}\n')
            f.write('  </g>\n')
            
            f.write('</svg>\n')
        print(f"Successfully generated: {self.filename}")


# --- Fractal Generator ---

def generate_sierpinski(x, y, size, depth=4):
    """Generates a Sierpinski triangle as a list of independent triangles (polylines).
    Returns a list of list of points (list of polylines).
    """
    lines = []

    def sierpinski(ax, ay, bx, by, cx, cy, d):
        if d == 0:
            # Add this triangle
            lines.append([(ax, ay), (bx, by), (cx, cy), (ax, ay)])
        else:
            # Midpoints
            abx, aby = (ax + bx) / 2, (ay + by) / 2
            bcx, bcy = (bx + cx) / 2, (by + cy) / 2
            cax, cay = (cx + ax) / 2, (cy + ay) / 2
            
            sierpinski(ax, ay, abx, aby, cax, cay, d - 1)
            sierpinski(abx, aby, bx, by, bcx, bcy, d - 1)
            sierpinski(cax, cay, bcx, bcy, cx, cy, d - 1)

    # Calculate vertices of the main equilateral triangle fitting in the box (x,y, size)
    # Height of equilateral triangle = size * sqrt(3)/2
    h_tri = size * math.sqrt(3) / 2
    
    # Vertices relative to (x,y) being top-left of the bounding box
    # Bottom Left
    p1x, p1y = x, y + h_tri
    # Bottom Right
    p2x, p2y = x + size, y + h_tri
    # Top Center
    p3x, p3y = x + size / 2, y

    sierpinski(p1x, p1y, p2x, p2y, p3x, p3y, depth)
    return lines


# --- Geometry Functions ---

def generate_rectangle_svg(width_mm, height_mm, filename="rectangle.svg"):
    """Generates a simple rectangle SVG."""
    # Validation
    if width_mm > 50.8 or height_mm > 50.8:
        print("Error: For rectangle task, dimensions must be <= 50.8 mm (2 inches).")
        return False
    if width_mm <= 0 or height_mm <= 0:
        print("Error: Dimensions must be positive.")
        return False

    # Add a small margin to the stock size for the SVG viewbox
    stock_w = width_mm + 10
    stock_h = height_mm + 10
    
    svg = SVGGenerator(filename, stock_w, stock_h)
    
    # Center the rectangle
    start_x = 5
    start_y = 5
    
    svg.add_rect(start_x, start_y, width_mm, height_mm, mode="CUT")
    svg.save()
    return True


def generate_box_svg(params, filename="box_square.svg"):
    """Generates the main box layout."""
    # Unpack params
    S_inner = params['S']
    H_inner = params['H']
    t = params['t']
    stock_w = params['stock_w']
    stock_h = params['stock_h']
    
    # Thickness Compensation ("Simple" Mode Strategy):
    # Base Panel: S x S (Inner dimensions).
    # North/South Walls (Front/Back): Width S. Attached to Base.
    # East/West Walls (Left/Right): Width S + 2t. Attached to Base.
    # This configuration ensures that when folded, the North/South walls 
    # sit between the East/West walls, maintaining the inner dimension S.
    
    base_w = S_inner
    base_h = S_inner
    
    # Walls
    wall_ns_w = S_inner
    wall_ns_h = H_inner
    
    wall_ew_w = S_inner + 2 * t
    wall_ew_h = H_inner
    
    # Tabs
    tab_w = params.get('tab_width', 20)
    tab_d = params.get('tab_depth', 15)
    
    # Layout Check
    # Total Width = West_H + Base_W + East_H = H + S + H
    # Total Height = North_H + Base_H + South_H = H + S + H
    
    total_layout_w = wall_ew_h + base_w + wall_ew_h
    total_layout_h = wall_ns_h + base_h + wall_ns_h
    
    if total_layout_w > stock_w or total_layout_h > stock_h:
        print(f"Error: Layout size ({total_layout_w:.1f} x {total_layout_h:.1f} mm) larger than stock ({stock_w} x {stock_h} mm).")
        return False
        
    svg = SVGGenerator(filename, stock_w, stock_h)
    
    # Center the layout
    cx = stock_w / 2
    cy = stock_h / 2
    
    # Coordinates for Base
    bx = cx - base_w / 2
    by = cy - base_h / 2
    
    # 1. Base (S x S) - Central region with folds.
    # Top Fold (Base-North)
    svg.add_line(bx, by, bx + base_w, by, mode="ENGRAVE") 
    # Bottom Fold (Base-South)
    svg.add_line(bx, by + base_h, bx + base_w, by + base_h, mode="ENGRAVE")
    
    # --- North Wall (Top) ---
    # Attached to Base Top Edge (by).
    
    # Helper to draw tab on vertical edge
    # x, y1, y2, direction (+1 for right, -1 for left)
    def draw_tab_vertical(vx, vy1, vy2, direction):
        # Center the tab
        h_edge = abs(vy2 - vy1)
        curr_tab_d = min(tab_d, h_edge) # Depth
        curr_tab_w = min(tab_w, h_edge) # Width along edge
        
        # Inset
        inset = (h_edge - curr_tab_w) / 2
        
        # Points for trapezoidal tab
        p1 = (vx, vy1)
        p2 = (vx, vy1 + inset)
        chamfer = 3
        p3 = (vx + direction * curr_tab_d, vy1 + inset + chamfer)
        p4 = (vx + direction * curr_tab_d, vy2 - inset - chamfer)
        p5 = (vx, vy2 - inset)
        p6 = (vx, vy2)
        
        # Draw the cut lines
        svg.add_line(p1[0], p1[1], p2[0], p2[1], "CUT")
        svg.add_line(p2[0], p2[1], p3[0], p3[1], "CUT")
        svg.add_line(p3[0], p3[1], p4[0], p4[1], "CUT")
        svg.add_line(p4[0], p4[1], p5[0], p5[1], "CUT")
        svg.add_line(p5[0], p5[1], p6[0], p6[1], "CUT")
        
        # Add Fold Line where tab meets wall
        svg.add_line(p2[0], p2[1], p5[0], p5[1], "ENGRAVE")

    # North Wall Draw
    # Top Edge:
    svg.add_line(bx, by - wall_ns_h, bx + wall_ns_w, by - wall_ns_h, "CUT")
    # Left Edge (with Leftward Tab)
    draw_tab_vertical(bx, by - wall_ns_h, by, -1)
    # Right Edge (with Rightward Tab)
    draw_tab_vertical(bx + wall_ns_w, by - wall_ns_h, by, 1)
    
    # --- South Wall (Bottom) ---
    # Attached to Base Bottom Edge (by + base_h).
    bs_y = by + base_h
    # Bottom Edge
    svg.add_line(bx, bs_y + wall_ns_h, bx + wall_ns_w, bs_y + wall_ns_h, "CUT")
    # Left Edge
    draw_tab_vertical(bx, bs_y, bs_y + wall_ns_h, -1)
    # Right Edge
    draw_tab_vertical(bx + wall_ns_w, bs_y, bs_y + wall_ns_h, 1)

    # --- West Wall (Left) ---
    # Attached to Base Left Edge.
    # Note: West Wall is wider (S + 2t). It extends 't' vertically past the base top/bottom
    # to cover the thickness of the North/South walls.
    
    # Fold Line
    svg.add_line(bx, by, bx, by + base_h, "ENGRAVE")
    
    ww_x_inner = bx
    ww_x_outer = bx - wall_ew_h # H wide
    
    ww_y_top = by - t
    ww_y_bot = by + base_h + t
    
    # Draw perimeter
    # Top Edge (Horizontal)
    svg.add_line(ww_x_outer, ww_y_top, ww_x_inner, ww_y_top, "CUT")
    # Left Edge (Vertical)
    svg.add_line(ww_x_outer, ww_y_top, ww_x_outer, ww_y_bot, "CUT")
    # Bottom Edge (Horizontal)
    svg.add_line(ww_x_outer, ww_y_bot, ww_x_inner, ww_y_bot, "CUT")
    
    # Inner Edge (Vertical) - Splits at the base corner (clearance cuts)
    svg.add_line(ww_x_inner, ww_y_top, ww_x_inner, by, "CUT")
    svg.add_line(ww_x_inner, by + base_h, ww_x_inner, ww_y_bot, "CUT")
    
    # --- East Wall (Right) ---
    ew_x_inner = bx + base_w
    ew_x_outer = bx + base_w + wall_ew_h
    
    ew_y_top = ww_y_top
    ew_y_bot = ww_y_bot
    
    # Fold
    svg.add_line(ew_x_inner, by, ew_x_inner, by + base_h, "ENGRAVE")
    
    # Perimeter
    svg.add_line(ew_x_inner, ew_y_top, ew_x_outer, ew_y_top, "CUT") # Top
    svg.add_line(ew_x_outer, ew_y_top, ew_x_outer, ew_y_bot, "CUT") # Right
    svg.add_line(ew_x_outer, ew_y_bot, ew_x_inner, ew_y_bot, "CUT") # Bot
    
    # Connection cuts
    svg.add_line(ew_x_inner, ew_y_top, ew_x_inner, by, "CUT")
    svg.add_line(ew_x_inner, by + base_h, ew_x_inner, ew_y_bot, "CUT")

    # --- Features ---
    
    # 1. Text Top (on Bottom Panel / Base)
    if params.get('text_top'):
        # Center on base
        svg.add_text(bx + base_w/2, by + base_h/2, params['text_top'], font_size=8, mode="ENGRAVE")
        
    # 2. Text Front (on South Wall)
    if params.get('text_front'):
        # Center on South Wall
        svg.add_text(bx + base_w/2, bs_y + wall_ns_h/2, params['text_front'], font_size=8, mode="ENGRAVE")

    # 3. Logo "Columbia Digital Manufacturing"
    if params.get('include_logo'):
        # Put it on North Wall
        svg.add_text(bx + base_w/2, by - wall_ns_h/2 - 4, "Columbia", font_size=6, mode="ENGRAVE")
        svg.add_text(bx + base_w/2, by - wall_ns_h/2 + 4, "Digital Manufacturing", font_size=4, mode="ENGRAVE")
        
    # 4. Fractal (on West Wall)
    if params.get('include_fractal'):
        # Fit inside West Wall (H x S approx)
        # West Wall center:
        cx_w = (ww_x_inner + ww_x_outer) / 2
        cy_w = (ww_y_top + ww_y_bot) / 2
        
        # Dimensions available (with margin)
        avail_w = abs(ww_x_outer - ww_x_inner) - 10
        avail_h = abs(ww_y_bot - ww_y_top) - 10
        size = min(avail_w, avail_h)
        
        # Top-left of fractal box
        # generate_sierpinski draws from (x,y) downwards.
        h_tri = size * math.sqrt(3) / 2
        fx = cx_w - size/2
        fy = cy_w - h_tri/2
        
        polys = generate_sierpinski(fx, fy, size)
        for poly in polys:
            svg.add_polyline(poly, mode="ENGRAVE")

    svg.save()
    return True

def generate_lid_svg(params, filename="lid_square.svg"):
    """Generates a simple flat lid (shallow box) to fit over the main box."""
    t = params['t']
    S = params['S']
    
    # Logically the box is a square of width S (inner) + walls?
    # Actually, E/W walls are S+2t.
    # So the physical outer width is S + 2t.
    
    tolerance = 1.0 # 1mm clearance
    lid_S = S + 2*t + tolerance
    lid_H = 25 # Fixed lid height
    
    lid_params = params.copy()
    lid_params['S'] = lid_S
    lid_params['H'] = lid_H
    # Disable features for lid
    lid_params['text_top'] = "LID"
    lid_params['text_front'] = ""
    lid_params['include_logo'] = False
    lid_params['include_fractal'] = False
    
    # Recursively use box generator logic for the lid
    return generate_box_svg(lid_params, filename)

def generate_divider_svg(params, filename="divider_square.svg"):
    """Generates a divider insert plate with a central slot."""
    t = params['t']
    S = params['S']
    H = params['H']
    
    # Divider dimensions fit inside S x S.
    clearance = 0.5
    div_w = S - clearance
    div_h = H - clearance
    
    slot_pos_ratio = 0.5 # Default to center
    
    stock_w = params['stock_w']
    stock_h = params['stock_h']
    
    svg = SVGGenerator(filename, stock_w, stock_h)
    
    cx, cy = stock_w/2, stock_h/2
    x = cx - div_w/2
    y = cy - div_h/2
    
    # Outline
    svg.add_rect(x, y, div_w, div_h, "CUT")
    
    # Slot (Rectangle cut)
    slot_x = x + div_w * slot_pos_ratio
    slot_w = t
    slot_h = div_h / 2
    
    svg.add_rect(slot_x - slot_w/2, y, slot_w, slot_h, "CUT")
    
    svg.save()
    return True


# --- CLI & Main ---

def get_float(prompt, default=None):
    try:
        val = input(f"{prompt} [{default}]: ")
        if not val and default is not None:
            return default
        return float(val)
    except ValueError:
        print("Invalid number.")
        return get_float(prompt, default)

def get_bool(prompt):
    val = input(f"{prompt} (y/n) [n]: ").lower()
    return val == 'y'

def validate_inputs(params):
    # thickness t must be less than min(S, H) / 2
    t = params['t']
    if t <= 0: return False, "Thickness must be positive."
    if params['S'] <= 0 or params['H'] <= 0: return False, "Dimensions must be positive."
    
    limit = min(params['S'], params['H']) / 2
    if t >= limit:
        return False, f"Thickness {t} is too large for dimensions (limit < {limit})."
    
    if params.get('tab_width', 0) > params['S']:
        return False, "Tab width larger than side S."
    if params.get('tab_depth', 0) > params['H']:
        return False, "Tab depth larger than height H."
        
    return True, ""

def run_examples():
    """Generates 4 example files for the report."""
    print("Generating examples...")
    
    # Example 1: Standard Box
    p1 = {
        'S': 100, 'H': 80, 't': 3,
        'stock_w': 600, 'stock_h': 400,
        'tab_width': 15, 'tab_depth': 10,
        'include_logo': True, 'include_fractal': False,
        'text_top': "Ex1 Top", 'text_front': "Front"
    }
    generate_box_svg(p1, "example1_basic.svg")
    
    # Example 2: Fractal Box
    p2 = p1.copy()
    p2['include_fractal'] = True
    p2['text_top'] = "Fractal Box"
    generate_box_svg(p2, "example2_fractal.svg")
    
    # Example 3: Lid and Divider
    p3 = p1.copy()
    p3['S'] = 80 # Smaller to fit
    p3['text_top'] = "Parts Box"
    generate_box_svg(p3, "example3_main.svg")
    generate_lid_svg(p3, "example3_lid.svg")
    generate_divider_svg(p3, "example3_div.svg")
    
    # Example 4: Rectangle
    generate_rectangle_svg(50, 25, "example4_rectangle.svg")
    
    print("Examples generated.")

def main():
    print("=== Parametric Box Generator ===")
    print("1. Generate Rectangle SVG (Task 1)")
    print("2. Generate Box SVGs")
    print("3. Run Example Suite")
    
    choice = input("Select > ")
    
    if choice == '1':
        w = get_float("Width (mm)", 40.0)
        h = get_float("Height (mm)", 30.0)
        generate_rectangle_svg(w, h)
        
    elif choice == '2':
        params = {}
        params['stock_w'] = get_float("Stock Width", 600)
        params['stock_h'] = get_float("Stock Height", 400)
        params['t'] = get_float("Thickness t", 3.0)
        params['S'] = get_float("Side S", 100.0)
        params['H'] = get_float("Height H", 50.0)
        params['tab_width'] = get_float("Tab Width", 20.0)
        params['tab_depth'] = get_float("Tab Depth", 15.0)
        
        params['text_top'] = input("Text on Top [enter for none]: ")
        params['text_front'] = input("Text on Front [enter for none]: ")
        
        params['include_logo'] = get_bool("Include Logo?")
        params['include_fractal'] = get_bool("Include Fractal?")
        
        do_lid = get_bool("Include Lid?")
        do_div = get_bool("Include Divider?")
        
        valid, msg = validate_inputs(params)
        if not valid:
            print("Error:", msg)
            return
            
        generate_box_svg(params, "box_square.svg")
        if do_lid:
            generate_lid_svg(params, "lid_square.svg")
        if do_div:
            generate_divider_svg(params, "divider_square.svg")
            
    elif choice == '3':
        run_examples()
    else:
        print("Invalid choice.")

if __name__ == "__main__":
    main()
