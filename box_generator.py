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

def draw_rollover_wall(svg, x_attach, y_top, y_bot, h, t, direction, slots_y):
    """Refactored helper to draw a rollover double wall with locking tabs."""
    # direction: -1 for Left, 1 for Right
    
    # Outer Wall
    x_fold = x_attach + direction * h
    svg.add_line(x_fold, y_top, x_fold, y_bot, "ENGRAVE") # Fold between Outer/Inner
    
    # Outer Wall Top/Bot Edges
    svg.add_line(x_attach, y_top, x_fold, y_top, "CUT")
    svg.add_line(x_attach, y_bot, x_fold, y_bot, "CUT")
    
    # Inner Wall
    inner_w = h - t # Slightly less
    x_end = x_fold + direction * inner_w
    
    # Inner Wall Top/Bot Edges
    svg.add_line(x_fold, y_top, x_end, y_top, "CUT")
    svg.add_line(x_fold, y_bot, x_end, y_bot, "CUT")
    
    # Edge with Locking Tabs aligned to slots_y
    # tab_h should match slot size (20 from mailer code)
    slot_h = 20
    tab_depth = 4 # sticking out
    x_tab = x_end + direction * tab_depth
    
    curr_y = y_top
    
    # Create locking tabs for each slot position
    # slots_y is a list of center Y coordinates for the tabs
    sorted_slots = sorted(slots_y)
    
    for slot_y_center in sorted_slots:
        t_y1 = slot_y_center - slot_h/2 + 1
        t_y2 = slot_y_center + slot_h/2 - 1
        
        # Line to start of tab
        svg.add_line(x_end, curr_y, x_end, t_y1, "CUT")
        # Tab shape
        svg.add_line(x_end, t_y1, x_tab, t_y1, "CUT")
        svg.add_line(x_tab, t_y1, x_tab, t_y2, "CUT")
        svg.add_line(x_tab, t_y2, x_end, t_y2, "CUT")
        
        curr_y = t_y2
        
    # Final segment to bottom
    svg.add_line(x_end, curr_y, x_end, y_bot, "CUT")



def generate_mailer_box_svg(params, filename="mailer_box.svg"):
    """Generates a detailed FEFCO 0427 mailer box with rollover side walls."""
    # Params
    S = params['S'] # Inner width/length (Square base)
    H = params['H'] # Inner height
    t = params['t'] # Material thickness
    stock_w = params['stock_w']
    stock_h = params['stock_h']
    
    # Topology:
    # Central Column: Lid Flap -> Lid -> Back Wall -> Base -> Front Wall (Double?)
    # Actually, 0427 standard:
    # - Base is S x S.
    # - Back Wall (H) -> Lid (S + ~t) -> Lid Flap (H inwards? or Tuck flap)
    # - Side Wings (Attached to Base): Outer Wall (H) -> Inner Wall (H-t) -> Locking Tabs
    # - Front Wall (Attached to Base): Double thickness usually? 
    #   Wait, looking at the blueprint:
    #   Bottom panel is "Front Wall". Included tabs on sides?
    #   Usually 0427 Front Wall is formed by the Side Wings folding over and locking?
    #   No, 0427 has a Front Wall panel attached to the base.
    #   The Side Wings have "ears" that tuck into the Front Wall?
    #   Image Analysis:
    #   - Bottom Panel (Front Wall): Single panel. Height ~H.
    #   - Side Wings have "ears" on the FRONT end. These ears tuck into the Front Wall (Double wall front).
    #   - BUT the blueprint shows Side Wings rolling over (Double Side Walls).
    #   
    #   Let's implement the robust "Rollover Side" 0427:
    #   1. Base (S x S)
    #   2. Back Wall (Attached to Base Top). Height H.
    #   3. Lid (Attached to Back Top). Depth S. Width S (plus clearance).
    #      - Includes Dust Flaps on sides (Tuck into side walls).
    #      - Includes Front Tuck Flap (Tucks into Front Wall).
    #   4. Front Wall (Attached to Base Bottom). Height H.
    #      - Has slots to receive Side Wing lugs?
    #      - Or is it double?
    #      Let's go with a simpler interpretation of the blueprint which seems to be:
    #      - Side Walls are double (Rollover). They have tabs that lock into the BASE slots.
    #      - Front Wall is attached to Base.
    #      - Lid tucks into Front Wall.
    
    # Dimensions
    # Base
    Base_W = S
    Base_D = S
    
    # Walls
    Back_H = H
    Front_H = H
    
    # Lid
    Lid_D = S # Covers the base
    # Lid Width: Needs to cover the Side Walls.
    # Side Walls are Double. Thickness = 2*t?
    # If Outer Wall is vertical, Inner Wall is inside.
    # Outer dimension is Base_W + 2*t.
    # So Lid Width = Base_W + 2*t.
    Lid_W = S + 2*t
    
    # Side Wings (Rollover)
    # Outer Side Height = H
    # Inner Side Height = H - t (to sit on base)
    # Locking Lugs on Inner Side bottom edge.
    
    # Layout Calculation
    # Width = LeftSide (H + H) + Base (S) + RightSide (H + H) = S + 4H
    # Height = LidFlap + Lid + Back + Base + Front
    
    # Check Stock
    total_w = S + 4 * H + 20 # margins
    total_h = H + Lid_D + Back_H + Base_D + Front_H + 20
    
    if total_w > stock_w or total_h > stock_h:
        print(f"Error: Layout size ({total_w:.1f} x {total_h:.1f} mm) larger than stock ({stock_w} x {stock_h} mm).")
        return False
        
    svg = SVGGenerator(filename, stock_w, stock_h)
    
    cx = stock_w / 2
    cy = stock_h / 2
    
    # Y-Flow (Top to Bottom)
    # We center the Base
    y_base_top = cy - S/2
    y_base_bot = cy + S/2
    
    bx = cx - S/2
    by = y_base_top
    
    # --- 1. Base ---
    # Folds
    svg.add_line(bx, by, bx + S, by, "ENGRAVE") # Back
    svg.add_line(bx, by + S, bx + S, by + S, "ENGRAVE") # Front
    svg.add_line(bx, by, bx, by + S, "ENGRAVE") # Left
    svg.add_line(bx + S, by, bx + S, by + S, "ENGRAVE") # Right
    
    # Slots for Rollover Tabs
    # Located in Base, near Left/Right edges.
    # Tab Position: Centered or specific? 
    # Usually 2 slots per side or 1 long one.
    # Blueprint has rectangles. Let's put 2 slots per side.
    slot_w = 4 # Width of slot (thickness accommodation)
    slot_h = 20 # Length of slot
    slot_inset = 2 # From edge
    
    # Left Slots
    ls_x = bx + slot_inset
    svg.add_rect(ls_x, by + S/4 - slot_h/2, slot_w, slot_h, "CUT")
    svg.add_rect(ls_x, by + 3*S/4 - slot_h/2, slot_w, slot_h, "CUT")
    
    # Right Slots
    rs_x = bx + S - slot_inset - slot_w
    svg.add_rect(rs_x, by + S/4 - slot_h/2, slot_w, slot_h, "CUT")
    svg.add_rect(rs_x, by + 3*S/4 - slot_h/2, slot_w, slot_h, "CUT")
    
    # --- 2. Back Wall ---
    # [bx, bx+S] x [by-H, by]
    y_back_top = by - H
    svg.add_line(bx, y_back_top, bx + S, y_back_top, "ENGRAVE") # Fold to Lid
    svg.add_line(bx, y_back_top, bx, by, "CUT") # Left Edge
    svg.add_line(bx + S, y_back_top, bx + S, by, "CUT") # Right Edge
    
    # --- 3. Lid ---
    # [lid_x_left, lid_x_right] x [y_lid_top, y_back_top]
    # Lid Width = S + 2t. Centered on Back.
    lid_x_left = cx - Lid_W/2
    lid_x_right = cx + Lid_W/2
    y_lid_top = y_back_top - Lid_D
    
    # Connection Back->Lid
    svg.add_line(bx, y_back_top, lid_x_left, y_back_top, "CUT")
    svg.add_line(bx + S, y_back_top, lid_x_right, y_back_top, "CUT")
    
    # Lid Side Folds (for Dust Flaps)
    svg.add_line(lid_x_left, y_lid_top, lid_x_left, y_back_top, "ENGRAVE")
    svg.add_line(lid_x_right, y_lid_top, lid_x_right, y_back_top, "ENGRAVE")
    
    # --- 4. Lid Dust Flaps ---
    # Attached to Lid Sides.
    # Width H - 2 (clearance). Depth Lid_D (tapered).
    df_w = H - 2
    # Left
    svg.add_line(lid_x_left - df_w, y_lid_top + 5, lid_x_left - df_w, y_back_top - 5, "CUT") # Outer
    svg.add_line(lid_x_left - df_w, y_lid_top + 5, lid_x_left, y_lid_top, "CUT") # Top Diag
    svg.add_line(lid_x_left - df_w, y_back_top - 5, lid_x_left, y_back_top, "CUT") # Bot Diag
    
    # Right
    svg.add_line(lid_x_right + df_w, y_lid_top + 5, lid_x_right + df_w, y_back_top - 5, "CUT")
    svg.add_line(lid_x_right + df_w, y_lid_top + 5, lid_x_right, y_lid_top, "CUT")
    svg.add_line(lid_x_right + df_w, y_back_top - 5, lid_x_right, y_back_top, "CUT")
    
    # --- 5. Lid Tuck Flap ---
    # Attached to Lid Top (y_lid_top).
    # Height ~ H.
    # Tucks into Front Wall.
    tf_h = H - 5
    y_flap_end = y_lid_top - tf_h
    
    # Fold
    svg.add_line(lid_x_left, y_lid_top, lid_x_right, y_lid_top, "ENGRAVE")
    
    # Shape: Locking Ears ("Cherry Locks")
    # Tapers in, then bumps out, then round.
    # Simplified: Tapered rect.
    flap_w_start = Lid_W
    flap_w_end = Lid_W - 10
    fx_left = cx - flap_w_end/2
    fx_right = cx + flap_w_end/2
    
    svg.add_line(lid_x_left, y_lid_top, fx_left, y_flap_end, "CUT")
    svg.add_line(lid_x_right, y_lid_top, fx_right, y_flap_end, "CUT")
    svg.add_line(fx_left, y_flap_end, fx_right, y_flap_end, "CUT")
    
    # --- 6. Front Wall ---
    # Attached to Base Bottom (y_base_bot).
    # Height H. Width S.
    y_front_bot = y_base_bot + H
    
    # Sides (Locking slots for Side Wings?)
    # Usually standard 0427 Front Wall has folded ends that lock with side wings.
    # Simple version: Rectangular panel.
    svg.add_line(bx, y_base_bot, bx, y_front_bot, "CUT")
    svg.add_line(bx + S, y_base_bot, bx + S, y_front_bot, "CUT")
    svg.add_line(bx, y_front_bot, bx + S, y_front_bot, "CUT")
    
    # --- 7. Rollover Side Wings ---
    # Attached to Base Left/Right.
    # Structure: Outer Wall (width H) -> Fold -> Inner Wall (width H-t) -> Tabs.
    
    # --- 7. Rollover Side Wings (using helper) ---
    # Attached to Base Left/Right.
    
    # Calculate slot positions (same as derived in Base section)
    # s1_y = by + S/4
    # s2_y = by + 3*S/4
    s1_y = by + S/4
    s2_y = by + 3*S/4
    
    draw_rollover_wall(svg, bx, by, by + S, H, t, -1, [s1_y, s2_y]) # Left
    draw_rollover_wall(svg, bx + S, by, by + S, H, t, 1, [s1_y, s2_y]) # Right

    # --- Text ---
    if params.get('text_top'):
        svg.add_text(cx, (y_lid_top + y_back_top)/2, params['text_top'], font_size=10, mode="ENGRAVE")

    svg.save()
    return True



def generate_tray_svg(params, filename="tray_rollover.svg"):
    """Generates a lidless tray with double rollover side walls (FEFCO 0422 style)."""
    S = params['S']
    H = params['H']
    t = params['t']
    stock_w = params['stock_w']
    stock_h = params['stock_h']
    
    # Layout Check
    total_w = S + 4*H + 20
    total_h = S + 2*(H + H) + 20 # S + 4H roughly
    
    if total_w > stock_w or total_h > stock_h:
        print(f"Error: Layout size ({total_w:.1f} x {total_h:.1f} mm) larger than stock ({stock_w} x {stock_h} mm).")
        return False
    
    svg = SVGGenerator(filename, stock_w, stock_h)
    
    cx = stock_w / 2
    cy = stock_h / 2
    
    # Coordinates
    bx = cx - S/2
    by = cy - S/2
    
    # --- Base ---
    svg.add_line(bx, by, bx + S, by, "ENGRAVE") 
    svg.add_line(bx, by + S, bx + S, by + S, "ENGRAVE")
    svg.add_line(bx, by, bx, by + S, "ENGRAVE")
    svg.add_line(bx + S, by, bx + S, by + S, "ENGRAVE")
    
    # Slots (Same as mailer)
    slot_w = 4; slot_h = 20; slot_inset = 2
    s1_y = by + S/4
    s2_y = by + 3*S/4
    
    ls_x = bx + slot_inset
    svg.add_rect(ls_x, s1_y - slot_h/2, slot_w, slot_h, "CUT")
    svg.add_rect(ls_x, s2_y - slot_h/2, slot_w, slot_h, "CUT")
    
    rs_x = bx + S - slot_inset - slot_w
    svg.add_rect(rs_x, s1_y - slot_h/2, slot_w, slot_h, "CUT")
    svg.add_rect(rs_x, s2_y - slot_h/2, slot_w, slot_h, "CUT")
    
    # --- Front/Back Walls with Corner Flaps ---
    
    def draw_end_wall_assembly(y_hinge, is_top):
        # direction_y: -1 for Top (up), 1 for Bot (down)
        dy = -1 if is_top else 1
        y_edge = y_hinge + dy * H
        
        # Left Corner Flap
        flap_w = H - 2*t # clearance
        flap_x_outer = bx - flap_w
        
        # Fold Main Wall - Flap (Actually fold base-wall is done)
        # We need folds between Wall and Flaps?
        # Yes, corner flaps are attached to Left/Right of the Wall panel.
        
        # Wall Panel Vertical Edges (Folds to flaps)
        svg.add_line(bx, y_hinge, bx, y_edge, "ENGRAVE")
        svg.add_line(bx + S, y_hinge, bx + S, y_edge, "ENGRAVE")
        
        # Draw Left Flap
        # y positions are same as wall: y_hinge to y_edge
        svg.add_line(flap_x_outer, y_hinge, bx, y_hinge, "CUT") # Side (near base)
        svg.add_line(flap_x_outer, y_hinge, flap_x_outer, y_edge, "CUT") # Outer
        svg.add_line(flap_x_outer, y_edge, bx, y_edge, "CUT") # Top/Bot edge
        
        # Draw Right Flap
        flap_x_outer_r = bx + S + flap_w
        svg.add_line(bx + S, y_hinge, flap_x_outer_r, y_hinge, "CUT")
        svg.add_line(flap_x_outer_r, y_hinge, flap_x_outer_r, y_edge, "CUT")
        svg.add_line(flap_x_outer_r, y_edge, bx + S, y_edge, "CUT")
        
        # Main Wall Top/Bot Edge (Outer edge)
        svg.add_line(bx, y_edge, bx + S, y_edge, "CUT")

    # Draw Back Assembly (Top)
    draw_end_wall_assembly(by, True)
    
    # Draw Front Assembly (Bottom)
    draw_end_wall_assembly(by + S, False)
    
    # --- Rollover Wings ---
    draw_rollover_wall(svg, bx, by, by + S, H, t, -1, [s1_y, s2_y]) # Left
    draw_rollover_wall(svg, bx + S, by, by + S, H, t, 1, [s1_y, s2_y]) # Right
    
    svg.save()
    return True

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
        'stock_w': 600, 'stock_h': 600,
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
    
    # Mailer Example
    p_mail = p1.copy()
    p_mail['S'] = 150
    p_mail['H'] = 50
    p_mail['text_top'] = "Mailer Box"
    generate_mailer_box_svg(p_mail, "example_mailer.svg")
    
    # Tray Example
    generate_tray_svg(p_mail, "example_tray.svg")
    
    print("Examples generated.")

def main():
    print("=== Parametric Box Generator ===")
    print("1. Generate Rectangle SVG (Task 1)")
    print("2. Generate Box SVGs (Standard Open Bin)")
    print("3. Generate Mailer Box SVG (FEFCO 0427)")
    print("4. Generate Tray SVG (Lidless 0427)")
    print("5. Run Example Suite")
    
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
        params = {}
        params['stock_w'] = get_float("Stock Width", 600)
        params['stock_h'] = get_float("Stock Height", 600) # Need more height for mailer
        params['t'] = get_float("Thickness t", 3.0)
        params['S'] = get_float("Side S", 150.0)
        params['H'] = get_float("Height H", 50.0)
        params['text_top'] = input("Text on Lid [enter for none]: ")
        
        generate_mailer_box_svg(params, "mailer_box.svg")
        
    elif choice == '4':
        params = {}
        params['stock_w'] = get_float("Stock Width", 600)
        params['stock_h'] = get_float("Stock Height", 600) 
        params['t'] = get_float("Thickness t", 3.0)
        params['S'] = get_float("Side S", 150.0)
        params['H'] = get_float("Height H", 50.0)
        generate_tray_svg(params, "tray_rollover.svg")
        
    elif choice == '5':
        run_examples()
    else:
        print("Invalid choice.")


if __name__ == "__main__":
    main()
