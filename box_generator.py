#!/usr/bin/env python3
"""
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
    """Generates the acrylic finger-joint box layout."""
    # Unpack params
    S = params['S']
    H = params['H']
    t = params['t']
    stock_w = params['stock_w']
    stock_h = params['stock_h']

    # --- New Finger Joint Logic ---

    # We need to draw 5 parts: Base, Front, Back, Left, Right.
    # Layout Strategy:
    # Row 1: Front, Back
    # Row 2: Left, Right, Base
    
    # Margin between parts
    margin = 5
    
    # Dimensions:
    # Base: S x S (Inner)
    # Front/Back: (S + 2t) x H (Covers corners)
    # Left/Right: S x H (Fits between Front/Back)
    
    # Parts List with Sizes (W, H)
    parts = {}
    
    # 1. Base (S x S)
    parts['Base'] = {'w': S, 'h': S}
    
    # 2. Front (S+2t x H)
    parts['Front'] = {'w': S + 2*t, 'h': H}
    parts['Back']  = {'w': S + 2*t, 'h': H}
    
    # 3. Side (S x H)
    parts['Left']  = {'w': S, 'h': H}
    parts['Right'] = {'w': S, 'h': H}

    # Layout Calculation
    max_row_h = H
    total_w_row1 = parts['Front']['w'] + parts['Back']['w'] + parts['Left']['w'] + parts['Right']['w'] + 4*margin
    # Row 1: Front, Back
    # Row 2: Left, Right, Base, Divider?
    
    x_cursor = margin
    y_cursor = margin
    
    # Front
    parts['Front']['x'] = x_cursor
    parts['Front']['y'] = y_cursor
    x_cursor += parts['Front']['w'] + margin
    
    # Back
    parts['Back']['x'] = x_cursor
    parts['Back']['y'] = y_cursor
    
    # Next Row
    y_cursor += H + margin
    x_cursor = margin
    
    # Left
    parts['Left']['x'] = x_cursor
    parts['Left']['y'] = y_cursor
    x_cursor += parts['Left']['w'] + margin
    
    # Right
    parts['Right']['x'] = x_cursor
    parts['Right']['y'] = y_cursor
    x_cursor += parts['Right']['w'] + margin
    
    # Base
    parts['Base']['x'] = x_cursor
    parts['Base']['y'] = y_cursor
    x_cursor += parts['Base']['w'] + margin

    # Divider
    if 'Divider' in parts:
        parts['Divider']['x'] = x_cursor
        parts['Divider']['y'] = y_cursor
        # Check stock width
        if x_cursor + parts['Divider']['w'] > stock_w:
            # Move to Row 3 if needed
            x_cursor = margin
            y_cursor += max(parts['Base']['h'], parts['Left']['h']) + margin # This should be max height of previous row
            parts['Divider']['x'] = x_cursor
            parts['Divider']['y'] = y_cursor
            
    # Check Layout Bounds
    max_x = 0
    max_y = 0
    for p in parts.values():
        max_x = max(max_x, p['x'] + p['w'])
        max_y = max(max_y, p['y'] + p['h'])
    
    if max_x + margin > stock_w or max_y + margin > stock_h:
        print(f"Error: Layout size ({max_x+margin:.1f} x {max_y+margin:.1f} mm) larger than stock ({stock_w} x {stock_h} mm).")
        return False
        
    svg = SVGGenerator(filename, stock_w, stock_h)
    
    # Finger parameters
    tab_size = 10 # nominal 10mm tabs
    
    def draw_finger_edge(x1, y1, x2, y2, thickness, parity, mode="CUT"):
        """
        Draws a finger-jointed line from (x1, y1) to (x2, y2).
        parity: 1 = start with TAB (stick out), -1 = start with SLOT (dent in).
        Direction: currently horizontal or vertical only.
        thickness: t (depth of tab/slot).
        """
        dx = x2 - x1
        dy = y2 - y1
        length = math.sqrt(dx*dx + dy*dy)
        if length == 0: return

        # Calculate tabs (ensure odd number for symmetry)
        n_tabs = max(3, int(length // tab_size))
        if n_tabs % 2 == 0: n_tabs += 1
        actual_tab_w = length / n_tabs
        
        # Geometry vectors
        ux = dx / length
        uy = dy / length
        
        # Initialize traversal
        curr_x, curr_y = x1, y1
        side = parity # 1 = Out, -1 = In
        points = [(x1, y1)]

        # Finger Joint Logic:
        # Parity 1 (Protruding Tab): [Tab (+t), Gap (0), Tab (+t)...]
        # Parity -1 (Recessed Slot): [Gap (0), Tab (+t), Gap (0)...]
        # Note: "Left" of vector (x1->x2) defines the protrusion direction.
        
        for i in range(n_tabs):
            # Calculate segment endpoint
            p_end_center = (x1 + (i+1)*actual_tab_w*ux, y1 + (i+1)*actual_tab_w*uy)
            
            current_offset_mag = thickness if side == 1 else 0
            
            # Calc start/end of this segment at the offset
            seg_start_x = curr_x + (-uy * current_offset_mag) # Move Left
            seg_start_y = curr_y + ( ux * current_offset_mag)
            
            seg_end_x = p_end_center[0] + (-uy * current_offset_mag)
            seg_end_y = p_end_center[1] + ( ux * current_offset_mag)
            
            # If this is not the first point, adding a perpendicular connector from previous
            if i > 0:
                svg.add_line(points[-1][0], points[-1][1], seg_start_x, seg_start_y, mode)
                
            # If first point and offset is non-zero, we surely need to connect from x1,y1?
            # Actually we usually start FROM the corner.
            # If Parity 1 (Start Tab +t): Do we start at corner or corner+t?
            # Corner is part dimension.
            # Base S x S.
            # Corner is at 0. Tab starts at corner?
            # Yes, corner is inside the tab? Or tab starts after corner?
            # Standard: Corner IS part of a tab (stronger).
            # So we start at (x1,y1) + offset.
            # But the "Line" starts at x1,y1.
            # So we need to draw (x1,y1) -> Start.
            if i == 0:
                 svg.add_line(x1, y1, seg_start_x, seg_start_y, mode)
            
            # Add segment
            svg.add_line(seg_start_x, seg_start_y, seg_end_x, seg_end_y, mode)
            points.append((seg_end_x, seg_end_y))
            
            # Flip side for next segment
            side = -side
            curr_x, curr_y = p_end_center[0], p_end_center[1]
            
        # Connect to final x2,y2
        svg.add_line(points[-1][0], points[-1][1], x2, y2, mode)


    # --- Draw Parts ---

    # 1. Base (S x S)
    # Parity: 1 (Tabs on all sides to stick out)
    bx, by = parts['Base']['x'], parts['Base']['y']
    bw, bh = parts['Base']['w'], parts['Base']['h']
    
    # Top (L->R): Parity 1
    draw_finger_edge(bx, by, bx+bw, by, t, 1, "CUT")
    # Right (T->B): Parity 1
    draw_finger_edge(bx+bw, by, bx+bw, by+bh, t, 1, "CUT")
    # Bottom (R->L): Parity 1
    draw_finger_edge(bx+bw, by+bh, bx, by+bh, t, 1, "CUT")
    # Left (B->T): Parity 1
    draw_finger_edge(bx, by+bh, bx, by, t, 1, "CUT")
    
    # Engraving on Base
    if params.get('text_top'):
        svg.add_text(bx + bw/2, by + bh/2, params['text_top'], font_size=8, mode="ENGRAVE")

    # Edge Logic for Front/Back (Bottom):
    # Base fits INSIDE the walls.
    # Base has Tabs (Parity 1) on all sides.
    # Front/Back Bottom Edge must accept Base Tabs with compatible Slot pattern.
    # The Mating Zone corresponds to the center 'S' width. 
    # Corners (size t) are solid to cover the Left/Right wall edges.
    
    def draw_wall_bottom_with_base_slots(px, py, pw, ph):
        # Path: Start (Bottom Left) -> t (Corner) -> S (Mating Slots) -> t (Corner) -> End
        p_y_bot = py + ph
        # 1. Left Corner
        svg.add_line(px, p_y_bot, px+t, p_y_bot, "CUT")
        # 2. Mating Zone (Parity -1 to Rx Base Tabs)
        draw_finger_edge(px+t, p_y_bot, px+t+S, p_y_bot, t, -1, "CUT")
        # 3. Right Corner
        svg.add_line(px+t+S, p_y_bot, px+t+S+t, p_y_bot, "CUT")


    # 2. Front (S+2t x H)
    # Location: Below Base? No, we set layout earlier.
    fx, fy = parts['Front']['x'], parts['Front']['y']
    fw, fh = parts['Front']['w'], parts['Front']['h']
    
    # Front Wall Draw
    # Top: Flat
    svg.add_line(fx, fy, fx+fw, fy, "CUT")
    # Sides: Mating with Left/Right Walls.
    draw_finger_edge(fx+fw, fy, fx+fw, fy+fh, t, 1, "CUT") # Right Side (Down)
    draw_wall_bottom_with_base_slots(fx, fy, fw, fh)       # Bottom
    draw_finger_edge(fx, fy+fh, fx, fy, t, 1, "CUT")       # Left Side (Up)
    
    # Divider Slot (Vertical) on Front Wall
    # Position is ratio along the Width (S).
    if params.get('divider_pos') is not None:
        # Front Wall Width is S+2t. The Inner S part is from fx+t to fx+t+S.
        pos_ratio = params['divider_pos']
        # Slot center x
        slot_cx = fx + t + S * pos_ratio
        
        # Slot dimensions (H/2 height, t width) centered vertically relative to H
        slot_h_cut = H / 2
        slot_y_start = fy + (H - slot_h_cut)/2
        
        svg.add_rect(slot_cx - t/2, slot_y_start, t, slot_h_cut, "CUT")

    # Logo on Front
    if params.get('include_logo'):
        cx_f = fx + fw/2
        cy_f = fy + fh/2
        svg.add_text(cx_f, cy_f - 4, "Columbia", font_size=6, mode="ENGRAVE")
        svg.add_text(cx_f, cy_f + 4, "Digital Manufacturing", font_size=4, mode="ENGRAVE")
    if params.get('text_front'):
         svg.add_text(fx + fw/2, fy + fh/2 - 15, params['text_front'], font_size=6, mode="ENGRAVE")


    # 3. Back Wall (Same as Front)
    bax, bay = parts['Back']['x'], parts['Back']['y']
    baw, bah = parts['Back']['w'], parts['Back']['h']
    
    svg.add_line(bax, bay, bax+baw, bay, "CUT") # Top
    draw_finger_edge(bax+baw, bay, bax+baw, bay+bah, t, 1, "CUT") # Right
    draw_wall_bottom_with_base_slots(bax, bay, baw, bah)          # Bottom
    draw_finger_edge(bax, bay+bah, bax, bay, t, 1, "CUT")         # Left
    
    # Divider Slot on Back Wall
    if params.get('divider_pos') is not None:
        pos_ratio = params['divider_pos']
        # Back is mirrored? Logic is same, X from Left is X from Left.
        # If Divider is at 0.3 from Left (West).
        # Front Wall: 0.3 from Left.
        # Back Wall: 0.3 from Left (if viewed from outside).
        # Inner S starts at t.
        slot_cx = bax + t + S * pos_ratio
        slot_h_cut = H / 2
        slot_y_start = bay + (H - slot_h_cut)/2
        svg.add_rect(slot_cx - t/2, slot_y_start, t, slot_h_cut, "CUT")


    # 4. Left Side Wall (S x H)
    lx, ly = parts['Left']['x'], parts['Left']['y']
    lw, lh = parts['Left']['w'], parts['Left']['h']
    
    # Mates base at bottom. Mates Front/Back at sides.
    # Bottom: Mates Base (Parity -1). Direct S length.
    # Sides: Mates Front/Back (Parity 1). So L/R needs Parity -1.
    
    svg.add_line(lx, ly, lx+lw, ly, "CUT") # Top Flat
    draw_finger_edge(lx+lw, ly, lx+lw, ly+lh, t, -1, "CUT") # Right Side (Matches Front Left 1)
    draw_finger_edge(lx+lw, ly+lh, lx, ly+lh, t, -1, "CUT") # Bot Side (Matches Base 1)
    draw_finger_edge(lx, ly+lh, lx, ly, t, -1, "CUT")       # Left Side (Matches Back Right 1)


    # 5. Right Side Wall
    rx, ry = parts['Right']['x'], parts['Right']['y']
    rw, rh = parts['Right']['w'], parts['Right']['h']
    
    svg.add_line(rx, ry, rx+rw, ry, "CUT")
    draw_finger_edge(rx+rw, ry, rx+rw, ry+rh, t, -1, "CUT")
    draw_finger_edge(rx+rw, ry+rh, rx, ry+rh, t, -1, "CUT")
    draw_finger_edge(rx, ry+rh, rx, ry, t, -1, "CUT")
    
    # 6. Divider (if enabled)
    if 'Divider' in parts:
        dx, dy = parts['Divider']['x'], parts['Divider']['y']
        dw, dh = parts['Divider']['w'], parts['Divider']['h'] # Total W, H
        
        # Logic from generate_divider_svg (Tabbed)
        clearance = 0.5
        body_w = S - clearance
        body_h = H - t - clearance
        
        tab_h = (H / 2) - clearance
        tab_w = t
        
        #
        
        # Vertical Alignment: 
        # Wall Slot Start (from Top 0): H/4.
        # Divider Top (flush with Wall Top).
        # So Divider Tab Start (from Top 0): H/4 + clearance/2.
        
        slot_y_start_rel = (H - (H/2)) / 2 # H/4
        tab_y_start_rel = slot_y_start_rel + (clearance / 2)
        tab_y_end_rel = tab_y_start_rel + tab_h
        
        # Origin (Top-Left of Body, not total). 
        # Let's map points relative to dx, dy (Top Left of Total Bounding Box).
        # Total Box: x=dx, y=dy.
        # Body X: dx + tab_w.
        # Body Y: dy (Flush Top? No, Divider height is H - t - clearance). 
        # Wall Slot is centered on H. Divider should sit on Base (at H-t).
        # But we want divider flush with TOP.
        # So Divider Top = Wall Top.
        # So Divider Body Y = dy.
        
        ox = dx + tab_w
        oy = dy
        
        # Points Logic (Clockwise from Top-Left of Body)
        # 1. Top Edge
        svg.add_line(ox, oy, ox + body_w, oy, "CUT")
        
        # 2. Right Side Upper
        svg.add_line(ox + body_w, oy, ox + body_w, oy + tab_y_start_rel, "CUT")
        
        # 3. Right Tab
        svg.add_line(ox + body_w, oy + tab_y_start_rel, ox + body_w + tab_w, oy + tab_y_start_rel, "CUT") # Out
        svg.add_line(ox + body_w + tab_w, oy + tab_y_start_rel, ox + body_w + tab_w, oy + tab_y_end_rel, "CUT") # Down
        svg.add_line(ox + body_w + tab_w, oy + tab_y_end_rel, ox + body_w, oy + tab_y_end_rel, "CUT") # Back
        
        # 4. Right Side Lower
        svg.add_line(ox + body_w, oy + tab_y_end_rel, ox + body_w, oy + body_h, "CUT")
        
        # 5. Bottom Edge
        svg.add_line(ox + body_w, oy + body_h, ox, oy + body_h, "CUT")
        
        # 6. Left Side Lower
        svg.add_line(ox, oy + body_h, ox, oy + tab_y_end_rel, "CUT")
        
        # 7. Left Tab
        svg.add_line(ox, oy + tab_y_end_rel, ox - tab_w, oy + tab_y_end_rel, "CUT") # Out (Left)
        svg.add_line(ox - tab_w, oy + tab_y_end_rel, ox - tab_w, oy + tab_y_start_rel, "CUT") # Up
        svg.add_line(ox - tab_w, oy + tab_y_start_rel, ox, oy + tab_y_start_rel, "CUT") # Back
        
        # 8. Left Side Upper
        svg.add_line(ox, oy + tab_y_start_rel, ox, oy, "CUT")

    svg.save()
    return True


# --- Wrapper & CLI ---

def print_software_description(params):
    """Prints the description of the software as required by rubric."""
    print("\n--- Software Description (Acrylic) ---")
    print("1. Calculation Steps:")
    print("   - Inputs: Box Side (S), Height (H), Material Thickness (t).")
    print("   - Design: 5-Part Finger Joint Box (Base, Front, Back, Left, Right).")
    print("   - Fingers: Calculated dynamically based on edge length (aiming for 10mm tabs).")
    print("   - Base Size: S x S (Tabs Out).")
    print("   - Front/Back Size: (S + 2t) x H. Covers corners.")
    print("2. Conditions:")
    print(f"   - Thickness t ({params.get('t')}mm) used for finger depth.")
    print("   - Layout checked against stock.")
    print("3. Formulas:")
    print("   - Tab Depth = t")
    print("   - Tab Parity ensures A-fits-B (1 vs -1).")
    print("----------------------------\n")

def get_float(prompt, default=None):
    try:
        prompt_str = f"{prompt} [{default}]: " if default is not None else f"{prompt}: "
        val = input(prompt_str)
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
    t = params['t']
    if t <= 0: return False, "Thickness must be positive."
    if params['S'] <= 0 or params['H'] <= 0: return False, "Dimensions must be positive."
    limit = min(params['S'], params['H']) / 2
    if t >= limit:
        return False, f"Thickness {t} is too large."
    return True, ""

def main():
    print("=== Parametric Acrylic Box Generator ===")
    print("1. Task 1: Generate Small Rectangle Check (<= 2x2 inches)")
    print("2. Task 2: Generate Acrylic Box (Finger Joints)")
    
    choice = input("\nSelect > ")
    
    if choice == '1':
        print("\n-- Generating Small Rectangle --")
        w = get_float("Width (mm) [Max 50.8]", 50.0)
        h = get_float("Height (mm) [Max 50.8]", 25.0)
        generate_rectangle_svg(w, h, "task1_rectangle.svg")
        
    elif choice == '2':
        print("\n-- Acrylic Box Configuration --")
        params = {}
        params['stock_w'] = get_float("Stock Width (mm)", 600)
        params['stock_h'] = get_float("Stock Height (mm)", 600)
        params['t'] = get_float("Material Thickness t (mm)", 3.0)
        
        print("\n-- Dimensions --")
        params['S'] = get_float("Box Side Length S (mm) [Inner]", 100.0)
        params['H'] = get_float("Box Height H (mm)", params['S'])
        
        print("\n-- Features --")
        do_div = input("Include Divider Insert? (y/n) [y]: ").lower() != 'n'
        if do_div:
            params['divider_pos'] = get_float("Divider Position Ratio (0.0-1.0, 0.5=Center)", 0.5)
        
        print("\n-- Engraving --")
        # Default Base Text
        default_base_text = "Digital Manufacturing"
        t_top = input(f"Text on Base [{default_base_text}]: ")
        params['text_top'] = t_top if t_top else default_base_text
        
        # Disable Front Text prompt (keep clean)
        params['text_front'] = ""
        
        params['include_logo'] = get_bool("Engrave Columbia Logo (Front Wall)?")
        # params['include_fractal'] removed per user request
        
        # Validation
        valid, msg = validate_inputs(params)
        if not valid:
            print(f"\nError: {msg}")
            return
            
        # Execution
        print("\nGenerating files base on Acrylic Finger Joints...")
        
        # Reuse 'box_main.svg' filename or specific
        success = generate_box_svg(params, "box_acrylic_parts.svg")
        
        if success:
            print_software_description(params)
            
    else:
        print("Invalid selection.")

if __name__ == "__main__":
    main()
