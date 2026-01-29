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
    # Row 1: Left, Front, Right, Back (Walls laid out in a strip? No, optimization isn't critical but clarity is.)
    # Let's do:
    # Base in center? Or simple grid.
    # Base (S x S)
    # Front (S x H)
    # Back (S x H)
    # Left (S x H)
    # Right (S x H)
    
    # Margin between parts
    margin = 5
    
    # Calculate Part Origins (Top-Left of each part's bounding box)
    # Base:
    # Side length S.
    # Walls have thickness t. 
    # To interlock:
    # Base is S x S inner? 
    # A standard finger joint box usually defines OUTSIDE dimensions or INSIDE.
    # If S is INSIDE:
    # Base Width = S + (2*t if joints overlap).
    # Actually, simpler: 
    # Base = S x S.
    # Front/Back sit ON SIDE of Base? Or Base sits inside?
    # Common sturdy design: Base is "caught" by walls. 
    # Sides overlap Base -> Base size is S x S. 
    # Front/Back Width = S + 2*t. Height = H.
    # Left/Right Width = S. Height = H.
    # Wait, if Front/Back overlap Left/Right corners?
    # Let's assume standard "pinwheel" or "symmetric" corners for simplicity?
    # Or "Front/Back overlap Left/Right".
    # Case: Front/Back Width = S + 2*t. Left/Right Width = S.
    # Base Width along Front/Back = S. (Tabs insert into F/B faces).
    # Base Depth along Left/Right = S. (Tabs insert into L/R faces).
    
    # Dimensions:
    # Base: S x S (with tabs sticking OUT).
    # Front/Back: (S + 2t) x H. Slots for Base. Slots for Side Walls.
    # Left/Right: S x H. Slots for Base. Tabs for Front/Back.
    
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
    
    # If row1 is too wide, split.
    # Let's do 2 Rows:
    # Row 1: Front, Back
    # Row 2: Left, Right, Base
    
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
    
    # Check Layout
    total_w_used = max(parts['Back']['x'] + parts['Back']['w'], parts['Base']['x'] + parts['Base']['w']) + margin
    total_h_used = y_cursor + S + margin
     
    if total_w_used > stock_w or total_h_used > stock_h:
        print(f"Error: Layout size ({total_w_used:.1f} x {total_h_used:.1f} mm) larger than stock ({stock_w} x {stock_h} mm).")
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

        # Determine number of tabs
        # ideally odd number to be symmetric? 
        n_tabs = max(3, int(length // tab_size))
        if n_tabs % 2 == 0: n_tabs += 1 # Ensure odd for symmetry
        
        actual_tab_w = length / n_tabs
        
        # Unit vector
        ux = dx / length
        uy = dy / length
        # Normal vector (pointing 'out' relative to standard CCW traversal, but we just need perpendicular)
        # For a horizontal line x1->x2 (ux=1, uy=0), 'out' is variable?
        # Let's explicitly define 'out' as +90 deg or consistent based on parity?
        # Simpler: 'parity=1' means the first segment sticks OUT to the 'Right' of the path vector.
        
        # Normal vector n = (-uy, ux) (Left turn)
        # If parity=1 (Tab), we go OUT (Left of path?).
        # Wait, usually for a box part boundary traversing CW:
        # Tab sticks OUT (into neighbor's space).
        # Slot cuts IN (into own space).
        
        nx = -uy
        ny = ux
        
        # Points list
        curr_x, curr_y = x1, y1
        side = parity # 1 = Out, -1 = In
        
        # If we start with parity=1 (Tab), the first "edge" is extended OUT.
        # Structure:
        # [Segment along edge] -> [90deg Out] -> [Top of Tab] -> [90deg In] -> ...
        # Actually standard finger joint:
        # 0 to w: Edge
        # w to 2w: Tab Top (at distance t)
        # ...
        # If parity=1 (Start with TAB):
        # We start "UP" at t? Or start at 0 and go UP?
        # "Start with Tab" usually means the first interval is a Tab.
        # So we travel along the protruded line.
        
        points = [] # list of (x,y)
        points.append((x1, y1))
        
        for i in range(n_tabs):
            # Each segment is actual_tab_w long.
            # State of this segment determined by side.
            # side 1: Tab (Protruding).
            # side -1: Slot (Recessed).
            
            # If we are side 1: We should be at distance 'thickness' from centerline?
            # Or is the centerline the mating line? Yes.
            # So Tab = (+nx*t, +ny*t). Slot = (-nx*t, -ny*t)? 
            # OR Slot = (0,0) (on line) and Tab = (Out)?
            # Usually: Line is the boundary.
            # Tab sticks OUT -> Parity 1 -> Go Out, Move, Go In.
            # Slot stays IN -> Parity -1 -> Stay on line (or cut IN?).
            # "Finger Joint": Teeth cross the line. 
            # Convention:
            # 1 (Tab): The material EXISTS outside the line. Cut path goes OUT.
            # -1 (Slot): The material is REMOVED inside the line. Cut path goes IN?
            # Actually, if we just want "interlock", one part has tabs OUT, other has tabs OUT (offset).
            # Let's define Parity 1 = Start with PROTRUSION (Tab).
            # Parity -1 = Start with INDENTATION (Slot).
            
            # Start of segment i
            p_start = (curr_x, curr_y)
            p_end_center = (x1 + (i+1)*actual_tab_w*ux, y1 + (i+1)*actual_tab_w*uy)
            
            # Vector for protrusion
            # if side == 1: displacement = +t * N
            # if side == -1: displacement = 0 (if we consider standard edge as base) OR -t * N?
            # Standard: Line is 0.
            # Tab goes to +t. Slot goes to -t? No, usually one part has tabs 0..t, other 0..t?
            # Mating: Part A (Tab) fills Part B (Slot).
            # If A has Tab (+t), B must have Slot (space).
            # Let's assume the "nominal edge" is the mating surface.
            # Tab: Cut path is at +t.
            # Slot: Cut path is at -t (or 0 if we subtract from stock?).
            # Let's use symmetric crossing: Tabs go +t/2, Slots go -t/2? No.
            # Robust: Line is geometric boundary.
            # Tab: Go Out t, traverse, go In.
            # Slot: Stay on line? And other part has Tabs?
            # If Base is S x S and Front sits on side...
            # Front needs SLOTS at bottom to accept Base TABS.
            # Base needs TABS at edges.
            
            # Logic:
            # Parity 1 (Tab): Line is at +thickness relative to (ux,uy) path?
            # Wait, easier loop:
            # We are at 'current level' (either 0 or t, or -t).
            # Transition draw perpendicular?
            
            # Let's stick to: Line is '0'.
            # Tab: Cut goes OUT (Left of vector) by 'thickness'.
            # Slot: Cut goes IN (Right of vector) by 'thickness' ? Or stays at 0?
            # If we want flush exterior:
            # Base S x S.
            # Front Bottom Edge: Needs to mate with Base Edge.
            # If Base has Tabs sticking OUT (S -> S+2t), Front Bottom needs Slots sticking UP (indenting H).
            # So "Slot" means Indent (Right of vector).
            # "Tab" means Outdent (Left of vector).
            
            indent = thickness
            
            # Determine displacement for this segment
            # side 1 -> Left (Out) -> disp = +indent * N
            # side -1 -> Right (In) -> disp = -indent * N (or 0?)
            # Usually strict box joint: One part goes 0..+t, other 0..-t?
            # Let's use:
            # Tab = displacement +t/2? No.
            # Let's assume "Line" is the mean.
            # Tab = Line + t (Left)
            # Slot = Line (Recessed? No, that's partial).
            # Slot = Line - t (Right) ? No, that eats into part.
            # Yes, Slot removes material. Tab adds material.
            # BUT: We need to respect the part's overall dimensions.
            # Part dimensions defined by the "zero line".
            # If Base is S x S. We want tabs to stick OUT of S. (Base total w ~ S+2t).
            # Then Base edge parity = 1 (Start Tab). "Slot" segments are at 0. "Tab" segments at +t.
            # Corresponding Wall: Edge parity = -1 (Start Slot). "Slot" segments at 0 (match Base Slot). "Tab" segments at -t (Indent to accept Base Tab).
            # Wait, if Base has Tab at +t, Wall needs Slot at +t? No, Wall needs space.
            # Wall has material. Slot REMOVES material.
            # So Wall Slot should go IN (Right of vector).
            # Does Wall start at '0' or 'S'?
            # If Wall sits ON SIDE of Base...
            
            # SIMPLIFIED MODEL:
            # Base (S x S) has TABS sticking OUT. (0 and +t). Parity = 1.
            # Wall (Bottom Edge) has SLOTS cutting IN. (0 and +t relative to wall bottom up?). 
            # Wall Bottom Edge Vector goes Right. Up is 'Left' (Into material). Down is 'Right' (Out).
            # We want slots to cut UP (Into material, Left). 
            # So Wall Slot = Left (+t). Wall Tab = 0.
            # This matches Base!
            # Base Tab (+t Out) meets Wall Slot (0? No).
            # Mismatch.
            
            # Let's use explicit 'Tab' vs 'Inverse Tab'.
            # Base Edge (Parity 1): Seg 1 is Tab (+t). Seg 2 is Base (0).
            # Wall Edge (Parity -1): Seg 1 is Slot (Indent, so cut at +t into material? No, cut at +t allows Tab to enter? No).
            # Wall mating edge:
            # Needs to be the NEGATIVE of Base Edge.
            # If Base Cut = Path P.
            # Wall Cut = Path P (reversed?).
            
            # Let's just calculate the offset for the segment:
            # Parity 1 (Protruding Tab Sequence): [Tab, Base, Tab, Base...]
            #   - Tab: Offset = +thickness (Left/Out)
            #   - Base: Offset = 0
            # Parity -1 (Recessed Slot Sequence): [Slot, Wall, Slot, Wall...]
            #   - Slot: Offset = 0 (Matches Base 'Base') ? No.
            #   - Wall: Offset = +thickness (Matches Base 'Tab'?) No.
            
            # Correct Logic:
            # Base Edge: [Tab (+t), Gap (0), Tab (+t)...]
            # Mating Wall Edge: [Slot (open), Finger (solid), Slot...]
            # Geometry must align.
            # If Base Edge is at Y=0. Tab goes Y=-t (Out). Gap is Y=0.
            # Wall Edge is at Y=0. Slot must be Y=0 (to accept Gap? No).
            # Slot must allow Tab (Y=-t) to enter.
            # So Wall edge must trace Y=-t and Y=0?
            # Yes! The cut path is IDENTICAL for both parts, just one is "Solid below", one is "Solid above".
            # So we just need to ensure the phases match.
            # If Base starts with Tab (+t), Wall must start with Slot (Space for Tab).
            # So Wall cut path matches Base cut path.
            
            # Draw Logic:
            # We draw the path relative to the line (x1,y1)->(x2,y2).
            # Parity 1: Starts at +t (Left). [ +t, 0, +t, 0 ]
            # Parity -1: Starts at 0. [ 0, +t, 0, +t ] (Complementary)
            
            # Note: "Left" of vector (x1->x2).
            # For Base Bottom edge (Right->Left): Left is Down (Out). Correct.
            # For Wall Bottom edge (Left->Right): Left is Up (In). 
            # If Wall Bottom has Parity -1 (Start 0).
            # Seq: 0, +t (In), 0, +t (In).
            # 0 aligns with Base Gap (0).
            # +t (In) aligns with Base Tab (+t Out). (Since Base Out = Down, Wall In = Up. Perfect overlap?)
            # Wait. Base Tab +t (Out/Down). Wall Slot +t (In/Up).
            # They occupy the same space? No. Solid meets Space.
            # Base is Solid "Above" cut. Wall is Solid "Above" cut.
            # If cut goes Down (Away from Base), Base gets Tab.
            # If cut goes Up (Into Wall), Wall gets Slot.
            # So we just need consistent "Left" offset logic.
            
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

    # 2. Front (S+2t x H)
    # Location: Below Base? No, we set layout earlier.
    fx, fy = parts['Front']['x'], parts['Front']['y']
    fw, fh = parts['Front']['w'], parts['Front']['h']
    
    # Bottom Edge (Mates with Base Top? No, Base Front side).
    # Base Front side (Bottom edge on screen) was R->L Parity 1.
    # To mate, Front Bottom Edge (L->R usually) needs Parity -1 (Start 0, then t, to match Base 1, 0...)
    # Base R->L (1, -1, 1...) -> Tab, Gap, Tab.
    # Front L->R must be: Gap, Tab, Gap. (Complementary).
    # Gap means Parity -1 (Start 0).
    # Wait, Front Wall Width is S+2t.
    # Base Width is S.
    # The CENTER S portion of Front Bottom mates with Base.
    # The CORNER t portions mate with Side Walls.
    # This is complex "Pinwheel" logic.
    
    # SIMPLIFIED LAYOUT:
    # All Walls sit ON TOP of Base (Base is largest S+2t x S+2t).
    # OR Base fits INSIDE.
    # Our dimensions: Base S x S. 
    # This means Base fits INSIDE.
    # Front Wall (S+2t) overlaps Left/Right Walls.
    
    # Let's fix the Mating Logic for "Base Fits Inside":
    # 1. Base: All edges TABS (Parity 1).
    # 2. Front/Back (Bottom): SLOTS (Parity -1). BUT Width is S + 2t.
    #    The slots for Base must be in the CENTER S region.
    #    The corners (size t) are solid?
    #    Actually, standard box joint generator:
    #    Base (S x S): Tabs.
    #    Front (S+2t x H): 
    #       - Bottom: Center S has Slots (to rx Base Tabs). Corners solid.
    #       - Sides: Tabs/Slots to rx Left/Right walls.
    # Let's handle the Bottom Edge of Front specifically.
    # Length S+2t.
    # 0..t: Corner (Solid).
    # t..S+t: Mating Zone (Slots).
    # S+t..S+2t: Corner (Solid).
    
    # Helper for Composite Edge?
    # Or just Draw Line, Draw Fingers, Draw Line.
    
    def draw_wall_bottom_with_base_slots(px, py, pw, ph):
        # Start (Bottom Left) -> t (Solid) -> S (Fingers) -> t (Solid) -> End
        p_y_bot = py + ph
        # 1. Left Corner
        svg.add_line(px, p_y_bot, px+t, p_y_bot, "CUT")
        # 2. Mating Zone (Parity -1 to Rx Base Tabs)
        draw_finger_edge(px+t, p_y_bot, px+t+S, p_y_bot, t, -1, "CUT")
        # 3. Right Corner
        svg.add_line(px+t+S, p_y_bot, px+t+S+t, p_y_bot, "CUT")


    # Front Wall Draw
    # Top: Flat
    svg.add_line(fx, fy, fx+fw, fy, "CUT")
    # Sides: Mating with Left/Right Walls.
    draw_finger_edge(fx+fw, fy, fx+fw, fy+fh, t, 1, "CUT") # Right Side (Down)
    draw_wall_bottom_with_base_slots(fx, fy, fw, fh)       # Bottom
    draw_finger_edge(fx, fy+fh, fx, fy, t, 1, "CUT")       # Left Side (Up)
    
    # Divider Slot (Vertical) on Front Wall
    # Divider runs Left-Right or Front-Back?
    # Usually Divider runs Left-Right (connecting Left and Right walls)? No, that's parallel to Front/Back.
    # Or Front-Back (connecting Front and Back)?
    # If Divider connects Front and Back, it is parallel to Left/Right.
    # Then Slots are needed on Front and Back.
    # Position is ratio along the Width (S).
    if params.get('divider_pos') is not None:
        # pos is 0..1 along S.
        # Front Wall Width is S+2t. The Inner S part is from fx+t to fx+t+S.
        pos_ratio = params['divider_pos']
        # X coord relative to Front Wall start
        # Slot center x
        slot_cx = fx + t + S * pos_ratio
        
        # Draw Slot (Width t, Height H/2 or similar?)
        # Let's do a simple vertical slot from Top down? Or internal?
        # Internal slot is stronger.
        # Slot height: e.g. H_inner - 20? 
        # Divider height is H-clearance.
        # Let's make a slot of height H/2 centered?
        # Or simpler: 2 slots (top and bottom)?
        # Let's just do one central vertical slot of length H/2.
        
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


    rx, ry = parts['Right']['x'], parts['Right']['y']
    rw, rh = parts['Right']['w'], parts['Right']['h']
    
    svg.add_line(rx, ry, rx+rw, ry, "CUT")
    draw_finger_edge(rx+rw, ry, rx+rw, ry+rh, t, -1, "CUT")
    draw_finger_edge(rx+rw, ry+rh, rx, ry+rh, t, -1, "CUT")
    draw_finger_edge(rx, ry+rh, rx, ry, t, -1, "CUT")

    svg.save()
    return True

def generate_divider_svg(params, filename="divider_square.svg"):
    """Generates a simple rectangular divider insert."""
    # Just a rectangle sized S-tolerance x H-tolerance?
    # Or specific slots?
    # Rubric: "a divider that inserts into the container".
    # For acrylic, loose insert is safest unless we add internal slots to walls.
    # Adding internal slots is complex (modify wall drawing).
    # Let's stick to simple rectangle for loose fit.
    
    S = params['S']
    H = params['H']
    t = params['t'] # not used for dim, but thick
    
    clearance = 0.5
    div_w = S - 2*t - clearance # Fits inside Inner S?
    # Wait, Inner Dim is S. Walls are outside base.
    # So Inner Box is S x S.
    # Divider Width = S - clearance.
    
    div_w = S - clearance
    div_h = H - clearance # Flush with top?
    
    stock_w = params['stock_w']
    stock_h = params['stock_h']
    
    svg = SVGGenerator(filename, stock_w, stock_h)
    svg.add_rect(stock_w/2 - div_w/2, stock_h/2 - div_h/2, div_w, div_h, "CUT")
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
        params['H'] = get_float("Box Height H (mm)", 60.0)
        
        print("\n-- Features --")
        do_div = get_bool("Include Divider Insert?")
        if do_div:
            params['divider_pos'] = get_float("Divider Position Ratio (0.0-1.0, 0.5=Center)", 0.5)
        
        print("\n-- Engraving --")
        # User request: "for the base I want it to say Digital Manufacturing"
        # We can set this as default or hardcode it. Let's provide it as default.
        default_base_text = "Digital Manufacturing"
        t_top = input(f"Text on Base [{default_base_text}]: ")
        params['text_top'] = t_top if t_top else default_base_text
        
        # User request: "for the 4 walls I want no text"
        # We will disable the user-specified front text prompt.
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
            if do_div:
                generate_divider_svg(params, "box_divider.svg")
            print_software_description(params)
            
    else:
        print("Invalid selection.")

if __name__ == "__main__":
    main()
