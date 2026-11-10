import math
from dataclasses import dataclass

@dataclass
class Rect:
    name: str
    b: float      # width (horizontal)
    h: float      # height (vertical)
    yc: float     # centroid height measured from bottom of beam
    n: float = 1.0  # modular ratio relative to the beam material
 
    @property
    def area(self):            # transformed area
        return self.n * self.b * self.h
 
    @property
    def Ix_own(self):          # about own horizontal centroidal axis
        return self.n * self.b * self.h ** 3 / 12
 
    @property
    def Iy_own(self):          # about vertical axis of symmetry
        return self.n * self.h * self.b ** 3 / 12

def section_properties(d, bf, tf, tw, bp, tp, n=1.0):
    """
        d  : overall depth of the I-beam
       bf : flange width,  tf : flange thickness,  tw : web thickness
       bp : plate width,   tp : plate thickness
       n  : modular ratio E_plate / E_beam (1.0 for same material)
    """
    hw = d - 2 * tf
    if hw <= 0:
        raise ValueError("Depth must be larger than twice the flange thickness.")
 
    parts = [
        Rect("Bottom flange", bf, tf, tf / 2),
        Rect("Web",           tw, hw, d / 2),
        Rect("Top flange",    bf, tf, d - tf / 2),
        Rect("Plate",         bp, tp, d + tp / 2, n),
    ]
 
    A = sum(p.area for p in parts)
    y_bar = sum(p.area * p.yc for p in parts) / A          # from bottom
    H = d + tp                                              # total height
 
    # Parallel axis theorem
    Ix = sum(p.Ix_own + p.area * (p.yc - y_bar) ** 2 for p in parts)
    Iy = sum(p.Iy_own for p in parts)

    c_top = H - y_bar
    c_bot = y_bar
 
    res = {
        "A": A,
        "y_bar_from_bottom": y_bar,
        "y_bar_from_top": c_top,
        "total_height": H,
        "Ix": Ix,
        "Iy": Iy,
    }

    return res


# Scale model: FD7TA 

# main beam (W12X14) + 0.375 Checker Plate

W12X14 = section_properties(302.5e-3, 100.8e-3, 5.72e-3, 5.08e-3, 100.8e-3, 9.525e-3)

# Scale model: FD5TA 

# main beam (W12X22) + 0.5 Checker Plate

W12X22 = section_properties(312.00e-3, 102e-3, 10.8e-3, 6.6e-3, 102e-3, 12.7e-3)

# Deck Stringer (W4X13) + 0.5 Checker Plate

W4X13 = section_properties(106e-3, 103e-3, 8.76e-3, 7.11e-3, 103e-3, 12.7e-3)

def main():
    W12X14 = section_properties(302.5e-3, 100.8e-3, 5.72e-3, 5.08e-3, 0,0)
    
    W12X22 = section_properties(312.00e-3, 102e-3, 10.8e-3, 6.6e-3, 102e-3, 12.7e-3)
    
    W4X13 = section_properties(106e-3, 103e-3, 8.76e-3, 7.11e-3, 103e-3, 12.7e-3)

    print(f'W12X22: Ix = {W12X22['Ix']:.2e} m4\nW4X13 : Ix = {W4X13['Ix']:.2e} m4')

    Itotal = 6*W12X22['Ix'] + 4*W4X13['Ix']

    print(f'W12X22/Itotal(%) = {W12X22['Ix']/Itotal*100:.2f}%')

    #print(f'Beam area FD7: {11*W12X14['A']: .3f} m2\nBeam area FD5: {6*W12X22['A'] + 3*W4X13['A']: .3f} m2')
    

if __name__ == "__main__":
    main()