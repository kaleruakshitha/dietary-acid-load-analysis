"""Feature-engineering formulas used by the analysis pipeline."""
from __future__ import annotations
import math

def _require_nonnegative(name: str, value: float) -> None:
    if not math.isfinite(value): raise ValueError(f"{name} must be finite")
    if value < 0: raise ValueError(f"{name} cannot be negative")

def calculate_pral(protein_g: float, phosphorus_mg: float, potassium_mg: float, magnesium_mg: float, calcium_mg: float) -> float:
    """Calculate potential renal acid load in mEq/day; protein g, minerals mg."""
    inputs = {"protein_g":float(protein_g),"phosphorus_mg":float(phosphorus_mg),"potassium_mg":float(potassium_mg),"magnesium_mg":float(magnesium_mg),"calcium_mg":float(calcium_mg)}
    for name,value in inputs.items(): _require_nonnegative(name,value)
    return (0.49*inputs["protein_g"]+0.037*inputs["phosphorus_mg"]-0.021*inputs["potassium_mg"]-0.026*inputs["magnesium_mg"]-0.013*inputs["calcium_mg"])

def calculate_egfr_2021(serum_creatinine_mg_dl:float,age_years:float,sex:str)->float:
    """Race-free 2021 CKD-EPI creatinine eGFR (mL/min/1.73 m²)."""
    scr=float(serum_creatinine_mg_dl);age=float(age_years)
    if not math.isfinite(scr) or scr<=0: raise ValueError("serum_creatinine_mg_dl must be finite and positive")
    if not math.isfinite(age) or age<18: raise ValueError("age_years must be finite and at least 18")
    normalized_sex=str(sex).strip().lower()
    if normalized_sex in {"female","f"}: kappa=0.7;alpha=-0.241;sex_factor=1.012
    elif normalized_sex in {"male","m"}: kappa=0.9;alpha=-0.302;sex_factor=1.0
    else: raise ValueError("sex must be female/f or male/m for this equation")
    ratio=scr/kappa
    return 142.0*min(ratio,1.0)**alpha*max(ratio,1.0)**-1.200*0.9938**age*sex_factor
