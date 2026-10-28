import math

def logK(Kd_uM, M_gmol, n=1.0, Cw_uM=0.0, per="kg", rho_kg_per_L=1/0.733):
    """log10 of the protein-water partition coefficient.

    Kd_uM   dissociation constant, micromolar
    M_gmol  molar mass of the protein actually used in the assay
    n       binding sites per protein molecule (from the SAME paper)
    Cw_uM   free water concentration; 0 = the trace-limit plateau
    per     "kg" or "L" of protein
    """
    C = math.log10(1e9 / M_gmol)              # umol of protein per kg of it
    if per == "L":
        C = math.log10(1e9 * rho_kg_per_L / M_gmol)
    return C + math.log10(n) - math.log10(Kd_uM) - math.log10(1 + Cw_uM / Kd_uM)

def Kd_from_logK(lK, M_gmol, n=1.0, per="kg", rho_kg_per_L=1/0.733):
    """The same equation solved the other way. Needs the same n."""
    C = math.log10(1e9 / M_gmol)
    if per == "L":
        C = math.log10(1e9 * rho_kg_per_L / M_gmol)
    return 10 ** (C + math.log10(n) - lK)

HSA, LFABP = 66463.0, 14208.0

print("constant log10(1e9/M):")
print(f"  albumin  {math.log10(1e9/HSA):.4f}   L-FABP {math.log10(1e9/LFABP):.4f}"
      f"   gap {math.log10(HSA/LFABP):+.3f} log = {HSA/LFABP:.2f}x")

print("\nChen 2025, PFOA on HSA, Kd 2.57 uM, n = 1")
print(f"  log K = {logK(2.57, HSA):.2f} L_w/kg      log D = {logK(2.57, HSA, per='L'):.2f} L_w/L")

print("\nKhazaee 2021, PFOA on human L-FABP, Kd 0.099 uM, single-site fit")
print(f"  log K = {logK(0.099, LFABP):.2f} L_w/kg")
print(f"  with albumin's mass by mistake: {logK(0.099, HSA):.2f}  <- {HSA/LFABP:.1f}x low")

print("\nround trip")
print(f"  {Kd_from_logK(logK(0.099, LFABP), LFABP):.3f} uM back from log K")

print("\nsaturation: the same Kd at three free concentrations (L-FABP, PFOA)")
for Cw in (0.0, 0.099, 1.0):
    print(f"  Cw = {Cw:5.3f} uM -> log K = {logK(0.099, LFABP, Cw_uM=Cw):.2f}")

print("\nseveral classes of site: Woodcroft 2010, PFOA on rat L-FABP")
for label, kds in [("displacement", [13.1, 23.8]), ("ITC", [3.1, 26.2, 52.6])]:
    K0 = (1e9 / 14272.0) * sum(1.0 / k for k in kds)
    print(f"  {label:13s} {len(kds)} sites -> log K0 = {math.log10(K0):.2f}")
