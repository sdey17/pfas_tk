"""How albumin's share of binding is computed, from eq 2 up."""

# ---------- STEP 0.  The one idea ------------------------------------
# Qin's eq 2:   n_tot = n_free + SUM_i n_bound,i
#                     = C_free x ( V_w + SUM_i D_i/w x V_i )
#
# Each phase i contributes ONE product to that sum:  D_i/w x V_i.
# Units: [L_water / L_phase] x [L_phase] = L_water.  So D x V is the
# "water-equivalent volume" of that phase -- the volume of water that
# would hold as much chemical as the phase does.
#
# Because every term shares the same C_free, the MOLES on phase i are
#       n_bound,i = C_free x D_i/w x V_i
# and the share is just that product over the sum of the products.
# C_free cancels, so the share does not depend on dose.

# ---------- STEP 1.  The volumes (one 384-well) ----------------------
V_MED, V_CELL = 4.00e-5, 1.71e-8          # L   SI Tables S7, S8
V = {
    "water, medium":   0.999    * V_MED,   # Vf_w,medium   = 0.999 L/L
    "albumin (medium protein)": 0.94e-3  * V_MED,   # 0.94   mL/L
    "lipid, medium":   1.47e-5  * V_MED,   # 0.0147 mL/L
    "protein, cell":   33.9e-3  * V_CELL,  # 33.9   mL/L of cell
    "lipid, cell":     5.3e-3   * V_CELL,  # 5.3    mL/L of cell
}
print("STEP 1.  phase volumes in one well")
for k, v in V.items():
    print(f"   {k:26s} {v:.4g} L")
print(f"\n   albumin volume / medium-lipid volume  = "
      f"{V['albumin (medium protein)']/V['lipid, medium']:.0f}x")
print(f"   albumin volume / cell-protein volume  = "
      f"{V['albumin (medium protein)']/V['protein, cell']:.0f}x")
print("   ^ this is why albumin wins: VOLUME, not affinity.\n")

# ---------- STEP 2.  Worked example, PFOA ----------------------------
# Table 1, PFOA row:  log D_BSA/w 3.98, log D_SP/w 3.14, log D_lip/w 3.52
D = {"albumin (medium protein)": 10**3.98, "lipid, medium": 10**3.52,
     "protein, cell": 10**3.14, "lipid, cell": 10**3.52}

print("STEP 2.  PFOA, the four D x V products")
print(f"   {'phase':26s} {'D [Lw/Lph]':>11s} {'V [L]':>10s} "
      f"{'D x V [Lw]':>11s} {'share':>7s}")
terms = {k: D[k] * V[k] for k in D}
tot = sum(terms.values())
for k in D:
    print(f"   {k:26s} {D[k]:11.4g} {V[k]:10.4g} "
          f"{terms[k]:11.4g} {terms[k]/tot*100:6.1f}%")
print(f"   {'':26s} {'':>11s} {'sum bound':>10s} {tot:11.4g} "
      f"{'100.0%':>7s}")
print(f"\n   albumin share = {terms['albumin (medium protein)']:.4g}"
      f" / {tot:.4g} = "
      f"{terms['albumin (medium protein)']/tot*100:.2f}%")
print("   The free term V_w is NOT in the denominator: this is a share")
print("   of what is BOUND, not of what was dosed.\n")

# ---------- STEP 3.  The independent cross-check ---------------------
# Delete the three non-albumin phases and see how far f_free moves.
Vw, Vtot = V["water, medium"], V_MED + V_CELL
f_sys = Vtot / (Vw + tot)
f_alb = Vtot / (Vw + terms["albumin (medium protein)"])
print("STEP 3.  cross-check: delete the other three phases entirely")
print(f"   f_free, all four phases      {f_sys*100:7.3f}%")
print(f"   f_free, albumin term alone   {f_alb*100:7.3f}%")
print(f"   difference                   {(f_alb-f_sys)*100:+7.3f} "
      "percentage points")
print("   Two routes to the same conclusion: the share says 99.2% of")
print("   bound moles are on albumin; the deletion says dropping the")
print("   rest moves the free fraction by 0.08 pp.\n")

# ---------- STEP 4.  Why the share ranges 92.8-99.2% -----------------
print("STEP 4.  the spread is set by D_lip/w and D_SP/w vs D_BSA/w")
T1 = {  # (logD_BSA/w, logD_SP/w, logD_lip/w)
 "PFOA":(3.98,3.14,3.52), "PFUnA":(4.75,3.97,4.54),
 "PFHxS":(3.52,2.95,4.13), "PFOSA":(4.33,3.52,4.94),
 "HFPO-DA":(2.28,2.77,2.41), "PFOS":(4.52,3.94,4.89),
}
print(f"   {'PFAS':9s} {'D_lip/D_BSA':>12s} {'D_SP/D_BSA':>11s} "
      f"{'share':>7s}  biggest competitor")
for c,(lb,ls,ll) in sorted(T1.items(), key=lambda x:-x[1][0]+x[1][2]):
    Dc = {"alb":10**lb, "lm":10**ll, "pc":10**ls, "lc":10**ll}
    tm = {"alb":Dc["alb"]*V["albumin (medium protein)"],
          "lm":Dc["lm"]*V["lipid, medium"],
          "pc":Dc["pc"]*V["protein, cell"],
          "lc":Dc["lc"]*V["lipid, cell"]}
    s = sum(tm.values())
    comp = max(("lm","pc","lc"), key=lambda k: tm[k])
    name = {"lm":"medium lipid","pc":"cell protein","lc":"cell lipid"}[comp]
    print(f"   {c:9s} {10**(ll-lb):12.2f} {10**(ls-lb):11.2f} "
          f"{tm['alb']/s*100:6.1f}%  {name}")
print("\n   Lipid can out-affinity albumin 4x (PFHxS, PFOSA) and still")
print("   lose, because it is 1/64 the volume.  HFPO-DA is the one case")
print("   where CELL PROTEIN is the competitor, not lipid.")
