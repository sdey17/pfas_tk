"""
How much of the mouse-to-human renal clearance gap is protein binding?

Section 6 adjudicates what PFAS binding to albumin actually is; section 10 shows
what the choice of unbound fraction does to the transport endpoint R. Neither
states the arithmetic that connects them, which is the thing a reader asking
"does binding affect elimination?" wants: binding is one multiplicative term in
the renal gap, and under the measured unbound fraction it is the LARGEST one.

Only unbound PFAS is filtered at the glomerulus, so the clearance free
filtration alone would produce is fu x GFR. Measured renal clearance divides
into three factors:

    CL_renal_mouse     fu_mouse     GFR_mouse         R_mouse
    --------------  =  --------  x  ---------    x   --------
    CL_renal_human     fu_human     GFR_human         R_human

                        BINDING      BODY SIZE        TRANSPORT

the last being the dimensionless handling ratio R = CL_renal/(fu x GFR) of
section 10, which is what a structure-activity model would have to learn.

Run:  python3 scripts/binding_elimination_decomposition.py
"""

import csv
import math
import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

GFR_MOUSE = 10_000.0   # mL/d/kg, back-calculated from Argoul 2026 Table 4
GFR_HUMAN = 2_570.0    # mL/d/kg, Han 2012 Table 4 (125 mL/min / 70 kg)
CL_RENAL_HUMAN = 0.03  # mL/d/kg, Han 2012 Table 4

# The three human unbound fractions this review found in circulation.
FU_HUMAN = [
    (0.00061, "Fischer, measured at physiological ligand:protein"),
    (0.02,    "Han 2012's stated assumption"),
    (0.10,    "PBPK models reading '>90% bound' as 'about 90%'"),
]


def main():
    p = os.path.join(HERE, "db", "qsar", "qsar_endpoint_table.csv")
    a = {r["chemical"]: r for r in csv.DictReader(open(p))}["PFOA"]
    cl_m = float(a["cl_renal_mL_kg_day"])
    fu_m = float(a["fu_pct"]) / 100

    print(__doc__.split("Run:")[0].rstrip())
    print("=" * 78)
    print("PFOA, female mouse against human")
    print("=" * 78)
    print(f"""
  mouse  CL_renal {cl_m:>6.2f} mL/d/kg   fu {fu_m:.5f}   GFR {GFR_MOUSE:>7,.0f} mL/d/kg
  human  CL_renal {CL_RENAL_HUMAN:>6.2f} mL/d/kg   fu (three values below)
                                         GFR {GFR_HUMAN:>7,.0f} mL/d/kg

  The renal clearance gap is fixed at {cl_m/CL_RENAL_HUMAN:.0f}x whatever fu is assumed.
  What the assumption changes is how that gap is APPORTIONED.
""")
    print(f"  {'human fu':>9s} {'binding':>9s} {'body size':>10s} {'transport':>10s}"
          f"   dominant term")
    for fu_h, label in FU_HUMAN:
        g_cl = cl_m / CL_RENAL_HUMAN
        g_fu, g_gfr = fu_m / fu_h, GFR_MOUSE / GFR_HUMAN
        g_tr = g_cl / (g_fu * g_gfr)
        terms = {"binding": g_fu, "body size": g_gfr, "transport": g_tr}
        top = max(terms, key=terms.get)
        print(f"  {fu_h:>9.5f} {g_fu:>8.1f}x {g_gfr:>9.1f}x {g_tr:>9.1f}x   {top}")
    print()
    for fu_h, label in FU_HUMAN:
        print(f"    fu {fu_h:<8g} {label}")

    print("""
  Read down the binding column. Under the MEASURED unbound fraction, binding is
  the largest single term in the species difference: humans hold PFOA about 14x
  more tightly than mice do, which alone keeps 14x less of it in front of the
  glomerulus. Body size adds 3.9x. What is left for transport biology is 2.2x.

  Under the ASSUMED value the apportionment inverts. Binding becomes a term
  that works the WRONG WAY (0.4x or 0.1x, i.e. humans would bind LESS tightly
  than mice), and transport has to carry 71x or 354x to make the arithmetic
  balance. That is not a small difference of emphasis: it decides whether the
  species difference is a binding problem or a transporter problem, and
  therefore which experiment is worth doing next.

  Two cautions before this is used. The mouse and human unbound fractions come
  from different methods -- equilibrium dialysis in Argoul, solid-phase
  microextraction in Fischer -- at different ligand:protein ratios, which is
  precisely the artefact section 6.1 diagnoses, so the 14x is the number least
  secure here. And the mouse limb inherits the cocktail study's roughly 2.5x low
  clearance bias (section 11.3a), which lands on the transport term: correcting
  it would raise the 2.2x to about 5x, still far below 71x.
""")

    print("=" * 78)
    print("Why binding alone does not explain long human half-lives")
    print("=" * 78)
    fu_h = FU_HUMAN[0][0]
    filt = fu_h * GFR_HUMAN
    print(f"""
  A common shortcut says PFAS persist because they are bound to albumin and so
  cannot be filtered. Binding is necessary to that story but nowhere near
  sufficient, and the numbers say why.

  At fu = {fu_h}, free filtration alone would clear {filt:.2f} mL/d/kg of human plasma.
  With Vd = 170 mL/kg that is a half-life of {math.log(2)*170/filt:.2f} days -- weeks, not years.
  The measured human renal clearance is {CL_RENAL_HUMAN} mL/d/kg, {filt/CL_RENAL_HUMAN:.0f}x LOWER still.

  So binding sets how much PFOA is presented to the glomerulus, and tubular
  REABSORPTION then recovers about {100*(1-CL_RENAL_HUMAN/filt):.1f}% of what is filtered. Binding is the
  smaller of the two factors in absolute terms, and it is the one that is
  species-variable; reabsorption is the one that is near-complete in every
  species measured. A model that gets binding right and reabsorption wrong will
  still miss the half-life by two orders of magnitude.
""")


if __name__ == "__main__":
    main()
