import math
BR = 4.1774

# Beesoon Figures 1-2, as published.  (Kd uM, n)
PFOS = {"linear": (0.08, 2.0), "3m": (400, 0.5), "4m": (80, 0.3),
        "5m": (90, 0.6)}
PFOA = {"linear": (100, 2.0), "3m": (400, 0.5), "5m": (300, 0.7),
        "4m": (None, None)}          # 4m-PFOA: p > 0.05, not estimable

# ECF manufacturing composition, from the paper's own Introduction
MIX = {"PFOS": 0.70, "PFOA": 0.80}   # linear fraction

def logK(n, kd):
    return BR + math.log10(n) - math.log10(kd)

print("At low occupancy each isomer contributes n_i/Kd_i to the bound")
print("pool, so a mixture's apparent n/Kd is SUM f_i (n_i/Kd_i) and")
print("log K_mix = 4.1774 + log10 of that sum.\n")

for name, iso in (("PFOS", PFOS), ("PFOA", PFOA)):
    fl = MIX[name]
    branched = [k for k in iso if k != "linear"]
    fb = (1 - fl) / len(branched)
    print(f"{name}: technical ECF is {fl*100:.0f}% linear, "
          f"{(1-fl)*100:.0f}% branched ({len(branched)} isomers "
          f"at {fb*100:.1f}% each)")
    print(f"   {'isomer':8s} {'f':>6s} {'Kd uM':>8s} {'n':>5s} "
          f"{'n/Kd':>12s} {'share of bound':>15s}")
    contrib = {}
    for k in ["linear"] + branched:
        kd, n = iso[k]
        f = fl if k == "linear" else fb
        c = 0.0 if kd is None else f * n / kd
        contrib[k] = c
        kds = "n.e." if kd is None else f"{kd:.4g}"
        ns = "n.e." if n is None else f"{n:.1f}"
        print(f"   {k:8s} {f:6.2f} {kds:>8s} {ns:>5s} "
              f"{(0 if kd is None else n/kd):12.5g}", end="")
        print(f" {'':>15s}" if kd is None else "", end="")
        print()
    tot = sum(contrib.values())
    for k in contrib:
        print(f"      {k:8s} contributes {contrib[k]/tot*100:6.2f}% "
              f"of the bound pool")
    lk_mix = BR + math.log10(tot)
    lk_lin = logK(*reversed(list(reversed(iso['linear']))[::-1])) \
        if False else logK(iso["linear"][1], iso["linear"][0])
    print(f"   log K, pure linear          {lk_lin:.3f}")
    print(f"   log K, technical mixture    {lk_mix:.3f}")
    print(f"   DIFFERENCE                  {lk_mix-lk_lin:+.3f} log\n")

print("Context: the sheet's pooled SD inside one compound is 0.87 log.")
print("Both differences are far below it, so whether another lab used")
print("pure linear or technical material cannot be read off its log K.")
print()
print("Note this conclusion is INDEPENDENT of the mM/M question: every")
print("Kd scales by the same 1000x, so n/Kd ratios and therefore the")
print("mixture-vs-linear gap are unchanged.")
