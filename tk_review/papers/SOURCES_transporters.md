# SOURCES — PFAS renal and hepatic membrane transporters (in vitro kinetics)

Retrieval log for the transporter-kinetics strand. Date of all attempts: **2026-10-01**.
Agent: transporter in-vitro kinetics researcher. Output database: `db/transporter_kinetics.csv`.

Notes on this log
- `papers/` is shared with other agents in this review; where a file with the same content was
  already present under a different name, that pre-existing file is named here and no duplicate
  was created.
- Several Toxicological Sciences (Oxford) PMC records refuse XML/PDF download; those were
  retrieved as rendered PMC HTML via Firecrawl and saved as `.txt` (no PDF obtainable).
- `pmc.ncbi.nlm.nih.gov/articles/<PMCID>/pdf/` returns an anti-bot HTML interstitial for every
  PMC record tried from this container, so **no PMC PDFs could be downloaded**. Europe PMC
  `fullTextXML` works only for records in the OA subset.
- `n_perfluorinated_carbons` convention used in the CSV: a PFCA with *n* total carbons has
  *n*−1 perfluorinated carbons (PFOA = 7); a PFSA with *n* carbons has *n* (PFOS = 8). This is
  the same convention as the η_pfc used by Ryu et al. 2024.

---

## Successful retrievals (full text or full tables obtained)

| File in papers/ | Full citation | PMID | DOI | URL used | Result |
|---|---|---|---|---|---|
| `Weaver 2010 rat renal OAT PFCA chain length.txt` (pre-existing); also retrieved independently | Weaver YM, Ehresman DJ, Butenhoff JL, Hagenbuch B (2010). Roles of rat renal organic anion transporters in transporting perfluorinated carboxylates with different chain lengths. *Toxicol Sci* 113(2):305–314. | 19915082 | 10.1093/toxsci/kfp275 | https://pmc.ncbi.nlm.nih.gov/articles/PMC2807038/ (Firecrawl) | SUCCESS — full text + Tables 1 and 2 + Yang 2009 Ki values quoted in Discussion. Europe PMC fullTextXML and NCBI efetch both REFUSED ("publisher does not allow downloading of the full text in XML form"); PDF blocked. |
| `Zhao 2015 NTCP ASBT perfluoroalkyl sulfonates.txt` | Zhao W, Zitzow JD, Ehresman DJ, Chang SC, Butenhoff JL, Forster J, Hagenbuch B (2015). Na+/Taurocholate cotransporting polypeptide and apical sodium-dependent bile acid transporter are involved in the disposition of perfluoroalkyl sulfonates in humans and rats. *Toxicol Sci* 146(2):363–373. | 26001962 | 10.1093/toxsci/kfv102 | https://pmc.ncbi.nlm.nih.gov/articles/PMC4607751/ (Firecrawl) | SUCCESS — full text + Table 2 (hNTCP/rNTCP kinetics) + MRP2/BCRP/BSEP inhibition + ASBT and OSTα/β results. XML/PDF blocked. |
| `Zhao 2017 OATP disposition PFAA human rat.txt` (pre-existing) | Zhao W, Zitzow JD, Weaver Y, Ehresman DJ, Chang SC, Butenhoff JL, Hagenbuch B (2017). Organic anion transporting polypeptides contribute to the disposition of perfluoroalkyl acids in humans and rats. *Toxicol Sci* 156(1):84–95. | 28013215 | 10.1093/toxsci/kfw236 | https://pmc.ncbi.nlm.nih.gov/articles/PMC6075085/ (Firecrawl targeted query) | SUCCESS — Table 3 (human OATP1B1/1B3/2B1 Km, Vmax, Vmax/Km) and Table 4 (rat Oatp1a1/1a5) quoted verbatim; rOatp1b2/rOatp2b1 qualitative results. XML/PDF blocked. |
| `Louisse 2023 OAT4 PFAS substrates.txt` (pre-existing) | Louisse J, Dellafiora L, van den Heuvel JJMW, Rijkers D, Leenders L, Dorne JCM, Punt A, Russel FGM, Koenderink JB (2023). Perfluoroalkyl substances (PFASs) are substrates of the renal human organic anion transporter 4 (OAT4). *Arch Toxicol* 97(3). | 36436016 | 10.1007/s00204-022-03428-6 | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9968691/fullTextXML | SUCCESS — OA; full text and Table 1 (Km/Vmax/efficiency for 6 PFAS) and Table 2 (docking scores). Springer OA PDF download returned an HTML stub (FAIL for PDF). |
| `Ryu 2024 14 PFAS OAT interactions.txt` (pre-existing) | Ryu S, Yamaguchi E, Sadegh Modaresi SM, et al. (2024). Evaluation of 14 PFAS for permeability and organic anion transporter interactions: implications for renal clearance in humans. *Chemosphere* 361:142390. | 38801906 | 10.1016/j.chemosphere.2024.142390 | NCBI efetch db=pmc id=11774580 | SUCCESS — full text + Table 1 (Papp, ECCS class, OAT1–4/OCT2 substrate calls) + narrative uptake-ratio fold values. Note: the user's citation list named this paper "Niu Z"; the first author is **Ryu S** and the senior author is Slitt/Bello group, not Niu. |
| `Sakolish 2026 in vitro in silico renal clearance IVIVE.txt` | Sakolish C, Lin HC, Moyer HL, et al. (2026). An in vitro–in silico workflow for predicting renal clearance of environmental chemicals and drugs. *Toxicology* (online 2025). | 41224188 | 10.1016/j.tox.2025.154336 | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12951621/fullTextXML | SUCCESS — OA; Table 1 (28 PFAS tested) and Table 2 (relative renal clearance factors and predicted human half-life categories). |
| `Lin 2024 in vitro in silico renal clearance PFAS.txt` (pre-existing) | Lin HC, Sakolish C, Moyer HL, Carmichael PL, Baltazar MT, Ferguson SS, Stanko JP, Hewitt P, Rusyn I, Chiu WA (2024). An in vitro–in silico workflow for predicting renal clearance of PFAS. *Toxicol Appl Pharmacol* 491:117015. | 38917890 | 10.1016/j.taap.2024.117015 | PMC11585971 (pre-existing file) | SUCCESS — IVIVE performance statements extracted. |
| `Worley 2015 PBPK kidney transporters PFOA rat.txt` (pre-existing) | Worley RR, Fisher J (2015). Application of physiologically-based pharmacokinetic modeling to explore the role of kidney transporters in renal reabsorption of perfluorooctanoic acid in the rat. *Toxicol Appl Pharmacol* 289(3):428–441. | 26522833 | 10.1016/j.taap.2015.10.017 | NCBI efetch db=pmc id=4662604 | SUCCESS — Tables 2–6, including Km_baso = 27.20 µg/mL, Km_apical = 52.3 µg/mL and the fitted male/female relative activity factors. |
| `Kimura 2020 PFCA inhibition BSP uptake OATP Caco-2.txt` | Kimura O, Fujii Y, Haraguchi K, Kato Y, Ohta C, Koga N, Endo T (2020). Effects of perfluoroalkyl carboxylic acids on the uptake of sulfobromophthalein via organic anion transporting polypeptides in human intestinal Caco-2 cells. *Biochem Biophys Rep* 24:100807. | 32964147 | 10.1016/j.bbrep.2020.100807 | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7490525/fullTextXML | SUCCESS — OA; Table 1 and Dixon-plot Ki values (PFOA 62.2, PFNA 35.3, PFDA 43.2 µM). |
| `Niu 2026 three renal transporters PFAS.txt` (pre-existing) | Niu S, Ma Y, Tewari A, Ng C (2026). Transport of perfluoroalkyl substances (PFAS) by three renal transporters: implications for PFAS bioaccumulation mechanisms. *Environment & Health* 4(9):2150–2157. | 42775203 | 10.1021/envhealth.6c00088 | PMC13595387 (pre-existing file) | SUCCESS (text) / PARTIAL (numbers) — abstract and Results text give Km ranges (BCRP 1.29–7.92 µmol/L; P-gp 1.41–2.78 µmol/L) and all Vmax/Km values, but the per-congener Km and Vmax are presented **only inside Figures 1f and 2b,d**, so individual Km values could not be tabulated. |
| `Niu 2023 PFAS renal transporters review.txt` (pre-existing) | Niu S, Cao Y, Chen R, Bedi M, Sanders AP, Ducatman A, Ng C (2023). A state-of-the-science review of interactions of PFAS with renal transporters in health and disease: implications for population variability in PFAS toxicokinetics. *Environ Health Perspect* 131(7). | 37418334 | 10.1289/EHP11885 | PMC10328216 (pre-existing file) | SUCCESS — review; used for context only, all kinetic numbers in the CSV were taken from the primary papers. |
| `Furukawa 2024 ABCB1 ABCB4 PFOA mouse.txt` (pre-existing) | Furukawa K, Okamoto-Matsuda K, et al. (2024). Role of ABCB1 and ABCB4 in renal and biliary excretion of perfluorooctanoic acid in mice. *Environ Health Prev Med* 29. | 38522903 | 10.1265/ehpm.23-00284 | Europe PMC fullTextXML (PMC10965415) | SUCCESS — OA. |
| `Eldasher 2013 hepatic renal Bcrp expression PFOA mice.txt` | Eldasher LM, Wen X, Little MS, et al. (2013). Hepatic and renal Bcrp transporter expression in mice treated with perfluorooctanoic acid. *Toxicology* 306:108–113. | 23435180 | 10.1016/j.tox.2013.02.009 | NCBI efetch db=pmc id=3645321 | SUCCESS (XML) — expression study, not a kinetic study; not used in the CSV. |
| `ATSDR 2021 Toxicological Profile Perfluoroalkyls.txt` (pre-existing) | ATSDR (2021). Toxicological Profile for Perfluoroalkyls. | — | — | pre-existing in papers/ | SUCCESS — Section 3 transporter narrative, NHANES serum geometric means. |
| `OEHHA 2021 PFOA PFOS PHG first draft.txt`, `OEHHA 2024 PFOA PFOS PHG.txt` (pre-existing) | OEHHA (2021, 2024). Public Health Goals for PFOA and PFOS in Drinking Water. | — | — | pre-existing in papers/ | SUCCESS — **Appendix 6 Table A6.3** "Experimental evidence for PFOA and PFOS transporters" is the only accessible source for the Katakura 2007 results; Table A6.4 net tubular reabsorption by species. |
| `EPA 2023 IRIS PFHxA toxicological review.txt`, `EPA 2024 PFOA human health toxicity assessment.txt` (pre-existing) | US EPA (2023, 2024). | — | — | pre-existing in papers/ | SUCCESS — secondary statements on Oatp1a1 male/female fold differences and Zhao 2017 PFOA/PFNA/PFHpA OATP substrate calls. |

## Abstract-level retrievals (abstract contains the kinetic numbers; full text not obtained)

All retrieved 2026-10-01 from PubMed via `mcp__PubMed__get_article_metadata`. No PDF obtainable
for any of these (all paywalled, none in PMC).

| Citation | PMID | DOI | What was extracted | Result |
|---|---|---|---|---|
| Yang CH, Glover KP, Han X (2009). Organic anion transporting polypeptide (Oatp) 1a1-mediated perfluorooctanoate transport and evidence for a renal reabsorption mechanism of Oatp1a1 in renal elimination of perfluorocarboxylates in rats. *Toxicol Lett* 190(2):163–171. | 19616083 | 10.1016/j.toxlet.2009.07.011 | Oatp1a1–PFOA Km = 162.2 ± 20.2 µM (CHO); Ki,app values C6–C10; no inhibition by C4/C5 at 1 mM; in vitro–in vivo correlation with male-rat total clearance | PARTIAL SUCCESS (abstract). Full text FAILED — ScienceDirect paywall (a Firecrawl scrape of a guessed PII returned the wrong article). The five Ki,app values were recovered verbatim from the Weaver 2010 Discussion (PMC2807038). |
| Yang CH, Glover KP, Han X (2010). Characterization of cellular uptake of perfluorooctanoate via OATP1A2, OAT4 and URAT1 … *Toxicol Sci* 117(2):294–302. | 20639259 | 10.1093/toxsci/kfq219 | OAT4–PFOA Km 172.3 ± 45.9 µM (pH 6) and 310.3 ± 30.2 µM (pH 7.4); URAT1–PFOA Km 64.1 ± 30.5 µM (no extracellular Cl⁻); OATP1A2 does NOT transport PFOA; chain-length-dependent inhibition C4–C12 | PARTIAL SUCCESS (abstract). Full text FAILED — Oxford Academic paywall (Firecrawl returned metadata only, `citation_fulltext_world_readable` empty). Vmax values (OAT4 37, URAT1 0.32 nmol/min/mg protein; 24-fold efficiency difference) taken from the Louisse 2023 Introduction as a secondary quotation. |
| Nakagawa H, Hirata T, Terada T, Jutabha P, Miura D, Harada KH, Inoue K, Anzai N, Endou H, Inui K, Kanai Y, Koizumi A (2008). Roles of organic anion transporters in the renal excretion of perfluorooctanoic acid. *Basic Clin Pharmacol Toxicol* 103(1):1–8. | 18373647 | 10.1111/j.1742-7843.2007.00155.x | hOAT1 48.0 ± 6.4; rOAT1 51.0 ± 12.0; hOAT3 49.1 ± 21.4; rOAT3 80.2 ± 17.8 µM; hOAT2 and rOAT2 no transport | **FULL TEXT SUPPLIED BY THE REVIEWER 2026-11-19.** Table 1 adds the four Vmax values (479.3 / 421.3 / 399.1 / 365.6 pmol/mg protein/min), ingested by `scripts/add_nakagawa2008_oat_vmax.py`. **Table 1's third column is mislabelled** — headed "Km/Vmax (ul/mg protein/min)", but those units belong to Vmax/Km; the true clearances are 9.99 / 8.26 / 8.13 / 4.56 uL/min/mg and the ranking inverts. Human/rat clearance ratio 1.21 (OAT1) and 1.78 (OAT3) in one expression system — a key aim-1 result. |
| Nakagawa H, Terada T, Harada KH, Hitomi T, Inoue K, Inui K, Koizumi A (2009). Human organic anion transporter hOAT4 is a transporter of perfluorooctanoic acid. *Basic Clin Pharmacol Toxicol* 105(2):136–138. | 19371258 | 10.1111/j.1742-7843.2009.00409.x | hOAT4 transports PFOA; **no kinetic parameters reported** | PARTIAL FAILURE — PubMed has no abstract for this short communication; no PMC record; Wiley paywall. Characterisation ("active; uptake reported but no kinetic characteristics") taken from OEHHA 2021/2024 PHG Table A6.3. |
| Han X, Nabb DL, Russell MH, Kennedy GL, Rickard RW (2012). Renal elimination of perfluorocarboxylates (PFCAs). *Chem Res Toxicol* 25(1):35–46. | 21985250 | 10.1021/tx200363w | Review abstract: organic anion transport proteins drive sex-, species- and chain-length-dependent tubular reabsorption; a basolateral efflux pathway is proposed as the missing step | PARTIAL SUCCESS (abstract). Full text FAILED — ACS paywall, no PMC record. No new primary kinetic numbers (it is a review). |
| Kimura O, Fujii Y, Haraguchi K, Kato Y, Ohta C, Koga N, Endo T (2017). Uptake of perfluorooctanoic acid by Caco-2 cells: involvement of organic anion transporting polypeptides. *Toxicol Lett* 277:18–23. | 28552774 | 10.1016/j.toxlet.2017.05.012 | Km 8.3 µM; uptake clearance Vmax/Km = 55.0 µL/mg protein/min; non-saturable clearance 18.1 µL/mg protein/min; BSP Ki 23.1 µM; Na⁺-independent; not inhibited by probenecid or PAH | PARTIAL SUCCESS (abstract carries all numbers). Full text FAILED — Elsevier paywall, no PMC record. |
| Han X, Yang CH, Snajdr SI, Nabb DL, Mingoia RT (2008). Uptake of perfluorooctanoate in freshly isolated hepatocytes from male and female rats. *Toxicol Lett* 181(2):81–86. | 18662756 | 10.1016/j.toxlet.2008.07.002 | Male Km 88.0 ± 9.1, Vmax 5.61 ± 0.88; female Km 76.1 ± 12.0, Vmax 3.59 ± 0.29 nmol/min/10⁶ cells; clearances 64.8 vs 47.6 µL/min/10⁶ cells; BSP Ki 85.9 vs 29.3 µM | PARTIAL SUCCESS (abstract). Full text FAILED — Elsevier paywall. |
| Kudo N, Katakura M, Sato Y, Kawashima Y (2002). Sex hormone-regulated renal transport of perfluorooctanoic acid. *Chem Biol Interact* 139(3):301–316. | 11879818 | 10.1016/s0009-2797(02)00006-6 | Male t½ 70× female; castration raised CLR 14-fold; oatp1 and OAT-K mRNA male-predominant; OAT2 mRNA in males only 13% of females; ovariectomy raised OAT3 mRNA | PARTIAL SUCCESS (abstract). Full text FAILED — Elsevier paywall. |
| Louisse J, Pedroni L, van den Heuvel JJMW, et al. (2024). In vitro and in silico characterization of the transport of selected PFCAs and PFSAs by human OAT1, OAT2 and OAT3. *Toxicology* 509:153961. | 39343156 | 10.1016/j.tox.2024.153961 | OAT1 substrates: PFHpA, PFOA only. OAT2: none. OAT3 substrates: PFHpA, PFOA, PFNA, PFHxS. HEK293-transduced cells | PARTIAL SUCCESS (abstract via Firecrawl on doi.org, 2026-10-01); **Table 1 SUBSEQUENTLY OBTAINED** and held at `db/primary_2026/louisse2024_human_oat_km.csv`. Six genuinely new Km+Vmax pairs (OAT1 PFHpA 159/3.4, PFOA 185/3.5; OAT3 PFHpA 367/3.2, PFOA 90/1.5, PFNA 42/1.0, PFHxS 778/5.8); the four OAT4 rows it also tabulates are Louisse 2023's own, digit-identical, and are NOT re-counted. Folded into the record 2026-11-18 by `scripts/add_louisse2024_oat_km.py`. |
| Argoul CML, Toutain PL, Picard-Hagen N, Klukovits A, Koncsos G, Gayrard V, Lacroix MZ (2026). Comparative in vitro study of structure related uptake of ten perfluoroalkyl substances by human and rat renal transporters. *Arch Toxicol*. | 42446666 | 10.1007/s00204-026-04508-7 | Direct human-vs-rat comparison across 8 human and 4 rat renal transporters; OAT1 specifically mediated PFHxA uptake; human OAT3/OAT4 substrates = PFHxA, PFOA, PFNA, PFBS, PFHxS; rat Oat1 and Oatp1a1 transported the same short/medium-chain set; GenX and PFO2OA are substrates of hOAT3, rOat1 and rOatp1a1; long-chain PFDA, PFOS, PFDS were substrates of none | **FULL TEXT SUPPLIED BY THE REVIEWER 2026-11-19.** Ingested by `scripts/add_argoul2026_renal_transporters.py`: 66 rows, 6 transporters x 10 PFAS plus 6 that transported nothing. Uptake ratios only, NO Km/Vmax. Table 1 (human net tubular reabsorption) stored separately at `db/primary_2026/argoul2026_table1_human_reabsorption.csv`. Supplemental Tables A1-A9 NOT supplied. |
| Ebert A, Allendorf F, Berger U, Goss KU, Ulrich N (2020). Membrane/water partitioning and permeabilities of perfluoroalkyl acids and four of their alternatives … *Environ Sci Technol* 54(8):5051–5061. | 32212724 | 10.1021/acs.est.0c00175 | Measured anionic membrane permeabilities "high enough to explain observed cellular uptake by passive diffusion with no need to postulate the existence of active uptake processes"; COSMO-RS log RMSE 0.61 / 0.46 | **FULL TEXT SUPPLIED BY THE REVIEWER 2026-11-19.** Ingested by `scripts/add_ebert2020_membrane.py` into `db/membrane_partitioning.csv` (9 rows) — NOT the transporter tables, since Kmem/w is a different quantity. **The abstract overstates the body**: the no-active-uptake claim rests on control-cell uptake rates, i.e. the passive baseline the transporter papers already subtract, so it contradicts no Km here. Per-compound SI Tables S3.1.1 / S3.2.1 NOT supplied — a live request, and the figures must not be digitised in their place. |
| Liu YZ, Pan LH, Bai Y, Yang K, Dong PP, Fang ZZ (2020). Per- and polyfluoroalkyl substances exert strong inhibition towards human carboxylesterases. *Environ Pollut* 263(Pt A):114463. | 32283456 | 10.1016/j.envpol.2020.114463 | CES1 IC50: PFDoA 10.6, PFTA 13.4, PFOcDA 12.6 µM; CES2 IC50: 9.56, 17.2, 8.73 µM; noncompetitive Ki 27.7 / 26.9 / 11.9 (CES1) and 4.04 / 29.1 / 27.4 µM (CES2); IVIVE thresholds 0.40–2.91 µM plasma | PARTIAL SUCCESS (abstract carries all values). Full text FAILED — Elsevier paywall. **This is an enzyme, not a transporter** — included because it was on the requested list. |
| Kummu M, Sieppi E, Koponen J, Laatio L, Vähäkangas K, Kiviranta H, Rautio A, Myllynen P (2015). Organic anion transporter 4 (OAT4) modifies placental transfer of PFOS and PFOA in human placental ex vivo perfusion system. *Placenta* 36(10):1185–1191. | 26303760 | 10.1016/j.placenta.2015.07.119 | OAT4 protein negatively correlated with transfer index of PFOA (R²=0.92, p=0.043) and PFOS (R²=0.99, p=0.007) at 120 min; ABCG2 no correlation; n=4 for the correlation | PARTIAL SUCCESS (abstract). Full text FAILED — Elsevier paywall. |
| Buist SCN, Klaassen CD (2004). Rat and mouse differences in gender-predominant expression of organic anion transporter (Oat1-3; Slc22a6-8) mRNA levels. *Drug Metab Dispos* 32(6):620–625. | 15155553 | 10.1124/dmd.32.6.620 | Mouse renal Oat1 male-predominant; mouse renal Oat2 NO sex predominance (vs strongly female-predominant in rat); mouse hepatic Oat3 male-predominant but less so than rat | PARTIAL SUCCESS (abstract). Full text FAILED — ASPET paywall. **This is the central direct rat-vs-mouse Oat comparison.** |
| Buist SCN, Cherrington NJ, Choudhuri S, Hartley DP, Klaassen CD (2002). Gender-specific and developmental influences on the expression of rat organic anion transporters. *J Pharmacol Exp Ther* 301(1):145–151. | 11907168 | 10.1124/jpet.301.1.145 | Rat renal Oat1 male > female at days 40–45; renal Oat2 much higher in females; hepatic Oat3 male-predominant | PARTIAL SUCCESS (abstract). Full text FAILED — ASPET paywall. |
| Buist SCN, Cherrington NJ, Klaassen CD (2003). Endocrine regulation of rat organic anion transporters. *Drug Metab Dispos* 31(5):559–564. | 12695343 | 10.1124/dmd.31.5.559 | Androgens (not GH) cause the renal Oat1 sex difference; the female GH pattern causes the renal Oat2 sex difference | PARTIAL SUCCESS (abstract). Full text FAILED — ASPET paywall. |
| Cheng X, Maher J, Chen C, Klaassen CD (2005). Tissue distribution and ontogeny of mouse organic anion transporting polypeptides (Oatps). *Drug Metab Dispos* 33(7):1062–1073. | 15843488 | 10.1124/dmd.105.003640 | **Mouse kidney: Oatp1a1, Oatp3a1 and Oatp4c1 all higher in males than females**; Oatp1a1/1a6/3a1/4c1 are the highly expressed renal Oatps | PARTIAL SUCCESS (abstract). Full text FAILED — ASPET paywall. **This is the mouse counterpart to the rat Oatp1a1 sex-difference argument.** |
| Cheng X, Maher J, Lu H, Klaassen CD (2006). Endocrine regulation of gender-divergent mouse organic anion-transporting polypeptide (Oatp) expression. *Mol Pharmacol* 70(4):1291–1297. | 16807376 | 10.1124/mol.106.025122 | Androgens increased Oatp1a1 mRNA in mouse liver AND kidney; male-pattern GH increased hepatic but not renal Oatp1a1; "In kidney, gender-divergent Oatp expression is exclusively caused by stimulation by androgens" | PARTIAL SUCCESS (abstract). Full text FAILED — ASPET paywall. ⚠ **Internal contradiction flagged**: one sentence of this abstract calls renal Oatp1a1/3a1 "female-predominant", which conflicts both with Cheng 2005 (PMID 15843488) and with this paper's own androgen-stimulation conclusion. |

## Failures (no usable content obtained)

| Target | Why it was wanted | Attempt | Result |
|---|---|---|---|
| Katakura M, Kudo N, Tsuda T, Hibino Y, Mitsumoto A, Kawashima Y (2007). Rat organic anion transporter 3 and organic anion transporting polypeptide 1 mediate perfluorooctanoic acid transport. *J Health Sci* 53(1):77–83. | The *Xenopus* oocyte evidence for rat Oat3 and Oatp1a1, plus the Mrp2-deficient-rat and rat Npt2 negative results | PubMed search (`Katakura perfluorooctanoic rat kidney transport`) and Europe PMC author search `AUTH:"Katakura M" AND perfluorooctanoic` | **FAILED** — not indexed in PubMed or Europe PMC; *Journal of Health Science* (J-STAGE) has no DOI located. Findings recovered second-hand from OEHHA 2021/2024 PHG Table A6.3 and from the Weaver 2010 Discussion. |
| Dankers AC et al. (2013), described in the user's list as "Toxicol Sci — OAT4 and URAT1 and PFOS placental" | Requested paper | Targeted PubMed Boolean search | **FAILED** — no such paper found. The only Dankers 2013 *Toxicol Sci* paper indexed (PMID 24014645, DOI 10.1093/toxsci/kft198) is "Endocrine disruptors differentially target ATP-binding cassette transporters in the blood-testis barrier and affect Leydig cell testosterone secretion in vitro" — ABC transporters in the blood–testis barrier, **not** OAT4/URAT1 and not placental. The requested paper appears to be mis-cited in the source citation list. |
| "Niu Z et al. 2023/2024" | Requested papers | PubMed and Europe PMC | **RESOLVED AS MIS-CITATION** — the two relevant papers are Niu **S** et al. 2023 *EHP* (review, PMID 37418334) and Niu **S** et al. 2026 *Environ Health* (PMID 42775203); the 2024 *Chemosphere* PFAS/OAT screen is Ryu S et al. (PMID 38801906), not Niu. |
| Chou WC, Lin Z (2019/2021) PBPK with transporters | Requested | PubMed Boolean search for Chou/Lin PBPK PFOS transporter | **FAILED** — the search returned no matching record from this container. Not covered in the CSV. |
| Fabrega F et al. | Requested | not searched (de-prioritised; PBPK not transporter kinetics) | **NOT ATTEMPTED** |
| Han X, Snow TA, Kemper RA, Jepson GW (2003). *Chem Res Toxicol* — binding, clearance and transport of PFOA in rat and human plasma | Requested; plasma protein binding | Covered by the parallel protein-binding agent (`papers/` already contains several albumin-binding papers) | **NOT ATTEMPTED by this agent** (out of transporter scope; flagged for the protein-binding strand) |
| Per-congener Km and Vmax values of Niu 2026 for P-gp, BCRP and OAT1 | Needed for the efflux row of the database | Full text read (PMC13595387) | **PARTIAL FAILURE** — values exist only inside Figures 1f and 2b,d; only ranges and Vmax/Km ratios are stated in text. |
| PDF files generally | Deliverable requested `.pdf` + `.txt` pairs | `pmc.ncbi.nlm.nih.gov/articles/<PMCID>/pdf/` with browser UA; `link.springer.com/content/pdf/...`; `sciencedirect.com/.../pdfft` | **FAILED for all transporter papers** — PMC returns a ~1.8 KB anti-bot HTML interstitial; Springer and Elsevier return HTML landing pages. Full text was therefore saved as `.txt` only. The pre-existing regulatory PDFs (ATSDR, EPA, OEHHA) in `papers/` were downloaded by other agents and are present. |

---

## Reconciliation of this log against the files, 2026-11-18

Stage 0 of `report/TRANSPORTERS_ROADMAP.md`. Every claim above was re-walked
against `db/transporter_kinetics.csv` and `db/primary_2026/*.csv`. Four
findings, two of which correct this log or that roadmap rather than the data.

**1. Louisse 2024 was the only real sync gap, and it was smaller than stated.**
Table 1 holds fifteen rows, ten with a Km, but four of those ten are Louisse
2023's OAT4 values re-tabulated — digit-identical to the rows already in the
record from the 2023 primary (60/4.5, 47/4.5, 58/8.5, 92/7.3), and the
re-tabulation even loses precision, printing the standard error as 15 for all
four where the primary gives 33 for PFHxS. Six values were genuinely new, not
ten. Entering the other four would have turned one measurement into two.

**2. Han 2012 Table 6 needs no fold-in: it is a compilation, not a
measurement.** All thirteen of its values are already in the record, each
attributed to the paper that measured it — Weaver 2010, Nakagawa 2008 and Yang
2009 — which is the correct attribution. Table 6 also silently picks between
conflicting primaries where the record keeps both: for rat Oat1/PFOA it prints
Nakagawa's 51 where Weaver measured 43.2, for rat Oat3/PFOA Nakagawa's 80.2
against Weaver's 65.7, and for rat Oatp1a1/PFOA Yang 2009's 162.2 against
Weaver's 126.4. `db/primary_2026/han2012_table6_transporter_km.csv` also
mis-assigns one value: the 64.1 µM in its third OAT4 row is URAT1's, from Yang
2010, because Han's Table 6 formats OAT4 and URAT1 as a single row. The record
has it correctly under URAT1; the extract's note has been corrected.

**3. Yang 2010 was already synced.** All four of its rows, the three Km values
and the OATP1A2 negative result, were already in the record. The roadmap's
claim that it "needs the same check" was wrong.

**4. The record holds 27 sources, not 29.** The roadmap's §1 figure was an
error in that survey, not a change from this reconciliation: the pre-Stage-0
file has 27 distinct study strings and 28 DOIs, and still does.

Row count went 217 → 230. The two second-hand Vmax values recorded at Yang
2010 above are now gated out of any computed clearance by
`scripts/build_transporter_working_table.py`, so Yang 2010's own Vmax remains
a live request.

---

## Four priority full texts received, 2026-11-19

Requests 1, 2, 3 and 6 of `report/TRANSPORTERS_ROADMAP.md` §6, supplied by the
reviewer. Rows added to `db/transporter_kinetics.csv`: 230 -> 301.

| Paper | Ingested by | What it added | Still missing |
|---|---|---|---|
| **Argoul 2026** | `add_argoul2026_renal_transporters.py` | 66 rows: 60 per-compound substrate calls across 6 active transporters (28 positive) + 6 transporters that transported nothing. Uptake ratios only | Supplemental Tables A1–A9 (per-compound uptake values) |
| **Niu 2026** | `fix_niu2026_from_figures.py` | 9 Km+Vmax pairs read off Figures 2f and 3c,f, replacing 8 placeholder rows; plus 3 rows without constants | SI; the two Cl-ether congeners have no constants (low solubility) |
| **Nakagawa 2008** | `add_nakagawa2008_oat_vmax.py` | 4 Vmax values from Table 1, completing the four Km values already held | nothing — the paper is fully read |
| **Ebert 2020** | `add_ebert2020_membrane.py` | 9 rows to `db/membrane_partitioning.csv`; a different quantity | **SI Tables S3.1.1 and S3.2.1**, the per-compound log Kmem/w and log Pion |

**Two defects found inside the new papers**, both by checking each against its
own internal redundancy:

1. **Niu 2026 transposes P-gp and BCRP** in its Vmax/Km paragraph and abstract.
   The figures' own constants divide to give the attribution reversed, on all
   six efflux values; the abstract's Km ranges and the OAT1 arithmetic both
   confirm the figures are right. The abstract's headline claim that the
   highest Vmax/Km is PFOA by P-gp should read **BCRP**.
2. **Nakagawa 2008's Table 1** labels a Km/Vmax column with the units of
   Vmax/Km.

**One defect found in our own tooling**: the working-table builder was parsing
numbers out of placeholder prose, turning `"reported only in Fig 1f"` into a
Km of 1. Ten phantom Km values resulted, all from the Niu placeholder rows now
replaced. Guarded in `build_transporter_working_table.py`.

**Requests that remain open**, in order: Yang 2010 (its own Vmax, now the top
request — see roadmap §3.1b), Yang 2009, Ebert's SI, and **rat renal
transporter expression levels**, which roadmap §8 argues is now the single
most valuable missing measurement in this strand.
