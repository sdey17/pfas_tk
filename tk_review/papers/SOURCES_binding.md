# SOURCES_binding.md — PFAS–protein binding literature retrieval log

Scope: in vitro measurements of PFAS binding to serum albumin and other toxicokinetically
relevant proteins (globulins, L-FABP, TTR/TBG, thyroid hormone receptor, PPARα,
α2u-globulin, β-lactoglobulin, lipoproteins, phospholipid membranes).

Retrieval dates: all attempts 2026-10-01.
Retrieval route unless stated: Europe PMC REST API
(`/search` for metadata+abstract, `/{PMCID}/fullTextXML` for full text with tables).
Files written to `tk_review/papers/`.

**Important caveat on PDFs:** every attempt to download a PDF from
`https://www.ncbi.nlm.nih.gov/pmc/articles/{PMCID}/pdf/` and from the Europe PMC
`fulltextRepo` endpoint returned a non-PDF body (HTML interstitial / blocked) through
the agent proxy. **No PDFs were obtained in this work stream.** Where full text was
available I saved a `.txt` extraction from the JATS full-text XML instead, with the
tables flattened and appended at the end of each file under
`===== TABLES (flattened) =====`. Paywalled papers are recorded abstract-only, and the
`protein_binding.csv` rows derived from them are marked in `source_table` as
"abstract" or "secondary (Starnes & Belcher 2026 Table 1)".

---

## A. Full text retrieved and saved (.txt with flattened tables)

| File (.txt) | Full citation | PMID | PMC | DOI | URL | Result |
|---|---|---|---|---|---|---|
| `Starnes 2024 cross-species serum albumin DSF.txt` | Starnes HM, Jackson TW, Rock KD, Belcher SM. Quantitative cross-species comparison of serum albumin binding of per- and polyfluoroalkyl substances from five structural classes. Toxicol Sci. 2024;199:132–149. | 38518100 | PMC11057469 | 10.1093/toxsci/kfae028 | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11057469/fullTextXML | SUCCESS (Tables 1–5 extracted) |
| `Starnes 2026 albumin PFAS binding methods review.txt` | Starnes HM, Belcher SM. Toward Comprehensive In Vitro Evaluation of Serum Albumin Binding of Per- and Polyfluoroalkyl Substances. J Xenobiot. 2026;16:54. | 41874125 | PMC13010723 | 10.3390/jox16020054 | https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13010723/fullTextXML | SUCCESS (Table 1 = master cross-study Ka compilation, 33 rows) |
| `Starnes 2025 65 PFAS HSA DSF machine learning.txt` | Starnes HM, Green AJ, Reif DM, Belcher SM. An in vitro and machine learning framework for quantifying serum albumin binding of per- and polyfluoroalkyl substances. Toxicol Sci. 2025;203:67–78. | 39298512 | PMC11664106 | 10.1093/toxsci/kfae124 | Europe PMC fullTextXML | SUCCESS (Table 1 = binary binding calls; per-congener Kd in SI, not retrieved) |
| `Alesio 2024 BSA PFAS differential scanning fluorimetry.txt` | Alesio J, Bothun GD. Differential scanning fluorimetry to assess PFAS binding to bovine serum albumin protein. Sci Rep. 2024;14:6501. | 38499613 | PMC10948889 | 10.1038/s41598-024-57140-9 | Europe PMC fullTextXML | SUCCESS (Table 2 = Ka at 61 °C) |
| `Li 2021 AFFF PFAS HSA equilibrium dialysis docking.txt` | Li W, Hu Y, Bischel HN. In-Vitro and In-Silico Assessment of PFAS in AFFF Binding to Human Serum Albumin. Toxics. 2021;9:63. | 33803062 | PMC8002870 | 10.3390/toxics9030063 | Europe PMC fullTextXML | PARTIAL — body text retrieved; the 26-congener Ka table is in Supporting Information and was NOT retrieved. Only PFPrS log Ka = 4.1 ± 0.2 is quoted in the body. |
| `Peng 2024 six PFAS HSA multispectroscopy site competition.txt` | Peng M, Xu Y, Wu Y, Cai X, Zhang W, Zheng L, Du E, Fu J. Binding Affinity and Mechanism of Six PFAS with Human Serum Albumin. Toxics. 2024;12:43. | 38250999 | PMC10819430 | 10.3390/toxics12010043 | Europe PMC fullTextXML | SUCCESS (Tables 2–4: Kb, n, thermodynamics, warfarin/ibuprofen/lidocaine competition, MM-GBSA) |
| `Wu 2024 legacy and novel PFAS HSA in vitro in silico.txt` | Wu Y, Bao J, Liu Y, Wang X, Lu X, Wang K. In Vitro and In Silico Analysis of the Bindings between Legacy and Novel PFAS and Human Serum Albumin. Toxics. 2024;12:46. | 38251003 | PMC10818824 | 10.3390/toxics12010046 | Europe PMC fullTextXML | SUCCESS (Tables 1–7: Ksv, Kq, KA, n, docking) |
| `Chen 2020 perfluoroalkyl carboxylic acids HSA spectroscopy.txt` | Chen H, Wang Q, Cai Y, Yuan R, Wang F, Zhou B. Investigation of the Interaction Mechanism of Perfluoroalkyl Carboxylic Acids with Human Serum Albumin by Spectroscopic Methods. Int J Environ Res Public Health. 2020;17:1319. | 32085632 | PMC7068604 | 10.3390/ijerph17041319 | Europe PMC fullTextXML | SUCCESS (body + tables) |
| `Moro 2022 branched short-chain PFAS HSA ITC.txt` | Moro G, Liberi S, Vascon F, et al. Investigation of the Interaction between Human Serum Albumin and Branched Short-Chain Perfluoroalkyl Compounds. Chem Res Toxicol. 2022;35:2049–2058. | 36148994 | PMC9682524 | 10.1021/acs.chemrestox.2c00211 | Europe PMC fullTextXML | SUCCESS (Table 1 = ITC KD, n, ΔH, ΔG at 298/310 K) |

## B. Already on disk (read, not re-retrieved — per assignment instruction)

| File | Citation | DOI | Result |
|---|---|---|---|
| `literature/fulltext/Fischer et al. 2024 Environ Sci Technol.txt` | Fischer FC, Ludtke S, Thackray C, Pickard HM, Haque F, Dassuncao C, Endo S, Schaider L, Sunderland EM. Binding of PFAS to Serum Proteins: Implications for Toxicokinetics in Humans. Environ Sci Technol. 2024;58:1055–1063. | 10.1021/acs.est.3c07415 | SUCCESS — **Table S2 fully recovered** (log D_C18/w, D_HSA/w, D_BSA/w, D_glob/w, D_serum/w for 32 PFAS). Column order verified by block counting (11 carboxylates / 7 sulfonates / 4 sulfonamides per column). |
| `literature/fulltext/Fischer 2025.txt` | Fischer FC, Thackray C, Ferguson N, et al. Understanding Mechanisms of PFAS Absorption, Distribution, and Elimination Using a PBTK Model. Environ Sci Technol. 2025;59:13240–13250. | 10.1021/acs.est.5c05473 | SUCCESS — read for TK linkage (sensitivity of elimination to albumin binding vs phospholipid binding vs transporters). |
| `tk_review/papers/Ryu 2024 unbound fractions PFAS rodent tissues(_raw_xml).txt` | Ryu S, Burchett W, Zhang S, et al. Unbound fractions of PFAS in human and rodent tissues: rat liver a suitable proxy for evaluating emerging PFAS? Environ Sci Technol. 2024;58:14641–14650. | 10.1021/acs.est.4c04050 | SUCCESS (file placed by a parallel work stream; read here for tissue f_unbound). PMID 39161261, PMC11825104. |

## C. Abstract-only (paywalled; metadata + abstract via Europe PMC `/search?resultType=core`)

All of the following returned a usable abstract containing the quantitative result quoted
in `protein_binding.csv`; full text and tables were NOT obtained.

| Citation | PMID | DOI | Result |
|---|---|---|---|
| Bischel HN, MacManus-Spencer LA, Luthy RG. Noncovalent interactions of long-chain perfluoroalkyl acids with serum albumin. Environ Sci Technol. 2010;44:5263–5269. | 20540534 | 10.1021/es101334s | ABSTRACT ONLY |
| Bischel HN, MacManus-Spencer LA, Zhang C, Luthy RG. Strong associations of short-chain perfluoroalkyl acids with serum albumin and investigation of binding mechanisms. Environ Toxicol Chem. 2011;30:2423–2430. | 21842491 | 10.1002/etc.647 | **PDF OBTAINED 2026-11-03** (`literature/fulltext/Bischel 2011 Environ Toxicol Chem short-chain PFAA serum albumin.pdf`). Table 1 read at source: 11 log K_PW with SE, reservoir ranges and fraction bound; every value already in the database confirmed. NOTE the pdftotext extract renders uM as mM throughout (it drops the mu, as it drops the plus-minus) -- read the PDF, not the .txt |
| MacManus-Spencer LA, Tse ML, Hebert PC, Bischel HN, Luthy RG. Binding of perfluorocarboxylates to serum albumin: a comparison of analytical methods. Anal Chem. 2010;82:974–981. | 20039637 | 10.1021/ac902238u | ABSTRACT ONLY |
| Hebert PC, MacManus-Spencer LA. Development of a fluorescence model for the binding of medium- to long-chain perfluoroalkyl acids to human serum albumin. Anal Chem. 2010;82. [page range not verified] | 20590160 | 10.1021/ac100721e | ABSTRACT ONLY |
| Allendorf F, Berger U, Goss KU, Ulrich N. Partition coefficients of four perfluoroalkyl acid alternatives between bovine serum albumin (BSA) and water in comparison to ten classical perfluoroalkyl acids. Environ Sci Process Impacts. 2019;21:1852–1863. | 31475719 | 10.1039/C9EM00290A | **PDF + SI OBTAINED 2026-11-03**. Table 1 and SI confirm the 24 ingested rows. Three findings: the ABSTRACT's range 2.8-4.8 is a typo (body and Table 1 give 2.5-4.8, PFBA 2.52); the duplicated Ka cells for HFPO-DA and DONA are the PAPER's error, not a transcription slip (11 of 14 reproduce from their own log K within 5%, those two sit 3.3x off); and SI-2.4 says BOTH batches are heat-shock fraction, A3803 differing only by charcoal defatting. SI Table SI-2.3 adds the measured analyte:BSA ratios, 4e-4 to 3.8e-2 |
| Ebert A, Allendorf F, Berger U, Goss KU, Ulrich N. Membrane/Water Partitioning and Permeabilities of Perfluoroalkyl Acids and Four of their Alternatives. Environ Sci Technol. 2020;54. [page range not verified] | 32212724 | 10.1021/acs.est.0c00175 | ABSTRACT ONLY |
| Forsthuber M, Kaiser AM, Granitzer S, et al. Albumin is the major carrier protein for PFOS, PFOA, PFHxS, PFNA and PFDA in human plasma. Environ Int. 2020;art. 105324. [volume not verified] | 32109724 | 10.1016/j.envint.2019.105324 | ABSTRACT ONLY |
| Han X, Snow TA, Kemper RA, Jepson GW. Binding of perfluorooctanoic acid to rat and human plasma proteins. Chem Res Toxicol. 2003;16:775–781. | 12807361 | 10.1021/tx034005w | ABSTRACT ONLY |
| Han X, Hinderliter PM, Snow TA, Jepson GW. Binding of perfluorooctanoic acid to rat liver-form and kidney-form alpha2u-globulins. Drug Chem Toxicol. 2004;27. [page range not verified] | 15573471 | 10.1081/dct-200039725 | ABSTRACT ONLY |
| Beesoon S, Martin JW. Isomer-specific binding affinity of PFOS and PFOA to serum proteins. Environ Sci Technol. 2015;49. [page range not verified] | 25826685 | 10.1021/es505399w | ABSTRACT ONLY |
| Luo Z, Shi X, Hu Q, Zhao B, Huang M. Structural evidence of perfluorooctane sulfonate transport by human serum albumin. Chem Res Toxicol. 2012. [volume/pages not verified] | 22482699 | 10.1021/tx300112p | ABSTRACT ONLY (PDB 4E99 site assignment taken from Li et al. 2021 body text) |
| Maso L, Trande M, Liberi S, et al. Unveiling the binding mode of perfluorooctanoic acid to human serum albumin. Protein Sci. 2021;30. [page range not verified] | 33550662 | 10.1002/pro.4036 | ABSTRACT ONLY |
| Zhang X, Chen L, Fei XC, Ma YS, Gao HW. Binding of PFOS to serum albumin and DNA: insight into the molecular toxicity of perfluorochemicals. BMC Mol Biol. 2009;10:16. | 19239717 | 10.1186/1471-2199-10-16 | ABSTRACT ONLY (PMC2656506 exists; full text not fetched) |
| Chen H, He P, Rao H, Wang F, Liu H, Yao J. Systematic investigation of the toxic mechanism of PFOA and PFOS on bovine serum albumin by spectroscopic and molecular modeling. Chemosphere. 2015;129:217–224. | 25497588 | 10.1016/j.chemosphere.2014.11.040 | ABSTRACT ONLY |
| Chen YM, Guo LH. Fluorescence study on site-specific binding of perfluoroalkyl acids to human serum albumin. Arch Toxicol. 2009;83:255–261. | 18854981 | 10.1007/s00204-008-0359-x | ABSTRACT ONLY (abstract contains all four Ka values) |
| Qin P, Liu R, Pan X, Fang X, Mou Y. Impact of carbon chain length on binding of perfluoroalkyl acids to bovine serum albumin determined by spectroscopic methods. J Agric Food Chem. 2010;58:5561–5567. | 20397730 | 10.1021/jf100412q | ABSTRACT ONLY |
| Chi Q, Li Z, Huang J, Ma J, Wang X. Interactions of PFOA and PFOS with serum albumins by native mass spectrometry, fluorescence and molecular docking. Chemosphere. 2018;198:442–449. | 29425944 | 10.1016/j.chemosphere.2018.01.152 | ABSTRACT ONLY |
| Fedorenko M, Alesio J, Fedorenko A, Slitt A, Bothun GD. Dominant entropic binding of PFASs to albumin protein revealed by 19F NMR. Chemosphere. 2021;263:128083. | 33297081 | 10.1016/j.chemosphere.2020.128083 | ABSTRACT ONLY (PMC8479757; XML not fetched) |
| Alesio JL, Slitt A, Bothun GD. Critical new insights into the binding of poly- and perfluoroalkyl substances (PFAS) to albumin protein. Chemosphere. 2022;287:131979. | 34450368 | 10.1016/j.chemosphere.2021.131979 | ABSTRACT ONLY — `fullTextXML` returned an empty (150-byte) body despite PMC8612954 being listed |
| Crisalli AM, Cai A, Cho BP. Probing the interactions of perfluorocarboxylic acids of various chain lengths with human serum albumin: calorimetric and spectroscopic investigations. Chem Res Toxicol. 2023;36:703–713. | 37001030 | 10.1021/acs.chemrestox.3c00011 | ABSTRACT ONLY — `fullTextXML` empty (author manuscript PMC11091765 not served) |
| Jackson TW, Scheibly CM, Polera ME, Belcher SM. Rapid characterization of human serum albumin binding for PFAS using differential scanning fluorimetry. Environ Sci Technol. 2021;55:12291–12301. | 34495656 | 10.1021/acs.est.1c01200 | ABSTRACT ONLY — `fullTextXML` empty (PMC8651256 not served) |
| Weiss JM, Andersson PL, Lamoree MH, Leonards PEG, van Leeuwen SPJ, Hamers T. Competitive binding of poly- and perfluorinated compounds to the thyroid hormone transport protein transthyretin. Toxicol Sci. 2009;109(2):206– . [end page not verified] | 19293372 | 10.1093/toxsci/kfp055 | ABSTRACT ONLY. Also attempted `WebFetch` of https://academic.oup.com/toxsci/article/109/2/206/1667072 — returned abstract-level content only; **per-compound IC50 table NOT obtained** |
| Ren XM, Zhang YF, Guo LH, Qin ZF, Lv QY, Zhang LY. Structure–activity relations in binding of perfluoroalkyl compounds to human thyroid hormone T3 receptor. Arch Toxicol. 2015;89. [page range not verified] | 24819616 | 10.1007/s00204-014-1258-y | ABSTRACT ONLY |
| Ren XM, Qin WP, Cao LY, et al. Binding interactions of perfluoroalkyl substances with thyroid hormone transport proteins and potential toxicological implications. Toxicology. 2016. [volume/pages not verified] | 27528273 | 10.1016/j.tox.2016.08.011 | ABSTRACT ONLY |
| Luebker DJ, Hansen KJ, Bass NM, Butenhoff JL, Seacat AM. Interactions of fluorochemicals with rat liver fatty acid-binding protein. Toxicology. 2002;176. [page range not verified] | 12093614 | 10.1016/s0300-483x(02)00081-1 | ABSTRACT ONLY |
| Zhang L, Ren XM, Guo LH. Structure-based investigation on the interaction of perfluorinated compounds with human liver fatty acid binding protein. Environ Sci Technol. 2013;47. [page range not verified] | 24006842 | 10.1021/es4026722 | **PDF OBTAINED 2026-11-03**. Table 1 read at source and ingested as primary by `scripts/add_zhang2013_hlfabp.py`: 19 rows, IC50 and Kd for 17 PFCs plus oleic and palmitic acid, probe constants and the 400 nM / 50 nM assay concentrations. Supersedes the transcription from Fan 2025 Table 1, which it reproduces exactly on all 10 shared values |
| Sheng N, Li J, Liu H, Zhang A, Dai J. Interaction of perfluoroalkyl acids with human liver fatty acid-binding protein. Arch Toxicol. 2016;90. [page range not verified] | 25370009 | 10.1007/s00204-014-1391-7 | ABSTRACT ONLY |
| Sheng N, Cui R, Wang J, Guo Y, Wang J, Dai J. Cytotoxicity of novel fluorinated alternatives to long-chain PFAS to human liver cell line and their binding capacity to human liver fatty acid binding protein. Arch Toxicol. 2018;92. [page range not verified] | 28864880 | 10.1007/s00204-017-2055-1 | ABSTRACT ONLY |
| Sheng N, Wang J, Guo Y, Wang J, Dai J. Interactions of PFOS and 6:2 chlorinated polyfluorinated ether sulfonate with human serum albumin: a comparative study. Chem Res Toxicol. 2020;33:1478–1486. | 32423201 | 10.1021/acs.chemrestox.0c00075 | ABSTRACT ONLY |
| Ishibashi H, Hirano M, Kim EY, Iwata H. In vitro and in silico evaluations of binding affinities of PFAS to Baikal seal and human PPARα. Environ Sci Technol. 2019;53. [page range not verified] | 30649875 | 10.1021/acs.est.8b07273 | ABSTRACT ONLY |
| Jia Y, Zhu Y, Xu D, et al. Insights into the competitive mechanisms of PFAS partition in liver and blood. Environ Sci Technol. 2022;56:6192–6200. | 35436088 | 10.1021/acs.est.1c08493 | ABSTRACT ONLY. `WebFetch` of https://pubs.acs.org/doi/10.1021/acs.est.1c08493 → **HTTP 403** (ACS blocks). **PDF OBTAINED 2026-11-03**, which supersedes that 403. Table 1 read at source: 11 hL-FABP Kd with dS, dH and dG, confirming every ingested value, and an HSA column that is Gao 2019 quoted by reference -- matching the Gao rows this database already held from Jia. STILL MISSING: Text S7, which holds the MST protein and ligand concentrations |
| Weiss-Errico MJ, Miksovska J, O'Shea KE. β-Cyclodextrin reverses binding of perfluorooctanoic acid to human serum albumin. Chem Res Toxicol. 2018;31. [page range not verified] | 29589912 | 10.1021/acs.chemrestox.8b00002 | ABSTRACT ONLY (qualitative site counts only; no Ka in abstract) |
| Wang Y, Zhang H, Kang Y, Cao J. Effects of perfluorooctane sulfonate on the conformation and activity of bovine serum albumin. J Photochem Photobiol B. 2016;159. [page range not verified] | 27031195 | 10.1016/j.jphotobiol.2016.03.024 | ABSTRACT ONLY (no Ka in abstract) |
| Sanchez Garcia D, Sjödin M, Hellstrandh M, et al. Cellular accumulation and lipid binding of perfluorinated alkylated substances (PFASs) — a comparison with lysosomotropic drugs. Chem Biol Interact. 2018;281. [page range not verified] | 29248446 | 10.1016/j.cbi.2017.12.021 | ABSTRACT ONLY |
| Han J, Fu J, Sun J, et al. Quantitative chemical proteomics reveals interspecies variations on binding schemes of L-FABP with perfluorooctanesulfonate. Environ Sci Technol. 2021;55. [page range not verified] | 34133149 | 10.1021/acs.est.1c00509 | ABSTRACT ONLY |
| Khazaee M, Ng CA. Evaluating parameter availability for PBPK modeling of perfluorooctanoic acid in zebrafish. Environ Sci Process Impacts. 2018;20. [page range not verified] | 29265128 | 10.1039/c7em00474e | ABSTRACT ONLY |
| Ng CA, Hungerbühler K. Bioconcentration of perfluorinated alkyl acids: how important is specific binding? Environ Sci Technol. 2013;47. [page range not verified] | 23734664 | 10.1021/es400981a | ABSTRACT ONLY |
| Ng CA, Hungerbühler K. Bioaccumulation of perfluorinated alkyl acids: observations and models. Environ Sci Technol. 2014;48. [page range not verified] | 24762048 | 10.1021/es404008g | ABSTRACT ONLY |
| Wilson DL, Chakraborty S, Sweety UH, et al. A molecular and cellular understanding of PFDA-exposure-associated outcomes on biological assemblies [β-lactoglobulin]. Environ Res. 2026;(in press). | 41748002 | 10.1016/j.envres.2026.124118 | ABSTRACT ONLY |

## D. Values taken at second hand from a compilation table (traced but not primary-verified)

The following appear in `protein_binding.csv` with
`source_table = "secondary (Starnes & Belcher 2026 Table 1)"`. They are **ranges across
congeners**, not per-congener values, and I did **not** reach the primary paper.

- Klevens HB, Ellenbogen E. Protein fluoroacid interaction: bovine serum albumin
  perfluoro-octanoic acid. Discuss Faraday Soc. 1954;18:277–288.
  doi:10.1039/df9541800277 — BSA/PFOA, equilibrium dialysis, Ka = 3.22 × 10² M⁻¹.
  **FAILURE: not indexed in PubMed/Europe PMC; not retrieved.**
- Wu LL, Gao HW, Gao NY, Chen FF, Chen L. Interaction of perfluorooctanoic acid with
  human serum albumin. BMC Struct Biol. 2009;9:31. doi:10.1186/1472-6807-9-31 —
  HSA/PFOA, equilibrium dialysis, Ka = 3.12 × 10⁴ M⁻¹. **Not separately retrieved.**
- Gao K, Zhuang T, Liu X, et al. Prenatal exposure to PFAS and association between
  placental transfer efficiencies and dissociation constant of serum protein–PFAS
  complexes. Environ Sci Technol. 2019. [volume/pages/DOI NOT verified — cited only as ref 55 of Starnes & Belcher 2026] —
  HSA, 14 PFAS, equilibrium dialysis, Ka = 0.24–2.63 × 10⁴ M⁻¹. **Not separately retrieved.**
- Peng SY, Yang YD, Tian R, Lu N. Critical new insights into the interactions of
  HFPO-DA (GenX) with albumin at molecular and cellular levels. J Environ Sci.
  2025;149:88–98. doi:10.1016/j.jes.2024.02.011 — HSA/BSA, 7 PFAS, fluorescence
  quenching, Ka = 0.14–2.92 × 10⁵ M⁻¹. **Not separately retrieved.**
- Kerstner-Wood C, Coward L, Gorman G. Protein binding of perfluorohexane sulfonate,
  perfluorooctane sulfonate and perfluorooctanoate to plasma (human, rat, and monkey)
  and various human-derived plasma protein fractions. Southern Research Institute
  Study ID 9921.7, 2003. **FAILURE: unpublished contract report, no DOI/PMID;
  not retrieved. This is the most-cited source for plasma f_unbound of PFHxS/PFOS/PFOA
  and remains a gap.**
- Cheng W, Ng CA. Predicting relative protein affinity of novel PFAS by an efficient
  molecular dynamics approach. Environ Sci Technol. 2018. [volume/pages/DOI NOT verified] — **FAILURE: the DOI I tried (10.1021/acs.est.7b03806)
  returned NOTFOUND; the correct record was not located. Not retrieved.**

## E. Search queries used (Europe PMC / PubMed)

1. PubMed: `(perfluoroalkyl OR perfluorooctanoic OR perfluorooctane sulfonate OR PFAS) AND (serum albumin binding OR protein binding)` — 487 hits.
2. Europe PMC: `PFAS serum albumin binding constant fluorescence quenching equilibrium dialysis` — 11 hits (found the Starnes papers, Alesio, Crisalli, Han 2021).
3. Europe PMC: `Starnes critical review per- and polyfluoroalkyl substances protein binding toxicokinetics` — 3 hits.
4. 13 targeted `DOI:"..."` lookups (section C/D DOIs).
5. 19 targeted title searches for the non-albumin proteins (TTR, TBG, TR, L-FABP, α2u-globulin, PPARα, β-lactoglobulin, haemoglobin, phospholipid membranes, lipoproteins).
6. 13 further title searches for Beesoon/Sheng/Zhang/Allendorf/Ng/Cheng/Khazaee/Maso/
   Kerstner-Wood/Hebert/β-lactoglobulin/haemoglobin/Jia.

## F. Explicit failures / gaps

- **No PDFs retrieved at all** (proxy blocks PMC and publisher PDF endpoints).
- **Haemoglobin**: no in vitro PFAS–haemoglobin binding constant found in any search.
  Searched "hemoglobin perfluorooctane sulfonate binding" — no primary quantitative hit.
- **Thyroxine-binding globulin (TBG)**: only relative potencies (Ren et al. 2016); no
  absolute Ka/Kd found.
- **Lipoproteins**: only the negative result of Forsthuber et al. 2020 (little or no
  affinity); no quantitative lipoprotein partition coefficients found.
- **Organic anion transporters as binding partners**: OAT work found in searches is
  transport kinetics (Km/IC50 for uptake), not equilibrium binding constants; deferred to
  the transporter work stream (files already in `papers/`: `Ryu 2024 14 PFAS OAT
  interactions.txt`, `Louisse 2023 OAT4 PFAS substrates.txt`, `Niu 2026 three renal
  transporters PFAS.txt`).
- **Surface plasmon resonance**: no SPR-derived PFAS–albumin binding constant was found
  in any of the searches run; Starnes & Belcher 2026 list SPR only as a method they did
  not find applied at scale. Recorded as a gap, not as zero evidence.
- Per-congener tables in Supporting Information were not reachable for
  Li et al. 2021, Starnes et al. 2025, Allendorf et al. 2019, Zhang et al. 2013,
  Weiss et al. 2009, and Jia et al. 2022.
  **Superseded 2026-11-03 for Allendorf (PDF + SI) and Zhang (PDF); Jia's
  main text is in hand but its Text S7 is not.** Li 2021, Starnes 2025 and
  Weiss 2009 remain open.
- **Why this sweep missed Fan et al. 2025** (Environ Int 201:109591), which
  measures rat albumin AND rat L-FABP in one laboratory and reports the
  stoichiometry: every query in section E was anchored on "serum albumin"
  or on a named protein or author, and Fan's title says only "liver and
  serum proteins". Nothing already held cites it. A title-anchored search
  cannot find a paper that names neither protein.
- **Khazaee et al. 2021 is five compounds on L-FABP, not seven.** PFHxA,
  PFOA, PFBS, PFHxS and PFOS are the L-FABP set; PFHpA and PFNA were run
  against I-FABP, a different protein. The database has always had this
  right; an analysis query that filtered on "fabp" swept both proteins
  into one comparison and produced a spurious chain-length inversion.
  Khazaee's L-FABP values are the lowest in the literature on all three
  compounds it reports (24x, 45x and 6.6x below the next), but there is
  no inversion and no basis for calling them unusable.
