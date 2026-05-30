# Reference: Human Liver Cell Types & Fibrosis-Associated Genes / Cell States

> Curated, literature-grounded reference used to (a) build `config/markers_liver.yaml`,
> (b) QA automated annotation, and (c) interpret fibrosis-associated genes/states.
> Marker sets are deliberately conservative (high-confidence canonical markers) and are
> meant to be *combined* (a cluster should express several markers of a type, not just one).
> Primary sources: Ramachandran et al. 2019 (*Nature*, GSE136103); MacParland et al. 2018;
> Aizarani et al. 2019; Andrews et al. 2022 (liver cell atlas); Dobie et al. 2019 (HSC
> zonation); Payen et al. 2021; the *Liver Cell Atlas* (livercellatlas.org). Gene symbols are
> human (HGNC).

---

## 1. Major liver cell compartments and canonical markers

| Compartment | Subtype / notes | Canonical markers (human) |
|---|---|---|
| **Hepatocytes** | parenchymal; zonated (periportal↔pericentral) | `ALB`, `APOA1`, `APOB`, `TTR`, `TF`, `SERPINA1`, `HP`, `CYP2E1` (pericentral), `CYP3A4`, `ASGR1`, `HNF4A` |
| **Cholangiocytes** | biliary epithelium; ductular reaction in fibrosis | `KRT19`, `KRT7`, `EPCAM`, `SOX9`, `CFTR`, `ANXA4`, `SPP1`, `FXYD2` |
| **Liver sinusoidal endothelial cells (LSEC)** | fenestrated sinusoidal endothelium | `CLEC4G`, `CLEC4M`, `STAB1`, `STAB2`, `FCGR2B`, `OIT3`, `CLEC1B`, `LYVE1` |
| **Vascular endothelial cells** | portal/central vein + arterial | `PECAM1` (CD31), `VWF`, `CD34`, `CLDN5`, `PLVAP`, `RAMP2`, `AQP1` |
| **Hepatic stellate cells (HSC) / mesenchyme** | pericytes of the sinusoid; quiescent vs activated | `PDGFRB`, `DCN`, `RGS5`, `DES`, `LRAT`, `RBP1`, `COLEC11`, `HGF`, `RELN`, `ANGPTL6` |
| **Portal fibroblasts / VSMC** | portal-tract mesenchyme; distinct from HSC | `PDGFRA`, `THY1` (CD90), `ELN`, `FBLN1`, `MSLN`, `SLIT2`, `GPC3`, `ACTA2` (VSMC) |
| **Kupffer cells / resident macrophages** | self-renewing tissue-resident | `CD68`, `CD163`, `MARCO`, `VSIG4`, `CD5L`, `TIMD4`, `C1QA`, `SLC40A1`, `VCAM1` |
| **Monocytes / monocyte-derived macrophages** | recruited; expand in injury | `CD14`, `LYZ`, `S100A8`, `S100A9`, `FCN1`, `VCAN`, `CCR2`, `FCGR3A` (CD16) |
| **Conventional / plasmacytoid DC** | antigen-presenting | cDC1 `CLEC9A`, `XCR1`; cDC2 `CD1C`, `FCER1A`; pDC `LILRA4`, `GZMB`, `IL3RA` |
| **T / NK cells** | cytotoxic + helper + tissue-resident NK | `CD3D`, `CD3E`, `CD8A`, `CD4`, `IL7R`, `NKG7`, `GNLY`, `KLRD1`, `CD7` |
| **B / plasma cells** | humoral | `CD79A`, `CD79B`, `MS4A1` (CD20), `IGHM`; plasma `MZB1`, `IGHG1`, `XBP1`, `JCHAIN` |
| **Mast cells** | granulocyte | `TPSAB1`, `TPSB2`, `CPA3`, `KIT`, `MS4A2` |
| **Cycling cells** | proliferating fraction (any lineage) | `MKI67`, `TOP2A`, `STMN1`, `PCNA` |

---

## 2. Disease-relevant fibrotic cell states (the validation targets)

These are the **fibrosis-associated states**, not just cell types — they emerge or expand in
the fibrotic/cirrhotic niche and are the headline findings the pipeline must recover and validate.

### 2a. Scar-associated macrophages (SAMac / "SAMΦ") — *macrophage compartment*
Pro-fibrogenic, monocyte-derived; topographically restricted to the fibrotic scar (Ramachandran
2019). Broadly conserved across organs as a `SPP1+`/`TREM2+` fibrogenic macrophage program.
- **Core markers:** `TREM2`, `CD9`, `SPP1`, `GPNMB`, `FABP5`, `CD63`, `LGALS3`
- **Supporting:** `CCR2`, `MNDA`, `FCN1` (monocyte origin), `LPL`, `ACP5`
- **Function / why it matters:** secretes `SPP1`, `PDGFB`, `TNFSF12` (TWEAK) → activates stellate
  cells; a prime source of secreted/surface biomarker candidates.

### 2b. Scar-associated / "lymphatic-like" endothelial cells (SAEndo) — *endothelial compartment*
Expand in the cirrhotic niche; promote leukocyte transmigration (Ramachandran 2019).
- **Core markers:** `ACKR1` (DARC), `PLVAP`, `VWA1`, `CD34`
- **Contrast:** healthy LSEC identity markers (`CLEC4G`, `STAB2`, `OIT3`) are **down** in SAEndo →
  loss of sinusoidal identity ("capillarization").
- **Function:** `VCAM1`/`ICAM1`-high; immune-recruiting; angiocrine `JAG1`/`DLL4` → NOTCH on HSC.

### 2c. Activated HSC → myofibroblast — *stellate / mesenchymal compartment*
The central collagen-producing, fibrogenic effector. A **continuum** from quiescent HSC to
activated myofibroblast (the screening-question Q4 cluster: COL1A1/COL3A1/ACTA2/TAGLN/PDGFRB/LUM/DCN).
- **Quiescent HSC:** `LRAT`, `RBP1`, `RGS5`, `DES`, `GFAP`, `ANGPTL6`, `HGF`, retinoid storage
- **Activated HSC / myofibroblast:** `ACTA2`, `TAGLN`, `COL1A1`, `COL1A2`, `COL3A1`, `TIMP1`,
  `TNC`, `LOX`, `LOXL2`, `MMP2`, `SPARC`, `POSTN`, `CTGF`(`CCN2`), `PDGFRB`(↑)
- **Pan-fibroblast / matrix:** `DCN`, `LUM`, `COL1A1`, `PDGFRB`
- **Translational HSC-derived biomarker:** `SMOC2` (tracks NAFLD/MASLD severity; secreted → plasma).

### 2d. Disambiguating the Q4 stromal cluster (HSC vs portal fibroblast vs myofibroblast vs mixed)
| State | Distinguishing markers (use *combinations*) |
|---|---|
| Activated **HSC** | `RGS5`, `DES`, `LRAT`, `RBP1`, `NGFR`, `PDGFRB`, retinoid/quiescence remnants |
| **Portal fibroblast** | `THY1`(CD90), `CD34`, `ELN`, `FBLN1`, `MSLN`, `SLIT2`, `GPC3`, `PDGFRA` |
| **Myofibroblast** (terminal) | `ACTA2`-high, `TAGLN`, `COL1A1`-high, `TIMP1`, `TNC`, `POSTN`, low quiescence markers |
| **Vascular smooth muscle** | `ACTA2`, `MYH11`, `NOTCH3`, `PLN` (perivascular, not scar) |
> Resolve with: sub-clustering of the mesenchymal compartment, the marker panels above,
> trajectory/diffusion ordering (quiescent→activated), spatial/zonation priors (HSC sinusoidal vs
> portal-fibroblast periportal), and cross-referencing the Ramachandran mesenchymal subsets.

---

## 3. Fibrosis pathways / mechanisms (for pathway & cell-cell-communication analysis)

| Pathway | Role in liver fibrosis | Representative genes / nodes |
|---|---|---|
| **TGF-β / SMAD** | master fibrogenic driver; HSC→myofibroblast | `TGFB1`, `TGFBR1/2`, `SMAD2/3`, `SERPINE1`, `CTGF` |
| **PDGF** | HSC proliferation & migration (macrophage→HSC) | `PDGFB`, `PDGFD`, `PDGFRB`, `PDGFRA` |
| **Notch** | HSC activation; endothelial→HSC angiocrine | `JAG1`, `DLL4`, `NOTCH2/3`, `HEY1`, `HES1` |
| **Hedgehog** | HSC activation, ductular reaction | `SHH`, `IHH`, `GLI1/2`, `PTCH1` |
| **YAP / TAZ (Hippo)** | mechanotransduction → myofibroblast | `YAP1`, `WWTR1`, `CTGF`, `CYR61`(`CCN1`) |
| **Wnt / β-catenin** | HSC activation; crosstalk with TGF-β | `WNT5A`, `CTNNB1`, `LGR5`, `AXIN2` |
| **TWEAK / Fn14** | macrophage→mesenchymal pro-fibrotic signal | `TNFSF12` (TWEAK), `TNFRSF12A` (Fn14) |
| **IL-1 / NF-κB / TNF** | inflammation-driven fibrogenesis | `IL1B`, `IL1R1`, `TNF`, `NFKB1`, `CCL2` |

**Key ligand–receptor axes in the fibrotic niche** (to look for in cell-cell communication):
`SPP1`(SAMac)→`CD44`/integrins(HSC) · `PDGFB`(SAMac/SAEndo)→`PDGFRB`(HSC) ·
`TNFSF12`(SAMac)→`TNFRSF12A`(HSC) · `JAG1`/`DLL4`(SAEndo)→`NOTCH3`(HSC) ·
`TGFB1`(multiple)→`TGFBR`(HSC) · `CCL2`→`CCR2` (monocyte recruitment).

---

## 4. Fibrosis-stage crosswalk (for metadata harmonization, Q1)

There is no single universal scale; map each study's scheme onto a common **ordinal axis (0–4)**
for harmonized analysis while **retaining the original label** for provenance.

| Harmonized axis | METAVIR | Ishak | NASH-CRN / Kleiner | Clinical binary | MASLD/MASH context |
|---|---|---|---|---|---|
| **0** none | F0 | 0 | 0 | non-cirrhotic / healthy | no fibrosis |
| **1** mild (portal) | F1 | 1–2 | 1 (1a/1b/1c) | non-cirrhotic | MASL or early MASH |
| **2** moderate (periportal/few septa) | F2 | 3 | 2 | non-cirrhotic | MASH, significant fibrosis |
| **3** severe (bridging septa) | F3 | 4 | 3 | non-cirrhotic | MASH, advanced fibrosis |
| **4** cirrhosis | F4 | 5–6 | 4 | **cirrhotic** | MASH cirrhosis |

Notes:
- **GSE136103 (primary)** provides only **healthy (≈0)** vs **cirrhotic (4)** → a *binary* axis;
  staged contrasts (e.g. **F2+**) come from the validation arm. State this caveat explicitly.
- **GSE244832 (validation)** spans Normal / MASL / MASH **F2–F4** → a graded axis.
- **Nomenclature (2023):** NAFLD→**MASLD**, NASH→**MASH**, "fatty liver"→**SLD** (Rinella et al.,
  multisociety Delphi consensus). Harmonize labels to MASLD/MASH and record the original term.
- "**F2+**" / "significant fibrosis" = axis ≥ 2; "**advanced fibrosis**" = axis ≥ 3.

---

## 5. Key references
- Ramachandran P, et al. *Resolving the fibrotic niche of human liver cirrhosis at single-cell level.* Nature 2019. doi:10.1038/s41586-019-1631-3 (GSE136103).
- MacParland SA, et al. *Single cell RNA sequencing of human liver reveals distinct intrahepatic macrophage populations.* Nat Commun 2018.
- Aizarani N, et al. *A human liver cell atlas reveals heterogeneity and epithelial progenitors.* Nature 2019.
- Andrews TS, et al. *Single-cell, single-nucleus, and spatial RNA sequencing of the human liver identifies cholangiocyte and mesenchymal heterogeneity.* Hepatol Commun 2022.
- Dobie R, et al. *Single-cell transcriptomics uncovers zonation of function in the mesenchyme during liver fibrosis.* Cell Rep 2019.
- Fabre T, et al. *Identification of a broadly fibrogenic macrophage subset induced by type 3 inflammation.* Sci Immunol 2023 (SPP1+/TREM2+ SAMac program; SCP2154).
- Rinella ME, et al. *A multisociety Delphi consensus statement on new fatty liver disease nomenclature.* Hepatology 2023. PMID 37363821.
- Larsen ABM, et al. *Stellate cell expression of SMOC2 is associated with NAFLD severity.* JHEP Rep 2022 (GSE207310).
