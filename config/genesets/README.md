# Gene sets

- **`fibrosis_core.gmt`** — curated, committed fibrosis gene sets (TGF-β, PDGF, ECM/collagen,
  matrix remodeling, myofibroblast activation, scar-associated macrophage/endothelium,
  Hedgehog/Notch, YAP/TAZ, Wnt, TWEAK/Fn14, inflammatory). Used by GSEA and per-cell scoring.
  Committed so results are deterministic and offline-reproducible.
- **`hallmark`** (referenced as `hallmark` in `config.pathway.gsea_sets`) — MSigDB Hallmark.
  Fetched by `gseapy` (`gseapy.get_library("MSigDB_Hallmark_2020")`) on first use and cached
  under `resources/`. For fully offline runs, drop a `h.all.*.symbols.gmt` here and point the
  config at it.

PROGENy (pathway activity) and CollecTRI (TF activity) networks are pulled by `decoupler`/`omnipath`
and cached locally on first use.
