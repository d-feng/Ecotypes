# TCGA-STAD EcoTyper vs Gastric Paper Ecotype Crosswalk

Paper reference:

- [Dissecting genetic and immune drivers of heterogeneous responses to neoadjuvant immunochemotherapy in gastric cancer](https://www.sciencedirect.com/science/article/pii/S1535610826000541)
- PubMed summary: [PMID 41720086](https://pubmed.ncbi.nlm.nih.gov/41720086/)

EcoTyper run used:

- `C:\Users\difen\Rcode\ecotyper\RecoveryOutput_STAD\bulk_stad_data`

## Important caveat

This is a soft biological mapping, not a direct label equivalence.

- The paper defines gastric-specific ecotypes `EC1-EC5` in a neoadjuvant treatment cohort.
- The EcoTyper run projects TCGA-STAD primary tumors onto the built-in pan-carcinoma `CE*` ecosystem reference.
- Therefore the crosswalk below should be interpreted as "closest recovered EcoTyper analogue" rather than an exact replication.

## Best-match crosswalk

| Paper ecotype | Paper description | Best EcoTyper match | Secondary candidates | Confidence | Rationale |
|---|---|---|---|---|---|
| `EC1` | T cell activation | `CE9` | `CE10`, `CE3` | Low-Medium | `CE9` had the strongest relative T/NK signal in the STAD run, with an immune-inflamed but still mixed myeloid context. |
| `EC2` | TLS / B cell-rich | `CE10` | `CE9`, `CE3` | Medium | `CE10` had the strongest B-cell / plasma-cell / dendritic signal and is the cleanest TLS-like immune-rich EcoTyper state. |
| `EC3` | Vascular normalization | `CE4` | `CE8`, `CE5` | Low | `CE4` had the highest relative endothelial abundance, but several STAD ecotypes were also endothelial-rich, so this match is not sharply separated. |
| `EC4` | ECM organization | `CE8` | `CE5`, `CE4` | Low | `CE8` had the strongest relative fibroblast/ECM component, with `CE5` and `CE4` close behind. |
| `EC5` | Immunosuppressive macrophage enrichment | `CE3` | `CE4`, `CE1`, `CE9` | Low-Medium | `CE3` had the highest relative monocyte/macrophage plus PMN/mast enrichment, making it the best myeloid-suppressive analogue. |

## Practical interpretation

- `CE10` looks most compatible with a lymphoid/TLS-rich gastric ecotype.
- `CE9` looks most compatible with a T-cell-inflamed immune ecotype.
- `CE4` and `CE8` look more stromal/vascular.
- `CE3` looks most compatible with macrophage-dominant suppression.

## Why the mapping is imperfect

- The publication's ecotypes are treatment-contextual and derived from gastric-specific multi-omic analyses.
- EcoTyper `CE*` labels are conserved pan-carcinoma ecosystem classes.
- Several STAD EcoTyper ecotypes are mixed states with simultaneous endothelial, fibroblast, epithelial, and immune enrichment.

## Next step for a stronger comparison

To improve this crosswalk, derive sample-level scores for the paper's gastric ecotype marker programs and compare them directly to EcoTyper `CE*` assignments in TCGA-STAD.
