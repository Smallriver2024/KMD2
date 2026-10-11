# KinaseMD 2.0 key functions

This repository is a compact collection of the main analysis functions used in KinaseMD 2.0.

## Included code

| Area | Files | Key functions |
| --- | --- | --- |
| PTM-SEA site and scoring logic | `ptmsea_core/sites.py`, `gct.py`, `ptmsea.py`, `stats.py`, `io.py`, `resources.py` | `prepare_sites`, `build_site_matrix`, `write_gct`, `invoke_ptmsea`, `compute_statistics`, `bh_adjust`, `read_tsv`, `validate_execution` |
| PTM-SEA evidence and release logic | `ptmsea_core/evidence.py`, `network.py`, `release.py`, `config.py`, `config/` | `build_site_evidence`, `annotate_network_support`, `assemble_phosphoproteomics_release`, `sample_groups`; minimal parameter and release contracts |
| Literature curation | `literature_curation/` | `normalize_kinase_field`, `normalize_target`, `normalize_agent`, `normalize_disease`, `apply_entity_level_edits`, `apply_final_target`, `classify_final` |
| Candidate drug resistance | `resistance/compute_drkinase_resistance.py` | `map_mutations_to_predicted_hotspots`, `load_response`, `build_cell_candidates`, `attach_mutation_evidence`, `build_primary_secondary` |
| Revision resistance statistics | `resistance_revision/update_resistance_statistics.py` | matched kinase/drug/cancer/sub-tissue grouping, independent-cell IC50 aggregation, two-sided Mann–Whitney U tests, comparison-wise BH-FDR |
| Drug and disease normalization | `resistance_revision/normalize_web.py`, `map_drugs.py`, `map_diseases.py` | conservative intervention classification and PubChem matching; exact ontology/synonym disease matching with reviewable unresolved cases |
| Mutation and drug-response inputs | `mutation_drug_response/` | `parse_domain_ranges`, `build_known_domains`, `build_predicted_domains`, `filter_non_overlapping_predictions`, `resolve_ccle_targets`, `parse_gdsc_targets`, `process_gdsc_excel` |
| GPS5.0 network rewiring | `network_rewiring/gps5_rewiring_core.py` | `classify_rewiring`, `classify_gps_output`, `deduplicate_rewiring_events` |
| Survival analysis | `survival_analysis/km_logrank_core.py`, `survival_matrix_core.py`, `render_km_ggsurvplot.R` | `split_lower_upper_half`, `status_to_event`, `logrank_test`, `analyse_expression_survival`, `parse_logrank_records`, `read_group_counts`, `make_result_row`, `plot_one` |

The mutation and drug-response modules illustrate the analytical chain from kinase substructure annotation and mutation mapping to harmonized IC50 measurements and candidate resistance associations.

The revision resistance-statistics script is retained because it is the actual analysis used to add group sizes, two-sided Mann–Whitney U p-values, and Benjamini–Hochberg FDR values to the revised download. It tests all eligible exact-variant comparisons within matched kinase/drug/cancer/sub-tissue strata before candidate selection. The normalization scripts preserve the conservative matching rules used for this release: drug/agent names are matched only through supported PubChem identifiers or retained for review, and diseases use exact ontology labels/synonyms or remain explicitly unresolved.

The network-rewiring module applies GPS5.0 reference and mutant scores with kinase-specific cutoffs to classify Gain, Loss, Up-regulation, and Down-regulation events, then deduplicates events. GPS5.0 performs the upstream sequence scoring. These categories represent predicted phosphorylation propensity, not measured phosphorylation abundance or calibrated effect sizes.

### GPS5.0 reference

The sequence-scoring tool used upstream of the rewiring classification is GPS 5.0. Please cite: Wang C, Xu H, Lin S, Deng W, Zhou J, Zhang Y, Shi Y, Peng D, and Xue Y. *GPS 5.0: An Update on the Prediction of Kinase-Specific Phosphorylation Sites in Proteins.* **Genomics, Proteomics & Bioinformatics** 18(1), 72–80 (2020). [https://doi.org/10.1016/j.gpb.2020.01.001](https://doi.org/10.1016/j.gpb.2020.01.001).

### PTM-SEA reference

The PTM-SEA analysis uses the official [ssGSEA2.0/PTM-SEA repository](https://github.com/broadinstitute/ssGSEA2.0). Please cite: Krug K, Mertins P, Zhang B, et al. *A Curated Resource for Phosphosite-specific Signature Analysis.* **Molecular & Cellular Proteomics** 18(3), 576–593 (2019). [https://doi.org/10.1074/mcp.TIR118.000943](https://doi.org/10.1074/mcp.TIR118.000943).

The survival module provides lower/upper expression grouping, two-group log-rank testing, result summarization, and Kaplan–Meier plotting with confidence intervals, median-survival guides, and risk tables. No clinical or expression data are included.

## PTM literature-extraction prompts

The LLM prompts and two-phase extraction/review design are maintained separately in the previously published [PhosphoRLIE skill](https://github.com/Smallriver2024/phosphorlie-skill). Use that repository as the authoritative prompt resource rather than duplicating prompt text here.

PhosphoRLIE provides an Extractor + Reviewer workflow for PTM events from biomedical literature, including the standalone extractor and reviewer prompts, JSON schema, and examples. For the exact prompt definitions, see its [prompts directory](https://github.com/Smallriver2024/phosphorlie-skill/tree/main/prompts) and [SKILL.md](https://github.com/Smallriver2024/phosphorlie-skill/blob/main/SKILL.md).
