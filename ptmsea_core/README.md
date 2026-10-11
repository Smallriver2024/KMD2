# PTM-SEA minimal call notes

This folder contains the selected PTM-SEA core functions plus the minimal
shared I/O and resource-validation helpers required by their relative imports.
`config/ptmsea_parameters.yaml` defines the public PTM-SEA parameter contract;
`config/release_config.yaml` defines the release fields and gates.

The usual flow is: prepare site rows (`sites.py`) -> build and write GCT
(`gct.py`) -> validate resources and call external ssGSEA2.0/PTM-SEA through
`invoke_ptmsea` (`ptmsea.py`) -> calculate statistics (`stats.py`) -> attach
site and PPI evidence (`evidence.py`, `network.py`) -> assemble release rows
(`release.py`). The external PTM-SEA CLI, GMT resources, source data and full
study configuration are intentionally not included.
