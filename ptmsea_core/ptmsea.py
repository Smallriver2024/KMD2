from __future__ import annotations

import subprocess
from pathlib import Path


def read_flat_yaml(path: str | Path) -> dict[str, str]:
    """Read the deliberately flat public PTM-SEA parameter schema."""
    result = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if ":" not in stripped or line[:1].isspace():
            raise ValueError("ptmsea_parameters.yaml must contain flat key: value entries")
        key, value = stripped.split(":", 1)
        result[key.strip()] = value.strip().strip("'\"")
    required = {
        "sample_norm_type", "weight", "correl_type", "statistic",
        "output_score_type", "nperm", "min_overlap", "global_fdr",
        "export_signat_gct", "run_parallel",
    }
    if required - set(result):
        raise ValueError(f"PTM-SEA parameters missing: {sorted(required - set(result))}")
    return result


def invoke_ptmsea(
    gct: str | Path, gmt: str | Path, output_prefix: str | Path,
    params: str | Path, cli_script: str | Path, wrapper: str | Path,
    *, seed: int, rscript: str = "Rscript", dry_run: bool = False,
) -> list[str]:
    for path in (gct, gmt, params, cli_script, wrapper):
        if not Path(path).exists():
            raise FileNotFoundError(path)
    read_flat_yaml(params)
    if seed < 1:
        raise ValueError("selected seed must be a positive integer")
    command = [
        rscript, str(wrapper), "--cli", str(cli_script), "--gct", str(gct),
        "--gmt", str(gmt), "--output", str(output_prefix),
        "--params", str(params), "--seed", str(seed),
    ]
    if not dry_run:
        Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(command, check=True)
    return command
