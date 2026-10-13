"""Read explicit paths and sample identities for ASM analysis stages."""

import hashlib
import json
import os
from pathlib import Path
import re
import sys


def config_path(required=True):
    value = os.environ.get("MINIPIG_ANALYSIS_CONFIG")
    for i, arg in enumerate(sys.argv[1:], 1):
        if arg == "--config" and i + 1 < len(sys.argv):
            value = sys.argv[i + 1]
        elif arg.startswith("--config="):
            value = arg.split("=", 1)[1]
    if not value:
        if required:
            raise ValueError("Supply --config analysis_inputs.json")
        return None
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        # Help must remain usable before inputs have been prepared.
        if not required and any(x in sys.argv for x in ("-h", "--help")):
            return None
        raise ValueError("Analysis configuration is missing: " + str(path))
    os.environ["MINIPIG_ANALYSIS_CONFIG"] = str(path)
    return path


def read_config(required=True):
    path = config_path(required)
    return json.loads(path.read_text(encoding="utf-8")) if path else {}


def read_stage_config():
    """Resolve native controller paths relative to the supplied JSON file."""
    base = config_path().parent
    path_fields = {
        "path",
        "gate",
        "index_path",
        "stage1_root",
        "stage2_root",
    }

    def absolute(value):
        return str((base / Path(value).expanduser()).resolve())

    def visit(value, key=""):
        if isinstance(value, dict):
            if key in {"registry", "checks", "root_checks"}:
                return {absolute(k): visit(v) for k, v in value.items()}
            if key == "locations":
                return {
                    stage: {name: absolute(path) for name, path in entries.items()}
                    for stage, entries in value.items()
                }
            return {k: visit(v, k) for k, v in value.items()}
        if isinstance(value, list):
            return [visit(v) for v in value]
        return absolute(value) if key in path_fields and isinstance(value, str) else value

    config = visit(read_config())
    for name in config.get("assets", {}):
        path = configured_asset(name)
        with path.open("rb") as handle:
            checksum = hashlib.file_digest(handle, "sha256").hexdigest()
        config["assets"][name] = {"path": str(path), "sha256": checksum}
        if "root_checks" in config:
            config["root_checks"][str(path)] = checksum
    return config


def configured_samples():
    cfg = read_config(False)
    ids = cfg.get(
        "ids", cfg.get("source_inputs", {}).get("ids", [f"S{i:02d}" for i in range(1, 11)])
    )
    if (
        len(ids) != 10
        or len(set(ids)) != 10
        or any(
            not isinstance(s, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", s)
            for s in ids
        )
    ):
        raise ValueError("ASM stages require ten unique sample IDs in a fixed order")
    return list(ids)


def has_asset(name):
    return name in read_config(False).get("assets", {})


def configured_asset(name):
    cfg = read_config()
    if name not in cfg.get("assets", {}):
        raise ValueError("Missing input file in configuration assets: " + name)
    spec = cfg["assets"][name]
    if isinstance(spec, str):
        spec = {"path": spec}
    path = Path(spec["path"]).expanduser()
    path = (config_path().parent / path).resolve() if not path.is_absolute() else path.resolve()
    if not path.is_file():
        raise ValueError("Input file is missing: " + str(path))
    if spec.get("sha256"):
        with path.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        if digest != spec["sha256"]:
            raise ValueError("Configured input checksum differs: " + name)
    return path
