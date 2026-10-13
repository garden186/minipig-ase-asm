"""Input and checkpoint validation for fragment counting."""

import re
import shutil
import sys
from pathlib import Path
from common import CONTRACT, ROLES, fingerprint, sha, read_json, write_json, capture


def index_for(path, bam=False):
    candidates = [path + s for s in ((".bai", ".csi") if bam else (".tbi", ".csi"))]
    if bam:
        candidates += [str(Path(path).with_suffix(".bai"))]
    for p in candidates:
        if Path(p).is_file():
            return p
    raise ValueError("No existing index found for " + path)


def preflight(data, args, output):
    missing = []
    for row in data:
        for role in ROLES:
            if not Path(row[role]).is_file():
                missing.append(row["sample"] + "-" + row["tissue"] + " " + role + ": " + row[role])
    if missing:
        write_json(
            output / "preflight.json",
            {"status": "FAIL_PREFLIGHT", "errors": missing, "bam_contents_read": False},
        )
        raise ValueError("Required inputs are missing; inspect preflight.json")
    tools, versions = {}, {}
    for name in ("samtools", "bcftools"):
        executable = shutil.which(getattr(args, name))
        if not executable:
            raise ValueError("Required tool is unavailable: " + name)
        tools[name] = str(Path(executable).resolve())
        versions[name] = capture(
            [executable, "--version"], output / (name + ".stderr.log"), 30
        ).splitlines()[:4]
    versions.update(python=sys.version, executable=sys.executable)
    reference_paths = {r["reference"] for r in data}
    if len(reference_paths) != 1:
        raise ValueError("A single explicit matching reference FASTA is required")
    reference = next(iter(reference_paths))
    lengths = {}
    with open(reference + ".fai", encoding="utf-8") as f:
        for line in f:
            a = line.split("\t")
            lengths[a[0]] = int(a[1])
    if any(str(c) not in lengths for c in range(1, 19)):
        raise ValueError(
            "Reference must contain the recorded Sscrofa11.1 autosomes named 1 through 18"
        )
    units, seen_headers, errors = [], set(), []
    for row in data:
        unit = row["sample"] + "-" + row["tissue"]
        try:
            paths = {role: row[role] for role in ROLES}
            for role in ("bam", "phased_vcf", "genotype_vcf"):
                paths[role + "_index"] = index_for(row[role], role == "bam")
            paths["reference_fai"] = reference + ".fai"
            inputs = {role: fingerprint(path) for role, path in paths.items()}
            if inputs["bam"]["size_bytes"] != int(row["expected_bam_bytes"]):
                raise ValueError("BAM size differs from the recorded source")
            capture(
                [tools["samtools"], "quickcheck", "-v", row["bam"]],
                output / (unit + ".quickcheck.log"),
            )
            header = capture(
                [tools["samtools"], "view", "-H", row["bam"]], output / (unit + ".header.log")
            )
            (output / (unit + ".bam_header.sam")).write_text(header, encoding="utf-8")
            dictionary, sorted_bam = {}, False
            for line in header.splitlines():
                fields = line.split("\t")
                tags = dict(x.split(":", 1) for x in fields[1:] if ":" in x)
                if fields[0] == "@HD":
                    sorted_bam = tags.get("SO") == "coordinate"
                if fields[0] == "@SQ":
                    dictionary[tags["SN"]] = int(tags["LN"])
            if not sorted_bam or any(
                dictionary.get(str(c)) != lengths[str(c)] for c in range(1, 19)
            ):
                raise ValueError("BAM sorting or reference chromosome lengths disagree")
            for role in ("phased_vcf", "genotype_vcf"):
                key = (role, row[role], row["sample"])
                if key in seen_headers:
                    continue
                header = capture(
                    [tools["bcftools"], "view", "-h", row[role]],
                    output / (row["sample"] + "." + role + ".log"),
                )
                sample_names = next(
                    line.split("\t")[9:]
                    for line in header.splitlines()
                    if line.startswith("#CHROM")
                )
                vcf_lengths = {
                    c: int(n) for c, n in re.findall(r"##contig=<ID=([^,>]+),length=(\d+)", header)
                }
                if sample_names != [row["sample"]] or any(
                    vcf_lengths.get(str(c)) != lengths[str(c)] for c in range(1, 19)
                ):
                    raise ValueError("VCF sample or reference chromosome lengths disagree: " + role)
                seen_headers.add(key)
            units.append(
                {
                    "sample": row["sample"],
                    "tissue": row["tissue"],
                    "inputs": inputs,
                    "manifest_paths": paths,
                    "bam_device": Path(inputs["bam"]["path"]).stat().st_dev,
                }
            )
        except Exception as exc:
            errors.append(unit + ": " + str(exc))
    result = {
        "status": "PASS_PREFLIGHT" if not errors else "FAIL_PREFLIGHT",
        "errors": errors,
        "units": units,
        "tools": tools,
        "versions": versions,
        "lengths": {str(c): lengths[str(c)] for c in range(1, 19)},
        "declared_assembly": CONTRACT["assembly"],
        "genetic_identity_verified": False,
        "reference_note": "Recorded reference path and matching chromosome dictionaries; not a whole-FASTA checksum or sample identity test",
    }
    write_json(output / "preflight.json", result)
    if errors:
        raise ValueError("Input preflight failed; inspect preflight.json")
    return result


def complete_region(path, task_id):
    p = Path(path)
    if not (p / "validation.json").exists() or not (p / "timing.json").exists():
        return False
    v = read_json(p / "validation.json")
    if v["status"] != "PASS_DESCRIPTIVE_CPG_ASSIGNMENT" or v["errors"]:
        return False
    if read_json(p / "task.json")["task_id"] != task_id:
        raise ValueError("Checkpoint task identity differs")
    if any(sha(p / n) != s for n, s in v["sha256"].items()):
        raise ValueError("A validated checkpoint was modified")
    return True


def main():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from count import main as count_main

    count_main()


if __name__ == "__main__":
    main()
