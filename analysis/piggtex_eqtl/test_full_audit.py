#!/usr/bin/env python3
"""Focused tests for full-universe PigGTEx gene-eQTL archive matching."""

import csv
import gzip
import io
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from audit_full_gene_eqtl_finemap import (  # noqa: E402
    clean_gene,
    default_evidence,
    load_candidates,
    merge_evidence,
    stream_finemapped,
    stream_significant,
    unordered_variant_key,
    write_tsv,
)


def add_gz_member(archive: tarfile.TarFile, name: str, text: str) -> None:
    compressed = io.BytesIO()
    with gzip.GzipFile(fileobj=compressed, mode="wb") as handle:
        handle.write(text.encode("utf-8"))
    payload = compressed.getvalue()
    member = tarfile.TarInfo(name)
    member.size = len(payload)
    archive.addfile(member, io.BytesIO(payload))


def main() -> None:
    assert clean_gene("ENSSSCG00000000001.4") == "ENSSSCG00000000001"
    assert unordered_variant_key("chr3_18354604_T_C") == "3_18354604_C_T"
    assert unordered_variant_key("3_18354604_C_T") == "3_18354604_C_T"
    assert unordered_variant_key("not_a_variant") is None
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        candidate = root / "candidate.tsv.gz"
        fields = [
            "internal_gene_id",
            "internal_gene_name",
            "internal_tissue",
            "piggtex_gene_id",
            "piggtex_tissue",
            "internal_variant_unordered_key",
            "is_primary_tissue_enriched_pair",
        ]
        with gzip.open(candidate, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerow(
                {
                    "internal_gene_id": "ENSSSCG00000000001",
                    "internal_gene_name": "TEST1",
                    "internal_tissue": "Heart",
                    "piggtex_gene_id": "ENSSSCG00000000001",
                    "piggtex_tissue": "Heart",
                    "internal_variant_unordered_key": "1_100_A_G",
                    "is_primary_tissue_enriched_pair": "1",
                }
            )
        significant = root / "significant.tar"
        with tarfile.open(significant, "w") as archive:
            add_gz_member(
                archive,
                "Heart.cis_qtl.txt.gz",
                "phenotype_id\tvariant_id\tpval_nominal\tslope\tslope_se\n"
                "ENSSSCG00000000001.2\t1_100_G_A\t0.001\t0.4\t0.1\n",
            )
        finemapped = root / "finemapped.tar.gz"
        with tarfile.open(finemapped, "w:gz") as archive:
            add_gz_member(
                archive,
                "Heart.susieinf.txt.gz",
                "gene_id\tvariant_id\tprob\tcs\n"
                "ENSSSCG00000000001\tchr1_100_A_G\t0.75\tL1\n",
            )
        rows, exact, pairs = load_candidates(candidate)
        evidence = [default_evidence() for _ in rows]
        stream_significant(significant, exact, pairs, evidence)
        fine_evidence = [default_evidence() for _ in rows]
        stream_finemapped(finemapped, exact, pairs, fine_evidence)
        evidence = merge_evidence(evidence, fine_evidence)
        assert evidence[0]["exact_significant_eqtl"] == 1
        assert evidence[0]["maximum_exact_finemapping_pip"] == 0.75
        write_tsv(root / "test.tsv", [{"status": "PASS"}], ["status"])
        assert (root / "test.tsv").read_text(encoding="utf-8").startswith("status")
        output = root / "parallel_output"
        completed = subprocess.run(
            [
                sys.executable,
                str(SCRIPT_DIR / "audit_full_gene_eqtl_finemap.py"),
                "--candidate-input",
                str(candidate),
                "--significant-eqtl",
                str(significant),
                "--finemapped-eqtl",
                str(finemapped),
                "--out-dir",
                str(output),
                "--archive-workers",
                "2",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        assert "\"status\": \"PASS\"" in completed.stdout
        validation = json.loads((output / "validation.json").read_text(encoding="utf-8"))
        assert validation["status"] == "PASS"
    print("PASS: full-audit helper tests")


if __name__ == "__main__":
    main()
