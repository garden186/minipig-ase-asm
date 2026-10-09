#!/usr/bin/env python3
"""WGS SNP balance and phase-orientation helpers for ASE v0.2.5.

The phased WGS VCF supplies genotype eligibility and phase.  REF/ALT balance
is evaluated only from the phASER-WGS allelic-count table generated from the
same high-confidence WGS BAM used by the phASER workflow.
"""

from __future__ import print_function

import csv
import gzip
from collections import Counter


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_variant_id(variant_id):
    parts = str(variant_id).rsplit("_", 3)
    if len(parts) != 4:
        return None
    contig, position, ref, alt = parts
    try:
        position = int(position)
    except ValueError:
        return None
    return contig, position, ref.upper(), alt.upper()


def _float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def is_biallelic_het(gt):
    alleles = str(gt or "").replace("|", "/").split("/")
    return len(alleles) == 2 and set(alleles) == {"0", "1"}


def load_phaser_wgs_counts(path, requested):
    """Return requested variantID rows from phASER-WGS allelic counts."""
    found = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "contig", "position", "variantID", "refAllele", "altAllele",
            "refCount", "altCount", "totalCount",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "phASER-WGS count table missing columns: {}".format(
                    ",".join(sorted(missing))
                )
            )
        for row in reader:
            variant_id = row["variantID"]
            if variant_id not in requested:
                continue
            parsed = requested[variant_id]
            observed = (
                row["contig"],
                int(row["position"]),
                row["refAllele"].upper(),
                row["altAllele"].upper(),
            )
            if observed != parsed:
                raise ValueError(
                    "phASER-WGS allele mismatch for {}: expected={} observed={}"
                    .format(variant_id, parsed, observed)
                )
            if variant_id in found:
                raise ValueError(
                    "Duplicate phASER-WGS row: {}".format(variant_id)
                )
            found[variant_id] = dict(row)
    return found


def load_phased_wgs_vcf(path, sample, requested):
    """Read requested biallelic SNPs, phased GT and PS from a sample VCF."""
    found = {}
    sample_index = None
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]
                if sample not in samples:
                    raise ValueError(
                        "Sample {} not found in phased WGS VCF".format(sample)
                    )
                sample_index = 9 + samples.index(sample)
                continue
            if line.startswith("#"):
                continue
            if sample_index is None:
                raise ValueError("VCF #CHROM header was not found")
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= sample_index:
                continue
            chrom, pos, _, ref, alts, _, filt, _, fmt = fields[:9]
            alt_values = alts.split(",")
            candidate_ids = [
                "{}_{}_{}_{}".format(chrom, pos, ref, alt)
                for alt in alt_values
            ]
            matching = [
                (index, variant_id)
                for index, variant_id in enumerate(candidate_ids, start=1)
                if variant_id in requested
            ]
            if not matching:
                continue
            keys = fmt.split(":")
            values = fields[sample_index].split(":")
            sample_data = dict(zip(keys, values))
            original_gt = sample_data.get("GT", "")
            for alt_index, variant_id in matching:
                gt = original_gt
                if len(alt_values) > 1:
                    separator = "|" if "|" in gt else "/"
                    alleles = gt.replace("|", "/").split("/")
                    if set(alleles) == {"0", str(alt_index)}:
                        gt = separator.join(
                            "1" if allele == str(alt_index) else "0"
                            for allele in alleles
                        )
                row = {
                    "variant_id": variant_id,
                    "contig": chrom,
                    "position": int(pos),
                    "ref": ref.upper(),
                    "alt": alt_values[alt_index - 1].upper(),
                    "filter": filt,
                    "gt": gt,
                    "gq": _float(sample_data.get("GQ")),
                    "ps": sample_data.get("PS", ""),
                }
                if variant_id in found:
                    raise ValueError(
                        "Duplicate phased WGS VCF row: {}".format(variant_id)
                    )
                found[variant_id] = row
    return found


def evaluate_wgs_snps(
    requested,
    vcf_rows,
    count_rows,
    min_gq,
    min_total,
    min_each,
    min_alt_fraction,
    max_alt_fraction,
):
    """Return variant-level WGS QC rows and the balanced variant set."""
    rows = []
    balanced = set()
    for variant_id in sorted(
        requested,
        key=lambda item: (
            requested[item][0], requested[item][1], item
        ),
    ):
        contig, position, ref, alt = requested[variant_id]
        vcf = vcf_rows.get(variant_id)
        counts = count_rows.get(variant_id)
        reasons = []
        is_snv = len(ref) == 1 and len(alt) == 1
        if not is_snv:
            reasons.append("NOT_SNV")
        if vcf is None:
            reasons.append("VCF_MISSING")
        else:
            if vcf["filter"] != "PASS":
                reasons.append("VCF_FILTER_FAIL")
            if not is_biallelic_het(vcf["gt"]):
                reasons.append("NOT_BIALLELIC_HET")
            if vcf["gq"] is None or vcf["gq"] < min_gq:
                reasons.append("GQ_LOW")

        ref_count = _int(counts.get("refCount")) if counts else None
        alt_count = _int(counts.get("altCount")) if counts else None
        total_count = _int(counts.get("totalCount")) if counts else None
        alt_fraction = (
            alt_count / float(total_count)
            if alt_count is not None and total_count is not None
            and total_count > 0 else None
        )
        genotype_ok = (
            is_snv
            and
            vcf is not None
            and vcf["filter"] == "PASS"
            and is_biallelic_het(vcf["gt"])
            and vcf["gq"] is not None
            and vcf["gq"] >= min_gq
        )
        if counts is None:
            reasons.append("PHASER_WGS_MISSING")
        else:
            if total_count is None or total_count < min_total:
                reasons.append("WGS_TOTAL_LOW")
            if (
                ref_count is None or alt_count is None
                or min(ref_count, alt_count) < min_each
            ):
                reasons.append("WGS_ALLELE_COUNT_LOW")
            if (
                alt_fraction is None
                or alt_fraction < min_alt_fraction
                or alt_fraction > max_alt_fraction
            ):
                reasons.append("WGS_ALT_FRACTION_IMBALANCED")
        coverage_ok = (
            counts is not None
            and total_count is not None and total_count >= min_total
            and ref_count is not None and alt_count is not None
            and min(ref_count, alt_count) >= min_each
        )
        balance_ok = (
            alt_fraction is not None
            and min_alt_fraction <= alt_fraction <= max_alt_fraction
        )
        if genotype_ok and coverage_ok and balance_ok:
            status = "WGS_BALANCED"
            balanced.add(variant_id)
        elif genotype_ok and coverage_ok:
            status = "WGS_UNBALANCED"
        else:
            status = "WGS_INSUFFICIENT"
        rows.append({
            "variant_id": variant_id,
            "contig": contig,
            "position": position,
            "ref": ref,
            "alt": alt,
            "vcf_filter": vcf["filter"] if vcf else "",
            "vcf_gt": vcf["gt"] if vcf else "",
            "vcf_gq": vcf["gq"] if vcf else None,
            "vcf_phase_set": vcf["ps"] if vcf else "",
            "wgs_ref_count": ref_count,
            "wgs_alt_count": alt_count,
            "wgs_total_count": total_count,
            "wgs_alt_fraction": alt_fraction,
            "wgs_snp_status": status,
            "wgs_snp_reason": ",".join(reasons) if reasons else "PASS",
        })
    return rows, balanced


def phase_agreement(component, variant_ids, vcf_rows):
    """Map phASER Hap A/B to one WGS phase-set orientation."""
    by_variant = {
        row["variant_id"]: row for row in component["parsed_variants"]
    }
    details = []
    phase_sets = set()
    for variant_id in sorted(set(variant_ids)):
        parsed = by_variant.get(variant_id)
        vcf = vcf_rows.get(variant_id)
        if parsed is None or vcf is None or "|" not in str(vcf["gt"]):
            continue
        alleles = str(vcf["gt"]).split("|")
        if len(alleles) != 2 or not set(alleles).issubset({"0", "1"}):
            continue
        lookup = {"0": vcf["ref"], "1": vcf["alt"]}
        left, right = lookup[alleles[0]], lookup[alleles[1]]
        if parsed["hapA"] == left and parsed["hapB"] == right:
            orientation = "DIRECT"
        elif parsed["hapA"] == right and parsed["hapB"] == left:
            orientation = "FLIPPED"
        else:
            orientation = "ALLELE_MISMATCH"
        if vcf.get("ps") not in ("", ".", None):
            phase_sets.add(str(vcf["ps"]))
        details.append({
            "variant_id": variant_id,
            "orientation": orientation,
            "phase_set": str(vcf.get("ps") or ""),
        })

    counts = Counter(row["orientation"] for row in details)
    comparable = counts["DIRECT"] + counts["FLIPPED"]
    if counts["ALLELE_MISMATCH"]:
        status = "ALLELE_MISMATCH"
        orientation = "UNRESOLVED"
    elif comparable < 2:
        status = "INSUFFICIENT"
        orientation = "UNRESOLVED"
    elif len(phase_sets) != 1:
        status = "MULTI_PHASE_SET"
        orientation = "UNRESOLVED"
    elif counts["DIRECT"] == comparable:
        status = "CONCORDANT"
        orientation = "DIRECT"
    elif counts["FLIPPED"] == comparable:
        status = "CONCORDANT"
        orientation = "FLIPPED"
    else:
        status = "DISCORDANT"
        orientation = "MIXED"
    return {
        "wgs_phase_status": status,
        "wgs_phase_orientation": orientation,
        "wgs_phase_set": ",".join(sorted(phase_sets)),
        "n_wgs_phase_compared_snps": comparable,
        "n_wgs_phase_direct_snps": counts["DIRECT"],
        "n_wgs_phase_flipped_snps": counts["FLIPPED"],
        "n_wgs_phase_mismatch_snps": counts["ALLELE_MISMATCH"],
    }
