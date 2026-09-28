#!/usr/bin/env python3
from __future__ import print_function

import csv
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest


HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "call_phaser_ase_candidates.py")
SPEC = importlib.util.spec_from_file_location("caller_v023", SCRIPT)
CALLER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CALLER)
PROJECT_TMP = os.path.abspath(
    os.environ.get("ASE_PROJECT_TMP", os.path.join(HERE, "tmp"))
)
os.makedirs(PROJECT_TMP, exist_ok=True)


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


class StandaloneCandidateIntegrationTest(unittest.TestCase):
    def test_extreme_binomial_pvalue_and_bh_do_not_underflow_to_zero(self):
        pvalue = CALLER.exact_binom_two_sided_half(6, 6666)
        self.assertEqual(pvalue, sys.float_info.min)
        self.assertGreater(pvalue, 0.0)
        rows = [
            {
                "measurement_id": "M_UNDERFLOW",
                "sample": "S1", "tissue": "Lung",
                "analysis_set": "MAIN", "ase_test_status": "TESTED",
                "ase_a_count": 6, "ase_b_count": 6666,
                "ase_p_exact": pvalue, "ase_q_bh": None
            },
            {
                "measurement_id": "M_REGULAR",
                "sample": "S1", "tissue": "Lung",
                "analysis_set": "MAIN", "ase_test_status": "TESTED",
                "ase_a_count": 20, "ase_b_count": 5,
                "ase_p_exact": 0.01, "ase_q_bh": None
            }
        ]
        CALLER.apply_unique_measurement_bh(rows)
        self.assertGreater(rows[0]["ase_q_bh"], 0.0)
        self.assertAlmostEqual(
            rows[0]["ase_q_bh"], sys.float_info.min * 2
        )

    def test_unique_measurement_fdr_includes_low_minor_and_extreme(self):
        rows = []
        for measurement_id, pvalue, counts in (
            ("M1", 0.01, (20, 5)),
            ("M1", 0.01, (20, 5)),
            ("M2", 0.02, (20, 0)),
            ("M3", 0.03, (20, 2))
        ):
            rows.append({
                "measurement_id": measurement_id,
                "sample": "S1", "tissue": "Heart",
                "analysis_set": "MAIN", "ase_test_status": "TESTED",
                "ase_a_count": counts[0], "ase_b_count": counts[1],
                "ase_p_exact": pvalue, "ase_q_bh": None
            })
        CALLER.apply_unique_measurement_bh(rows)
        self.assertEqual(rows[0]["measurement_annotation_multiplicity"], 2)
        self.assertEqual(rows[1]["measurement_annotation_multiplicity"], 2)
        self.assertAlmostEqual(rows[0]["ase_q_bh"], 0.03)
        self.assertAlmostEqual(rows[1]["ase_q_bh"], 0.03)
        self.assertAlmostEqual(rows[2]["ase_q_bh"], 0.03)
        self.assertAlmostEqual(rows[3]["ase_q_bh"], 0.03)

    def test_exon_union_unique_ambiguous_and_intronic_assignment(self):
        genes = {
            "G1": {
                "gene_id": "G1", "gene_name": "G1",
                "gene_biotype": "protein_coding", "contig": "1",
                "start": 50, "end": 350, "strand": "+",
                "exon_intervals": [(90, 210), (300, 320)]
            },
            "G2": {
                "gene_id": "G2", "gene_name": "G2",
                "gene_biotype": "protein_coding", "contig": "1",
                "start": 180, "end": 280, "strand": "-",
                "exon_intervals": [(190, 230)]
            }
        }
        variants = (
            "1_100_A_G", "1_200_C_T", "1_250_G_A", "1_310_T_C"
        )
        components = {
            "C1": {
                "component_id": "C1", "identity": ("1", 100, 310, variants),
                "contig": "1", "start": 100, "end": 310,
                "phaser_pi": "1", "variants": variants,
                "variant_set": set(variants), "haplotypeA": "",
                "haplotypeB": "", "genes": set(), "gene_variants": {},
                "n_bam_rows": 1
            }
        }
        mapped, assignment, invalid = CALLER.map_components_to_genes(
            components, genes
        )
        self.assertEqual(invalid, [])
        self.assertEqual(
            mapped["C1"]["gene_variants"]["G1"],
            ["1_100_A_G", "1_310_T_C"]
        )
        self.assertEqual(
            assignment["exonic"]["1_200_C_T"], ["G1", "G2"]
        )
        self.assertEqual(assignment["span"]["1_250_G_A"], ["G1", "G2"])

    def test_ambiguous_only_component_is_retained_for_strand_resolution(self):
        genes = {
            "GP": {
                "gene_id": "GP", "gene_name": "GP",
                "gene_biotype": "protein_coding", "contig": "1",
                "start": 90, "end": 110, "strand": "+",
                "exon_intervals": [(90, 110)]
            },
            "GM": {
                "gene_id": "GM", "gene_name": "GM",
                "gene_biotype": "protein_coding", "contig": "1",
                "start": 90, "end": 110, "strand": "-",
                "exon_intervals": [(90, 110)]
            }
        }
        variants = ("1_100_A_G",)
        components = {
            "C1": {
                "component_id": "C1",
                "identity": ("1", 100, 100, variants),
                "contig": "1", "start": 100, "end": 100,
                "phaser_pi": "", "variants": variants,
                "variant_set": set(variants), "haplotypeA": "A",
                "haplotypeB": "G", "genes": set(),
                "gene_variants": {}, "n_bam_rows": 1
            }
        }
        mapped, assignment, invalid = CALLER.map_components_to_genes(
            components, genes
        )
        self.assertEqual(invalid, [])
        self.assertIn("C1", mapped)
        self.assertEqual(mapped["C1"]["gene_variants"], {})
        self.assertEqual(
            set(mapped["C1"]["candidate_gene_variants"]),
            {"GP", "GM"}
        )
        self.assertEqual(
            assignment["exonic"]["1_100_A_G"], ["GM", "GP"]
        )

    def test_rescue_evidence_on_haplotype_conflict_is_not_accepted(self):
        fragment = CALLER.bam_engine.new_fragment_record()
        fragment["hap_support"].update(("A", "B"))
        fragment["transcript_strands"].add("-")
        fragment["gene_evidence"].add("GM")
        fragment["strand_resolved_gene_evidence"].add("GM")
        fragment["gene_variant_evidence"]["GM"].add("1_100_A_G")
        fragment["supported_variants"].add("1_100_A_G")
        measurement = {
            "sample": "S1", "tissue": "Heart", "bam": "S1_Heart"
        }
        component = {"component_id": "C1"}
        genes = {"GM": {"gene_name": "GM", "strand": "-"}}
        row, _ = CALLER.bam_engine.finalize_fragment_assignment(
            "q1", fragment, measurement, component, genes
        )
        self.assertEqual(row["final_status"], "HAPLOTYPE_CONFLICT")
        self.assertEqual(row["gene_assignment_method"], "STRAND_RESOLVED")
        self.assertEqual(row["strand_resolution_required"], 0)

    def legacy_v023_gtf_to_candidate_without_component_master(self):
        with tempfile.TemporaryDirectory(
            prefix="phaser_ase_v023_", dir=PROJECT_TMP
        ) as temp:
            gtf = os.path.join(temp, "test.gtf")
            haplotypes = os.path.join(temp, "haplotypes.tsv")
            connections = os.path.join(temp, "connections.tsv")
            allele_config = os.path.join(temp, "allele_config.tsv")
            counts = os.path.join(temp, "counts.tsv")
            out_dir = os.path.join(temp, "out")

            with open(gtf, "w") as handle:
                handle.write(
                    '1\ttest\tgene\t50\t350\t.\t+\t.\t'
                    'gene_id "G1"; gene_name "GENE1"; '
                    'gene_biotype "protein_coding";\n'
                )
                handle.write(
                    '1\ttest\ttranscript\t50\t350\t.\t+\t.\t'
                    'gene_id "G1"; transcript_id "T1"; '
                    'gene_name "GENE1"; gene_biotype "protein_coding"; '
                    'transcript_biotype "protein_coding"; '
                    'tag "Ensembl_canonical";\n'
                )
                handle.write(
                    '1\ttest\texon\t50\t350\t.\t+\t.\t'
                    'gene_id "G1"; transcript_id "T1"; '
                    'gene_name "GENE1"; gene_biotype "protein_coding"; '
                    'transcript_biotype "protein_coding";\n'
                )
                handle.write(
                    '1\ttest\tCDS\t100\t300\t.\t+\t0\t'
                    'gene_id "G1"; transcript_id "T1"; '
                    'gene_name "GENE1"; gene_biotype "protein_coding"; '
                    'transcript_biotype "protein_coding";\n'
                )

            write_tsv(haplotypes, [
                "contig", "start", "stop", "length", "variants",
                "variant_ids", "variant_alleles", "reads_hap_a",
                "reads_hap_b", "reads_total", "edges_supporting",
                "edges_total", "annotated_phase", "phase_concordant",
                "gw_phase", "gw_confidence"
            ], [[
                "1", "100", "300", "200", "3",
                "1_100_A_G,1_200_C_T,1_300_G_A", "A,C,G|G,T,A",
                "20", "5", "25", "3", "3", "000|111", "1", "000|111", "1"
            ]])
            write_tsv(connections, [
                "variant_a", "variant_b", "supporting_connections",
                "total_connections", "conflicting_configuration_p",
                "phase_concordant"
            ], [
                ["1_100_A_G", "1_200_C_T", "4", "4", "1", "1"],
                ["1_200_C_T", "1_300_G_A", "4", "4", "1", "1"],
                ["1_100_A_G", "1_300_G_A", "3", "3", "1", "1"]
            ])
            write_tsv(allele_config, [
                "variant_a", "rsid_a", "variant_b", "rsid_b",
                "configuration"
            ], [
                [
                    "1_100_A_G", "1_100_A_G", "1_200_C_T",
                    "1_200_C_T", "cis"
                ],
                ["MALFORMED_ONLY_ONE_FIELD"]
            ])
            counts_header = [
                "contig", "start", "stop", "variants", "variantCount",
                "variantsBlacklisted", "variantCountBlacklisted",
                "haplotypeA", "haplotypeB", "aCount", "bCount",
                "totalCount", "blockGWPhase", "gwStat", "max_haplo_maf",
                "bam", "aReads", "bReads"
            ]
            a_reads = [
                ",".join(str(i) for i in range(0, 10)),
                ",".join(str(i) for i in range(5, 15)),
                ",".join(str(i) for i in range(10, 20))
            ]
            b_reads = ["0,1", "1,2,3", "3,4"]
            rows = []
            for tissue in ("Heart", "Cranial"):
                rows.append([
                    "1", "100", "300",
                    "1_100_A_G,1_200_C_T,1_300_G_A", "3", "", "0",
                    "A,C,G", "G,T,A", "20", "5", "25", "0|1", "1", "0",
                    "S1_{}.phaseready".format(tissue),
                    ";".join(a_reads), ";".join(b_reads)
                ])
            write_tsv(counts, counts_header, rows)

            subprocess.check_call([
                sys.executable, SCRIPT,
                "--sample", "S1",
                "--gtf", gtf,
                "--variant-connections", connections,
                "--haplotypes", haplotypes,
                "--haplotypic-counts", counts,
                "--full-allele-config", allele_config,
                "--write-transcript-annotation", "1",
                "--out-dir", out_dir,
                "--force"
            ])

            with open(
                os.path.join(out_dir, "ase_candidate_individual.tsv")
            ) as handle:
                candidates = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(candidates), 2)
            heart = [x for x in candidates if x["tissue"] == "Heart"][0]
            cranial = [x for x in candidates if x["tissue"] == "Cranial"][0]
            self.assertEqual(heart["gene_id"], "G1")
            self.assertEqual(heart["count_validation"], "MATCH")
            self.assertEqual(heart["phase_qc"], "PHASE_HC")
            self.assertEqual(
                heart["candidate_class"], "TIER1A_HC_LOO_ROBUST"
            )
            self.assertEqual(cranial["candidate_class"], "EXCLUDED_TISSUE")

            with open(
                os.path.join(out_dir, "phaser_component_gene_map.tsv")
            ) as handle:
                mapping = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(mapping), 1)
            self.assertEqual(mapping[0]["gene_id"], "G1")
            self.assertEqual(mapping[0]["n_gene_variants"], "3")
            self.assertEqual(
                mapping[0]["gene_assignment_policy"],
                "UNIQUE_PROTEIN_CODING_EXON_UNION"
            )

            with open(
                os.path.join(out_dir, "ase_measurement_unique.tsv")
            ) as handle:
                measurements = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(measurements), 2)

            with open(
                os.path.join(
                    out_dir, "ase_variant_transcript_annotation.tsv"
                )
            ) as handle:
                annotations = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(annotations), 3)
            self.assertTrue(all(x["transcript_feature"] == "CDS" for x in annotations))

            with open(os.path.join(out_dir, "ase_candidate_qc.json")) as handle:
                qc = json.load(handle)
            self.assertEqual(
                qc["allele_config_qc"]["full"]["rows_malformed_skipped"], 1
            )


if __name__ == "__main__":
    unittest.main()
