#!/usr/bin/env python3

import csv
import gzip
import hashlib
import importlib.util
import json
import os
import shutil
import sys
import unittest


HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "call_phaser_ase_candidates.py")
sys.path.insert(0, HERE)
SPEC = importlib.util.spec_from_file_location("caller_v025", SCRIPT)
CALLER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CALLER)
import validate_fragment_outputs as VALIDATOR
PROJECT_TMP = os.path.abspath(
    os.environ.get("ASE_PROJECT_TMP", os.path.join(HERE, "tmp"))
)


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


class FakeRead:
    def __init__(
        self,
        name,
        start,
        length,
        bases_by_reference_position,
        reverse,
        read1=True,
    ):
        self.query_name = name
        self.reference_start = start
        self.reference_end = start + length
        self.mapping_quality = 255
        self.template_length = 160
        self.is_reverse = reverse
        self.is_read1 = read1
        self.is_read2 = not read1
        self.is_paired = True
        self.is_proper_pair = True
        self.is_unmapped = False
        self.is_secondary = False
        self.is_supplementary = False
        self.is_qcfail = False
        self.is_duplicate = False
        sequence = ["N"] * length
        for position, base in bases_by_reference_position.items():
            query_index = (position - 1) - start
            if 0 <= query_index < length:
                sequence[query_index] = base
        self.query_sequence = "".join(sequence)
        self.query_qualities = [40] * length

    def get_tag(self, tag):
        if tag == "AS":
            return 50
        raise KeyError(tag)

    def get_aligned_pairs(self, matches_only=True):
        return [
            (index, self.reference_start + index)
            for index in range(len(self.query_sequence))
        ]


class FakeAlignmentFile:
    registry = {}

    def __init__(self, path, mode, threads=0):
        self.path = path
        self.references = ("1",)
        self.reads = list(self.registry[path])

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        return False

    def fetch(self, contig, start, end):
        return iter(
            read for read in self.reads
            if read.reference_start < end and read.reference_end > start
        )


class FakePysam:
    AlignmentFile = FakeAlignmentFile


class V025BamIntegrationTest(unittest.TestCase):
    def setUp(self):
        os.makedirs(PROJECT_TMP, exist_ok=True)
        self.root = os.path.join(PROJECT_TMP, "v025_bam_integration")
        if os.path.exists(self.root):
            shutil.rmtree(self.root)
        os.makedirs(self.root)

    def tearDown(self):
        if os.path.exists(self.root):
            shutil.rmtree(self.root)

    def test_direct_and_strand_resolved_fragments_merge_into_counts(self):
        sample = "S1"
        variants = ("1_101_A_G", "1_121_C_T")
        digest = hashlib.sha1(
            ",".join(variants).encode("utf-8")
        ).hexdigest()[:10]
        component_id = "S1:1:101-121:C{}".format(digest)

        gtf = os.path.join(self.root, "genes.gtf")
        with open(gtf, "w") as handle:
            handle.write(
                '1\ttest\tgene\t90\t220\t.\t+\t.\t'
                'gene_id "GP"; gene_name "GenePlus"; '
                'gene_biotype "protein_coding";\n'
            )
            handle.write(
                '1\ttest\ttranscript\t90\t220\t.\t+\t.\t'
                'gene_id "GP"; transcript_id "TP"; '
                'gene_name "GenePlus"; gene_biotype "protein_coding"; '
                'transcript_biotype "protein_coding";\n'
            )
            handle.write(
                '1\ttest\texon\t90\t220\t.\t+\t.\t'
                'gene_id "GP"; transcript_id "TP"; '
                'gene_name "GenePlus"; gene_biotype "protein_coding"; '
                'transcript_biotype "protein_coding";\n'
            )
            handle.write(
                '1\ttest\tgene\t95\t110\t.\t-\t.\t'
                'gene_id "GM"; gene_name "GeneMinus"; '
                'gene_biotype "protein_coding";\n'
            )
            handle.write(
                '1\ttest\ttranscript\t95\t110\t.\t-\t.\t'
                'gene_id "GM"; transcript_id "TM"; '
                'gene_name "GeneMinus"; gene_biotype "protein_coding"; '
                'transcript_biotype "protein_coding";\n'
            )
            handle.write(
                '1\ttest\texon\t95\t110\t.\t-\t.\t'
                'gene_id "GM"; transcript_id "TM"; '
                'gene_name "GeneMinus"; gene_biotype "protein_coding"; '
                'transcript_biotype "protein_coding";\n'
            )

        haplotypes = os.path.join(self.root, "haplotypes.tsv")
        write_tsv(
            haplotypes,
            [
                "contig", "start", "stop", "length", "variants",
                "variant_ids", "variant_alleles", "reads_hap_a",
                "reads_hap_b", "reads_total", "edges_supporting",
                "edges_total", "annotated_phase", "phase_concordant",
                "gw_phase", "gw_confidence"
            ],
            [[
                "1", "101", "121", "20", "2", ",".join(variants),
                "A,C|G,T", "3", "1", "4", "1", "1", "00|11",
                "1", "00|11", "1"
            ]]
        )
        connections = os.path.join(self.root, "connections.tsv")
        write_tsv(
            connections,
            [
                "variant_a", "variant_b", "supporting_connections",
                "total_connections", "conflicting_configuration_p",
                "phase_concordant"
            ],
            [[variants[0], variants[1], "4", "4", "1", "1"]]
        )
        counts = os.path.join(self.root, "haplotypic_counts.tsv")
        write_tsv(
            counts,
            [
                "contig", "start", "stop", "variants", "variantCount",
                "variantsBlacklisted", "variantCountBlacklisted",
                "haplotypeA", "haplotypeB", "aCount", "bCount",
                "totalCount", "blockGWPhase", "gwStat",
                "max_haplo_maf", "bam", "aReads", "bReads"
            ],
            [[
                "1", "101", "121", ",".join(variants), "2", "", "0",
                "A,C", "G,T", "3", "1", "4", "0|1", "1", "0",
                "S1_Tissue.phaseready", "", ""
            ]]
        )
        phased_vcf = os.path.join(self.root, "wgs.phased.vcf")
        with open(phased_vcf, "w") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            handle.write(
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t"
                "FORMAT\tS1\n"
            )
            handle.write(
                "1\t101\t.\tA\tG\t60\tPASS\t.\tGT:GQ:PS\t"
                "0|1:60:101\n"
            )
            handle.write(
                "1\t121\t.\tC\tT\t60\tPASS\t.\tGT:GQ:PS\t"
                "0|1:60:101\n"
            )
        wgs_counts = os.path.join(self.root, "wgs.allelic_counts.tsv")
        write_tsv(
            wgs_counts,
            [
                "contig", "position", "variantID", "refAllele",
                "altAllele", "refCount", "altCount", "totalCount"
            ],
            [
                ["1", "101", variants[0], "A", "G", "20", "20", "40"],
                ["1", "121", variants[1], "C", "T", "18", "22", "40"],
            ],
        )

        bam_root = os.path.join(self.root, "bamroot")
        bam_dir = os.path.join(bam_root, "Tissue", sample)
        os.makedirs(bam_dir)
        bam_path = os.path.join(
            bam_dir, "S1_Tissue.phaseready.bam"
        )
        reads = [
            # reverse R1 -> '+' transcript; shared + unique evidence -> GP
            FakeRead(
                "gp_unique_A", 90, 50, {101: "A", 121: "C"},
                reverse=True
            ),
            # reverse R1 -> '+' transcript; shared exon alone -> GP rescue
            FakeRead(
                "gp_rescue_A", 90, 20, {101: "A"}, reverse=True
            ),
            # forward R1 -> '-' transcript; shared exon alone -> GM rescue
            FakeRead(
                "gm_rescue_B", 90, 20, {101: "G"}, reverse=False
            ),
            # forward R1 -> '-', but SNP 121 is uniquely GP(+) exonic.
            FakeRead(
                "gp_unique_strand_mismatch", 115, 20, {121: "C"},
                reverse=False
            ),
        ]
        with open(bam_path, "wb") as handle:
            handle.write(b"fake-bam-for-unit-test\n")
        with open(bam_path + ".bai", "wb") as handle:
            handle.write(b"fake-index-for-unit-test\n")
        FakeAlignmentFile.registry[bam_path] = reads
        CALLER.bam_engine.pysam = FakePysam

        out_dir = os.path.join(self.root, "out")
        rc = CALLER.main(
            [
                "--sample", sample,
                "--gtf", gtf,
                "--variant-connections", connections,
                "--haplotypes", haplotypes,
                "--haplotypic-counts", counts,
                "--phased-wgs-vcf", phased_vcf,
                "--phaser-wgs-allelic-counts", wgs_counts,
                "--bam-root", bam_root,
                "--write-fragments", "1",
                "--exclude-tissues", "",
                "--min-total", "1",
                "--min-each-haplotype", "0",
                "--out-dir", out_dir,
                "--force",
            ]
        )
        self.assertEqual(rc, 0)

        with open(
            os.path.join(out_dir, "audit", "ase_count_validation.tsv")
        ) as handle:
            validation = next(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(validation["component_id"], component_id)
        self.assertEqual(validation["count_match_exact"], "1")
        self.assertEqual(validation["reconstructed_phaserlike_a"], "3")
        self.assertEqual(validation["reconstructed_phaserlike_b"], "1")
        self.assertEqual(
            validation["n_annotated_exon_strand_mismatch"], "1"
        )

        with open(
            os.path.join(out_dir, "audit", "ase_strand_assignment_qc.tsv")
        ) as handle:
            gene_counts = {
                row["gene_id"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }
        self.assertEqual(gene_counts["GP"]["direct_a_count"], "1")
        self.assertEqual(gene_counts["GP"]["strand_resolved_a"], "1")
        self.assertEqual(gene_counts["GP"]["strict_a_count"], "2")
        self.assertEqual(gene_counts["GM"]["direct_total_count"], "0")
        self.assertEqual(gene_counts["GM"]["strand_resolved_b"], "1")
        self.assertEqual(gene_counts["GM"]["strict_b_count"], "1")
        for row in gene_counts.values():
            self.assertEqual(
                int(row["strict_total_count"]),
                int(row["direct_total_count"])
                + int(row["strand_resolved_total"])
            )

        with gzip.open(
            os.path.join(
                out_dir, "audit", "ase_fragment_assignment.tsv.gz"
            ),
            "rt",
        ) as handle:
            fragments = {
                row["qname"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }
        self.assertEqual(
            fragments["gp_unique_A"]["gene_assignment_method"],
            "DIRECT"
        )
        self.assertEqual(
            fragments["gp_unique_A"]["strand_resolution_required"],
            "0"
        )
        for qname in ("gp_rescue_A", "gm_rescue_B"):
            self.assertEqual(fragments[qname]["final_status"], "ASSIGNED")
            self.assertEqual(
                fragments[qname]["gene_assignment_method"],
                "STRAND_RESOLVED"
            )
            self.assertEqual(
                fragments[qname]["strand_resolution_required"],
                "1"
            )
        self.assertEqual(
            fragments["gp_unique_strand_mismatch"]["final_status"],
            "ANNOTATED_EXON_STRAND_MISMATCH"
        )
        self.assertEqual(
            fragments["gp_unique_strand_mismatch"][
                "strand_resolution_required"
            ],
            "0"
        )

        with open(
            os.path.join(out_dir, "ase_haplotype_block_results.tsv")
        ) as handle:
            candidates = {
                row["gene_id"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }
        self.assertEqual(set(candidates), {"GP", "GM"})
        self.assertEqual(candidates["GP"]["ase_a_count"], "2")
        self.assertEqual(candidates["GM"]["ase_b_count"], "1")
        self.assertEqual(candidates["GP"]["n_wgs_balanced_gene_snps"], "2")
        self.assertEqual(
            candidates["GP"]["candidate_class"], "NOT_PRIMARY"
        )

        with open(os.path.join(out_dir, "ase_v025_qc.json")) as handle:
            qc = json.load(handle)
        self.assertEqual(qc["version"], "0.2.5")
        self.assertEqual(qc["n_bam_measurements_processed"], 1)
        self.assertEqual(
            qc["n_count_validation_exact_match"], 1
        )
        self.assertEqual(
            VALIDATOR.main(
                [
                    "--out-dir", out_dir,
                    "--expect-sample", sample,
                    "--require-exact-phaser-match", "1",
                ]
            ),
            0,
        )


if __name__ == "__main__":
    unittest.main()
