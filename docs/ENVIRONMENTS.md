# Environment evidence

## Local curation tests

Python 3.12.14, NumPy 2.5.2, SciPy 1.18.1 and psutil 7.2.2 were available in the local test runtime. The numerical requirements file records those versions; it is not evidence of historical preprocessing versions. Run the scripts in separate processes to avoid collisions among preserved native module names.

Actual BAM-based ASE counting additionally requires pysam. Its original execution version has not been authenticated. Full preprocessing also requires fastp, BWA, SAMtools, BCFtools, GATK/Java, Docker/DeepVariant, WhatsHap, mosdepth, STAR, Trim Galore/Cutadapt, phASER, Bismark and Bowtie2. Executables are supplied by the user or configured explicitly. No software is installed silently.

## Historical version evidence

The archived final GO environment in `environments/go-linux-64.explicit.txt` records R 4.5.3 and clusterProfiler 4.18.4. It can be used with a compatible Linux Conda installation; the export contains public package URLs, not biological data. The public R entry point has not been rerun during Windows curation.

The manuscript declares fastp 1.3.3, BWA-MEM 0.7.18, SAMtools/BCFtools 1.23.1, GATK 4.6.2.0, DeepVariant 1.10.0, WhatsHap 2.8, RNA Trim Galore 2.2.0, STAR 2.7.10b and phASER 1.2.0. Not all declarations have been cross-checked against original run logs. The old GitHub WGS environment specified BWA 0.7.17; it is not adopted as proof of the manuscript's run environment.

Still to authenticate from original execution records: WGBS Trim Galore/Bismark/Bowtie2, mosdepth, pysam, stage-specific Python/numerical environments, container digest, reference FASTA checksum and complete tool-index identities. Current `--version` output documents a rerun environment; it does not establish which version generated past results.

The requested reference filenames are `Sus_scrofa.Sscrofa11.1.dna.toplevel.fa` and `Sus_scrofa.Sscrofa11.1.115.gtf`. The original annotation code pins GTF SHA256 `1d94edba462adb702ee9bf16a68e455684d48a5b07fc3b12871661bc40079cbe`; verify the actual supplied bytes before use.
