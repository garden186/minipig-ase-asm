"""Build the original Ensembl 115 canonical-transcript annotation index."""
import argparse
from pathlib import Path
from common import GTF_SHA256
from reference import build
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--gtf',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--expected-gtf-sha256',default=GTF_SHA256);p.add_argument('--memory-gib',type=float,default=8);p.add_argument('--free-disk-gib',type=float,default=20)
    a=p.parse_args();build(a.gtf,a.output,expected_hash=a.expected_gtf_sha256,max_gib=a.memory_gib,min_free_gib=a.free_disk_gib)
if __name__=='__main__':main()
