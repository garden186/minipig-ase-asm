"""Runtime configuration for asm allele prep."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
from runtime_common import *
from runtime_common import seal as seal_checkpoint
from source_identity import code_identity

IDS = configured_samples()


def release():
    return code_identity(ROOT)


def seal(parent, out, identity, summary):
    return seal_checkpoint(parent, out, identity, summary, "PASS_ALLELE_PREPARATION_STAGE")
