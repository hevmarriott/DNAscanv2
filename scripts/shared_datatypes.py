from typing import NewType, Dict
import pysam.libcfaidx

Variant = NewType('Variant', dict)
RefGenome = NewType('RefGenome', pysam.libcfaidx.FastaFile)
RefSeq = NewType('RefSeq', str)
StartPos = NewType('StartPos', Dict[str, int])
