"""
FakherMSA
===========
FFT tabanli (MAFFT mantigi) Coklu Dizi Hizalamasi kutuphanesi.

Hizli kullanim:
    >>> import fakhermsa as msa
    >>> aligned = msa.align(["MKQLED", "MKKLED", "MTQLED"])
    >>> for row in aligned:
    ...     print(row)

NOT: Bu kutuphaneyi kendi adiniza yayinlamak icin "YourName" / "fakhermsa"
ifadelerini kendi adinizla degistirin (README'deki adimlara bakin).
"""

from .msa import align, align_fasta, sum_of_pairs_score
from .fft_align import fft_similarity, best_offset, correlation_profile
from .distance import distance_matrix
from .guide_tree import upgma
from .progressive import align_two_profiles
from .io_fasta import read_fasta, write_fasta, format_alignment

__version__ = "0.1.0"
__all__ = [
    "align", "align_fasta", "sum_of_pairs_score",
    "fft_similarity", "best_offset", "correlation_profile",
    "distance_matrix", "upgma", "align_two_profiles",
    "read_fasta", "write_fasta", "format_alignment",
]
