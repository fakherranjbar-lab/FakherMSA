"""
msa.py
------
Tum adimlari birlestiren ust seviye arayuz.

Boru hatti (pipeline):
  diziler
    -> FFT ile mesafe matrisi      (distance.py / fft_align.py)
    -> UPGMA rehber agac           (guide_tree.py)
    -> ilerlemeli hizalama         (progressive.py)
    -> hizalanmis diziler
"""

from .distance import distance_matrix
from .guide_tree import upgma
from .progressive import progressive_align, _pair_score


def align(sequences):
    """
    Bir dizi listesini Coklu Dizi Hizalamasi (MSA) yapar.

    Parametre:
        sequences: string listesi, orn. ["MKQL...", "MKKL...", ...]
    Donus:
        hizalanmis diziler listesi (hepsi esit uzunlukta, '-' = bosluk),
        orijinal sira korunur.
    """
    if len(sequences) == 0:
        return []
    if len(sequences) == 1:
        return [sequences[0]]

    D = distance_matrix(sequences)
    merges = upgma(D)
    aligned = progressive_align(sequences, merges)
    return aligned


def align_fasta(input_path, output_path=None):
    """
    FASTA dosyasini hizalar. output_path verilirse sonucu yazar.
    Donus: (basliklar, hizalanmis_diziler)
    """
    from .io_fasta import read_fasta, write_fasta
    headers, seqs = read_fasta(input_path)
    aligned = align(seqs)
    if output_path:
        write_fasta(output_path, headers, aligned)
    return headers, aligned


def sum_of_pairs_score(aligned):
    """
    Hizalama kalitesi olcusu: Sum-of-Pairs (SP) skoru.
    Tum dizi ciftlerinin, tum sutunlardaki harf-cifti skorlarinin toplami.
    Daha yuksek = daha iyi hizalama.
    """
    if not aligned:
        return 0.0
    n = len(aligned)
    length = len(aligned[0])
    total = 0.0
    for col in range(length):
        for i in range(n):
            for j in range(i + 1, n):
                total += _pair_score(aligned[i][col], aligned[j][col])
    return float(total)
