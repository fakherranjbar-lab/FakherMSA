"""
distance.py
-----------
FFT benzerlik skorlarini kullanarak tum dizi ciftleri arasinda
bir MESAFE MATRISI olusturur. Bu matris, rehber agacin (guide tree)
girdisidir.

mesafe = 1 - benzerlik
"""

import numpy as np
from .fft_align import fft_similarity


def distance_matrix(sequences):
    """
    sequences: dizi (string) listesi
    Donus: NxN simetrik numpy mesafe matrisi (kosegen = 0)
    """
    n = len(sequences)
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            sim = fft_similarity(sequences[i], sequences[j])
            d = 1.0 - sim
            D[i, j] = D[j, i] = d
    return D
