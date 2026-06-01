"""
example.py
----------
Kutuphanenin tipik kullanimi. Calistirmak icin:

    cd FakherMSA
    PYTHONPATH=src python examples/example.py
"""

import fakhermsa as msa

sequences = [
    "MKQLEDKVEELLSKNYHLENEVARLKKLVGER",
    "MKQLEDKVEELLSKNYHLENEVARLKKLVGE",
    "MKKLEDKVEELLSKNAHLENEVARLKKLVGD",
    "MTQLEDRVEELLSQNYHLENQVARLKALVGN",
]

print("=== 1) Cift bazli FFT benzerligi ===")
print("seq1-seq2 benzerlik:", round(msa.fft_similarity(sequences[0], sequences[1]), 3))
print("seq1-seq4 benzerlik:", round(msa.fft_similarity(sequences[0], sequences[3]), 3))
off, peak = msa.best_offset(sequences[0], sequences[2])
print("seq1-seq3 en iyi kayma (offset):", off)

print("\n=== 2) Mesafe matrisi (FFT) ===")
import numpy as np
print(np.round(msa.distance_matrix(sequences), 3))

print("\n=== 3) Rehber agac (UPGMA birlestirme sirasi) ===")
print(msa.upgma(msa.distance_matrix(sequences)))

print("\n=== 4) Coklu Dizi Hizalamasi ===")
aligned = msa.align(sequences)
for i, row in enumerate(aligned):
    print(f"seq{i+1}: {row}")

print("\nSum-of-Pairs skoru:", msa.sum_of_pairs_score(aligned))
