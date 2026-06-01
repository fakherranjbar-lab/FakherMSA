"""
progressive.py
--------------
Ilerlemeli hizalama (progressive alignment): rehber agacin verdigi
sirayla, kumeleri (profilleri) ikiser ikiser birlestiririz.

Bir "profil" = ayni uzunlukta, hizalanmis dizilerin listesidir.
Iki profili birlestirmek icin profil-profil Needleman-Wunsch (global
dinamik programlama) kullaniyoruz. Sutun-sutun skoru, iki profilin
o sutunlarindaki tum harf ciftlerinin ortalama skorudur (sum-of-pairs).
"""

import numpy as np

GAP_PENALTY = -8        # bosluk (gap) cezasi
MATCH = 5               # ayni harf
MISMATCH = -2           # farkli harf


def _pair_score(a, b):
    """Iki harf arasi temel skor."""
    if a == '-' or b == '-':
        return GAP_PENALTY
    return MATCH if a == b else MISMATCH


def _column_score(profile_a, profile_b, i, j):
    """profile_a'nin i. sutunu ile profile_b'nin j. sutunu arasi ortalama skor."""
    col_a = [row[i] for row in profile_a]
    col_b = [row[j] for row in profile_b]
    total = 0.0
    for ca in col_a:
        for cb in col_b:
            total += _pair_score(ca, cb)
    return total / (len(col_a) * len(col_b))


def align_two_profiles(profile_a, profile_b):
    """
    Iki hizalanmis profili (dizi listesi) profil-profil NW ile birlestirir.
    Donus: birlesik profil (yeni dizi listesi, hepsi esit uzunlukta).
    """
    la, lb = len(profile_a[0]), len(profile_b[0])

    # DP tablosu ve geri-izleme tablosu
    dp = np.zeros((la + 1, lb + 1), dtype=float)
    trace = np.zeros((la + 1, lb + 1), dtype=int)  # 0=cap, 1=yukari, 2=sol

    for i in range(1, la + 1):
        dp[i, 0] = dp[i - 1, 0] + GAP_PENALTY
        trace[i, 0] = 1
    for j in range(1, lb + 1):
        dp[0, j] = dp[0, j - 1] + GAP_PENALTY
        trace[0, j] = 2

    for i in range(1, la + 1):
        for j in range(1, lb + 1):
            diag = dp[i - 1, j - 1] + _column_score(profile_a, profile_b, i - 1, j - 1)
            up = dp[i - 1, j] + GAP_PENALTY
            left = dp[i, j - 1] + GAP_PENALTY
            best = max(diag, up, left)
            dp[i, j] = best
            trace[i, j] = 0 if best == diag else (1 if best == up else 2)

    # Geri izleme ile sutun islemlerini topla
    i, j = la, lb
    ops = []
    while i > 0 or j > 0:
        t = trace[i, j]
        if i > 0 and j > 0 and t == 0:
            ops.append('diag'); i -= 1; j -= 1
        elif i > 0 and t == 1:
            ops.append('up'); i -= 1
        else:
            ops.append('left'); j -= 1
    ops.reverse()

    new_a = ['' for _ in profile_a]
    new_b = ['' for _ in profile_b]
    ai = bj = 0
    for op in ops:
        if op == 'diag':
            for k, row in enumerate(profile_a):
                new_a[k] += row[ai]
            for k, row in enumerate(profile_b):
                new_b[k] += row[bj]
            ai += 1; bj += 1
        elif op == 'up':              # profile_a ilerler, profile_b'ye bosluk
            for k, row in enumerate(profile_a):
                new_a[k] += row[ai]
            for k in range(len(profile_b)):
                new_b[k] += '-'
            ai += 1
        else:                          # profile_b ilerler, profile_a'ya bosluk
            for k in range(len(profile_a)):
                new_a[k] += '-'
            for k, row in enumerate(profile_b):
                new_b[k] += row[bj]
            bj += 1

    return new_a + new_b


def progressive_align(sequences, merges):
    """
    sequences: ham dizi listesi
    merges: guide_tree.upgma() ciktisi (birlestirme sirasi)
    Donus: (hizalanmis_diziler, eslestirme_indeksleri)
       hizalanmis_diziler, ORIJINAL dizi sirasina gore dondurulur.
    """
    # Her dizi tek elemanli bir profil olarak baslar.
    # profiles[id] = (dizi_listesi, orijinal_indeks_listesi)
    profiles = {i: ([seq], [i]) for i, seq in enumerate(sequences)}

    for i, j in merges:
        new_id = min(i, j)
        seqs_a, idx_a = profiles[i]
        seqs_b, idx_b = profiles[j]
        merged_seqs = align_two_profiles(seqs_a, seqs_b)
        profiles[new_id] = (merged_seqs, idx_a + idx_b)
        if max(i, j) in profiles and max(i, j) != new_id:
            del profiles[max(i, j)]

    # Geriye tek profil kalir
    final_seqs, final_idx = next(iter(profiles.values()))

    # Orijinal sıraya geri diz
    ordered = [None] * len(sequences)
    for pos, original_index in enumerate(final_idx):
        ordered[original_index] = final_seqs[pos]
    return ordered
