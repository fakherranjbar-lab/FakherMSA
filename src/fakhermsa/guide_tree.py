"""
guide_tree.py
-------------
Rehber agac (guide tree): hangi dizilerin once birlestirilecegini
belirleyen sira. En benzer (en kucuk mesafeli) ciftler once birlesir.

UPGMA (Unweighted Pair Group Method with Arithmetic mean) kullaniyoruz:
  1) Mesafe matrisindeki en kucuk degeri bul -> o iki kumeyi birlestir.
  2) Yeni kumenin diger kumelere mesafesini, eleman sayilarina gore
     agirlikli ortalama ile guncelle.
  3) Tek kume kalana kadar tekrarla.

Donus: birlestirme adimlarinin sirasi [(i, j), ...].
Bu sira progressive.py tarafindan kullanilir.
"""

import numpy as np


def upgma(D):
    """
    D: NxN mesafe matrisi (numpy)
    Donus: birlestirme adimlari listesi. Her adim (a, b) = a ve b kumelerini birlestir.
           Birlesim sonucu kume kimligi min(a, b) olur.
    """
    n = D.shape[0]
    Dm = D.astype(float).copy()
    sizes = {i: 1 for i in range(n)}     # her kumenin eleman sayisi
    active = list(range(n))              # hala aktif kume kimlikleri
    merges = []

    while len(active) > 1:
        # 1) En kucuk mesafeli cifti bul
        best = None
        for a in range(len(active)):
            for b in range(a + 1, len(active)):
                i, j = active[a], active[b]
                if best is None or Dm[i, j] < best[0]:
                    best = (Dm[i, j], i, j)
        _, i, j = best
        merges.append((i, j))

        new_id = min(i, j)
        old_id = max(i, j)
        # 2) Mesafeleri agirlikli ortalama ile guncelle
        for k in active:
            if k != i and k != j:
                d = (Dm[i, k] * sizes[i] + Dm[j, k] * sizes[j]) / (sizes[i] + sizes[j])
                Dm[new_id, k] = Dm[k, new_id] = d
        sizes[new_id] = sizes[i] + sizes[j]
        # 3) Eski kumeyi kaldir
        active.remove(old_id)

    return merges
