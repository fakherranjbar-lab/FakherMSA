"""
encoding.py
-----------
MAFFT'in temel fikri: amino asitleri DOGRUDAN harf olarak degil,
fizikokimyasal ozelliklerine gore SAYISAL sinyallere cevirmek.
Boylece iki diziyi karsilastirmak = iki sinyali karsilastirmak olur
ve bu islem FFT (Hizli Fourier Donusumu) ile cok hizli yapilabilir.

Her amino asit icin iki ozellik kullaniyoruz (orijinal MAFFT, Katoh 2002):
  1) Hacim (volume)   -> kalintinin kapladigi yer
  2) Polarite          -> su sevme/sevmeme egilimi

Degerleri normalize ediyoruz (ortalama 0, standart sapma 1) ki
FFT korelasyonu sapmasiz olsun.
"""

import numpy as np

# Yaklasik kalinti hacimleri (Angstrom^3) - literaturden
VOLUME = {
    'A': 88.6, 'R': 173.4, 'N': 114.1, 'D': 111.1, 'C': 108.5,
    'Q': 143.8, 'E': 138.4, 'G': 60.1, 'H': 153.2, 'I': 166.7,
    'L': 166.7, 'K': 168.6, 'M': 162.9, 'F': 189.9, 'P': 112.7,
    'S': 89.0, 'T': 116.1, 'W': 227.8, 'Y': 193.6, 'V': 140.0,
}

# Grantham (1974) polarite degerleri
POLARITY = {
    'A': 8.1, 'R': 10.5, 'N': 11.6, 'D': 13.0, 'C': 5.5,
    'Q': 10.5, 'E': 12.3, 'G': 9.0, 'H': 10.4, 'I': 5.2,
    'L': 4.9, 'K': 11.3, 'M': 5.7, 'F': 5.2, 'P': 8.0,
    'S': 9.2, 'T': 8.6, 'W': 5.4, 'Y': 6.2, 'V': 5.9,
}

# DNA/RNA dizileri icin basit yedek kodlama (4 nukleotidi sayisal eksene yayar)
NUCLEOTIDE_VOL = {'A': 1.0, 'C': -1.0, 'G': 0.5, 'T': -0.5, 'U': -0.5}
NUCLEOTIDE_POL = {'A': 0.5, 'C': -0.5, 'G': -1.0, 'T': 1.0, 'U': 1.0}


def _normalize(table):
    """Bir ozellik tablosunu ortalama 0, std 1 olacak sekilde normalize eder."""
    vals = np.array(list(table.values()), dtype=float)
    mean, std = vals.mean(), vals.std()
    if std == 0:
        std = 1.0
    return {k: (v - mean) / std for k, v in table.items()}


_NV = _normalize(VOLUME)
_NP = _normalize(POLARITY)
_NV_NUC = _normalize(NUCLEOTIDE_VOL)
_NP_NUC = _normalize(NUCLEOTIDE_POL)


def is_protein(seq):
    """Dizinin protein mi (amino asit) yoksa nukleotid mi oldugunu tahmin eder."""
    s = set(seq.upper()) - set('-')
    nuc = set('ACGTUN')
    # Karakterlerin %90'indan fazlasi ACGTU ise nukleotid kabul et
    if not s:
        return True
    return len(s & nuc) / len(s) < 0.9


def encode(seq):
    """
    Bir diziyi iki sayisal sinyale cevirir: (hacim_sinyali, polarite_sinyali).
    Bilinmeyen karakterler (X, bosluk vb.) 0.0 olur (notr).
    """
    seq = seq.upper()
    if is_protein(seq):
        vol_t, pol_t = _NV, _NP
    else:
        vol_t, pol_t = _NV_NUC, _NP_NUC
    v = np.array([vol_t.get(c, 0.0) for c in seq], dtype=float)
    p = np.array([pol_t.get(c, 0.0) for c in seq], dtype=float)
    return v, p
