"""
fft_align.py
------------
MAFFT'i diger MSA araclarindan AYIRAN imza adim budur:
iki diziyi karsilastirmak icin FFT (Hizli Fourier Donusumu) kullanmak.

Fikir (capraz-korelasyon teoremi):
  Iki sinyalin capraz-korelasyonu, birinin FFT'si ile digerinin
  FFT'sinin esleniginin (conjugate) carpiminin ters FFT'sine esittir.

  corr(k) = IFFT( FFT(a) * conj(FFT(b)) )

Burada k = "kayma" (offset/lag). corr(k)'nin tepe (peak) yaptigi k,
iki dizinin EN COK ortustugu (homolog oldugu) kaymayi verir.

Naif yontemle korelasyon O(n^2) iken, FFT ile O(n log n) olur.
Iste MAFFT'in "Fast" (hizli) olmasinin sebebi budur.
"""

import numpy as np
from .encoding import encode


def fft_cross_correlation(a, b):
    """
    Iki sayisal sinyalin tam capraz-korelasyonunu FFT ile hesaplar.
    Donen dizi uzunlugu len(a)+len(b)-1'dir; her indeks bir kaymaya karsilik gelir.
    """
    n = len(a) + len(b) - 1
    if n <= 0:
        return np.array([0.0])
    # Hizli FFT icin uzunlugu 2'nin kuvvetine yuvarla (sifir dolgu / zero-padding)
    size = 1 << ((n - 1).bit_length())
    fa = np.fft.rfft(a, size)
    fb = np.fft.rfft(b, size)
    corr = np.fft.irfft(fa * np.conj(fb), size)[:n]
    return corr


def correlation_profile(seq1, seq2):
    """
    Iki dizinin hacim + polarite korelasyon profilini birlestirir.
    Donus: her olasi kayma icin toplam korelasyon degeri (numpy dizisi).
    """
    v1, p1 = encode(seq1)
    v2, p2 = encode(seq2)
    return fft_cross_correlation(v1, v2) + fft_cross_correlation(p1, p2)


def best_offset(seq1, seq2):
    """
    Iki dizi arasinda korelasyonu en cok yapan kaymayi (homolog offset) dondurur.

    FFT capraz-korelasyonu DONGUSEL'dir: sonuc[k] = sum_m seq1[m+k] * seq2[m].
    Yani k = seq1'in seq2'ye gore SOLA kaymasi. NumPy'da negatif kaymalar
    dizinin sonuna sarar (wrap-around); bunlari k - N olarak yorumlariz.

    Ham korelasyon cok ortusen kaymalari kayirdigi icin her kaymayi
    ORTUSME UZUNLUGUNA boleriz (kenar yapay tepelerini engeller).

    Donus: (offset, tepe_korelasyon_degeri)
      offset > 0  -> seq1, seq2'ye gore saga kaymistir.
    """
    l1, l2 = len(seq1), len(seq2)
    v1, p1 = encode(seq1)
    v2, p2 = encode(seq2)
    n = l1 + l2 - 1
    size = 1 << ((n - 1).bit_length())
    corr = (np.fft.irfft(np.fft.rfft(v1, size) * np.conj(np.fft.rfft(v2, size)), size)
            + np.fft.irfft(np.fft.rfft(p1, size) * np.conj(np.fft.rfft(p2, size)), size))

    # Her dongusel indeksi gercek kaymaya (lag) cevir: > size/2 ise negatif
    lags = np.where(np.arange(size) <= size // 2,
                    np.arange(size), np.arange(size) - size)
    # Sadece anlamli kaymalari degerlendir (|lag| < max uzunluk)
    valid = np.abs(lags) < max(l1, l2)
    # lag = L icin ortusme sayisi: min(l2, l1 - L) - max(0, -L)
    overlap = np.minimum(l2, l1 - lags) - np.maximum(0, -lags)
    overlap = np.where(overlap > 0, overlap, 1)
    normalized = np.where(valid, corr / overlap, -np.inf)
    best_k = int(np.argmax(normalized))
    offset = int(lags[best_k])
    return offset, float(corr[best_k])


def fft_similarity(seq1, seq2):
    """
    FFT korelasyon tepesini, dizilerin "enerjisine" gore normalize ederek
    0 ile ~1 arasinda bir benzerlik skoru uretir.
    (Kosinus benzerligi mantigi: tepe / sqrt(enerji1 * enerji2))
    """
    v1, p1 = encode(seq1)
    v2, p2 = encode(seq2)
    corr = fft_cross_correlation(v1, v2) + fft_cross_correlation(p1, p2)
    peak = corr.max()
    energy1 = float(v1 @ v1 + p1 @ p1)
    energy2 = float(v2 @ v2 + p2 @ p2)
    norm = np.sqrt(energy1 * energy2) + 1e-9
    sim = peak / norm
    # Guvenlik icin [0, 1] araligina kelepcele
    return float(max(0.0, min(1.0, sim)))
