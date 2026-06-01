# FakherMSA

FFT tabanlı (MAFFT mantığı) **Çoklu Dizi Hizalaması (Multiple Sequence Alignment, MSA)** kütüphanesi.
Bioinformatik final projesi için **sıfırdan, kendi elimle** yazılmıştır. Tek bağımlılığı `numpy`'dir.

> **Atanan algoritma:** Öğrenci numarası `221201910` → `221201910 % 4 = 2` → **MAFFT**

---

## ⚠️ Önce: Kütüphaneyi kendi adınıza çevirin

Hoca "kendi adınıza bir kütüphane" istiyor (`pip install YosefMSA` örneği gibi). 3 yeri kendi adınızla değiştirin (örn. adınız *Ahmet* ise `AhmetMSA`):

1. **Klasör adı:** `FakherMSA/` → `AhmetMSA/`
2. **İç paket klasörü:** `src/fakhermsa/` → `src/ahmetmsa/` (küçük harf, boşluksuz)
3. **`pyproject.toml`** içinde:
   - `name = "FakherMSA"` → `"AhmetMSA"`
   - `fakhermsa = "fakhermsa.cli:main"` → `ahmetmsa = "ahmetmsa.cli:main"`
   - `authors`, `Homepage`, `Repository` alanlarını kendinize göre doldurun
4. Tüm `import fakhermsa` satırlarını `import ahmetmsa` yapın.

Hızlı toplu değiştirme (Linux/Mac):
```bash
grep -rl "fakhermsa" . | xargs sed -i 's/fakhermsa/ahmetmsa/g'
grep -rl "FakherMSA" . | xargs sed -i 's/FakherMSA/AhmetMSA/g'
mv src/fakhermsa src/ahmetmsa
```

---

## Kurulum

### Yerel (geliştirme) kurulumu
```bash
cd FakherMSA
pip install -e .
```

### GitHub'dan kurulum
```bash
pip install git+https://github.com/KULLANICI_ADINIZ/FakherMSA.git
```

### PyPI'dan kurulum (yayınladıktan sonra)
```bash
pip install FakherMSA
```

---

## Kullanım

### Python içinde
```python
import fakhermsa as msa

seqs = [
    "MKQLEDKVEELLSKNYHLENEVARLKKLVGER",
    "MKKLEDKVEELLSKNAHLENEVARLKKLVGD",
    "MTQLEDRVEELLSQNYHLENQVARLKALVGN",
]

aligned = msa.align(seqs)
for row in aligned:
    print(row)

print("SP skoru:", msa.sum_of_pairs_score(aligned))
```

### Komut satırından
```bash
fakhermsa girdi.fasta -o cikti.fasta --score
# veya kurulmadan:
PYTHONPATH=src python -m fakhermsa examples/sample.fasta --score
```

---

## Nasıl çalışır? (Boru hattı)

```
diziler
  │
  ├─ 1) Fizikokimyasal kodlama  (hacim + polarite)        → encoding.py
  ├─ 2) FFT çapraz-korelasyon  (hızlı benzerlik/offset)   → fft_align.py
  ├─ 3) Mesafe matrisi          (1 - benzerlik)            → distance.py
  ├─ 4) UPGMA rehber ağaç       (birleştirme sırası)       → guide_tree.py
  └─ 5) İlerlemeli hizalama     (profil-profil NW)         → progressive.py
  │
  ▼
hizalanmış diziler
```

**MAFFT'in imza adımı**, dizileri amino asitlerin *hacim* ve *polarite* değerlerinden
oluşan sayısal sinyallere çevirip, iki diziyi karşılaştırmayı **FFT** ile yapmaktır.
Çapraz-korelasyon teoremi sayesinde bu işlem O(n²) yerine **O(n log n)** olur — MAFFT'i
"Fast" yapan budur.

---

## Test

```bash
PYTHONPATH=src python tests/test_basic.py
# pytest varsa:
PYTHONPATH=src python -m pytest tests/ -v
```

---

## PyPI'a yükleme (özet)

```bash
pip install build twine
python -m build                 # dist/ altında .whl ve .tar.gz üretir
python -m twine upload dist/*   # PyPI hesabı/token gerekir
```

## Lisans
MIT (bkz. `LICENSE`).
