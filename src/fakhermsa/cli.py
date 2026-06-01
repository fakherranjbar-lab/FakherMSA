"""
cli.py
------
Komut satiri arayuzu. Kurulumdan sonra soyle calisir:

    fakhermsa girdi.fasta -o cikti.fasta

veya:

    python -m fakhermsa girdi.fasta
"""

import argparse
import sys

from .msa import align_fasta, sum_of_pairs_score
from .io_fasta import format_alignment


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="fakhermsa",
        description="FFT tabanli (MAFFT mantigi) Coklu Dizi Hizalamasi (MSA).",
    )
    parser.add_argument("input", help="Girdi FASTA dosyasi")
    parser.add_argument("-o", "--output", help="Cikti FASTA dosyasi (istege bagli)")
    parser.add_argument("--score", action="store_true",
                        help="Sum-of-Pairs skorunu da yazdir")
    args = parser.parse_args(argv)

    try:
        headers, aligned = align_fasta(args.input, args.output)
    except FileNotFoundError:
        print(f"Hata: dosya bulunamadi -> {args.input}", file=sys.stderr)
        return 1

    print(format_alignment(headers, aligned))
    if args.score:
        print("\nSum-of-Pairs skoru:", sum_of_pairs_score(aligned))
    if args.output:
        print(f"\nSonuc yazildi -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
