"""
io_fasta.py
-----------
Basit FASTA okuyucu/yazici. Dis kutuphane gerektirmez.
"""


def read_fasta(path):
    """
    FASTA dosyasini okur.
    Donus: (basliklar, diziler) -> iki ayri liste, ayni sirada.
    """
    headers, seqs = [], []
    current = []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current:
                    seqs.append(''.join(current))
                    current = []
                headers.append(line[1:].strip())
            else:
                current.append(line)
        if current:
            seqs.append(''.join(current))
    return headers, seqs


def write_fasta(path, headers, seqs, line_width=60):
    """Hizalanmis dizileri FASTA olarak yazar (bosluklar '-' korunur)."""
    with open(path, 'w', encoding='utf-8') as f:
        for h, s in zip(headers, seqs):
            f.write('>' + h + '\n')
            for k in range(0, len(s), line_width):
                f.write(s[k:k + line_width] + '\n')


def format_alignment(headers, seqs, name_width=12):
    """Hizalamayi okunabilir metin (Clustal benzeri) olarak bicimlendirir."""
    lines = []
    for h, s in zip(headers, seqs):
        name = (h[:name_width]).ljust(name_width)
        lines.append(f"{name} {s}")
    return '\n'.join(lines)
