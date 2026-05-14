from Bio.Seq import Seq
from Bio import SeqIO
from typing import Dict, Tuple

def get_reverse_complement(seq: str) -> str:
    """Возвращает обратно-комплементарную последовательность."""
    return str(Seq(seq).reverse_complement())

def parse_barcodes(fasta_path: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Парсит баркоды из FASTA файла."""
    fbs, rbs = {}, {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        # Очищаем последовательности от пробелов на всякий случай
        seq = str(record.seq).strip()
        if record.id.startswith("FB"):
            fbs[record.id] = seq
        elif record.id.startswith("RB"):
            rbs[record.id] = seq
    return fbs, rbs