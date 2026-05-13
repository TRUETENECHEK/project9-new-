from Bio.Seq import Seq
from Bio import SeqIO
from typing import Dict, Tuple

def get_reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())

def parse_barcodes(fasta_path: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    fbs, rbs = {}, {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        if record.id.startswith("FB"):
            fbs[record.id] = str(record.seq)
        elif record.id.startswith("RB"):
            rbs[record.id] = str(record.seq)
    return fbs, rbs