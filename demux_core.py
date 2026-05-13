import regex
import os
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator  # Секретное оружие для турбо-скорости
from typing import Dict, Tuple, Optional

# --- КОНСТАНТЫ ---
UNIVERSAL_FWD_ADAPTER = "TGTTGGTGCAGATATTGCGGG"
UNIVERSAL_REV_ADAPTER = "TTGCCTGTCGCTCTATCTTCAAA"
ANCHOR_ERRORS = "{s<=6,i<=2,d<=2}"
BARCODE_ERRORS = "{e<=1}"

# Окно поиска: ищем якоря только в крайних 80 нуклеотидах рида
SEARCH_WINDOW = 80


def get_reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


# --- 1. ПРЕКОМПИЛЯЦИЯ ЯКОРЕЙ (ГЛАВНЫЙ БУСТ СКОРОСТИ) ---
FWD_RC = get_reverse_complement(UNIVERSAL_FWD_ADAPTER)
REV_RC = get_reverse_complement(UNIVERSAL_REV_ADAPTER)

# Компилируем регулярки один раз в памяти
REG_FWD = regex.compile(f"(?b)({UNIVERSAL_FWD_ADAPTER}){ANCHOR_ERRORS}")
REG_REV_RC = regex.compile(f"(?b)({REV_RC}){ANCHOR_ERRORS}")
REG_REV = regex.compile(f"(?b)({UNIVERSAL_REV_ADAPTER}){ANCHOR_ERRORS}")
REG_FWD_RC = regex.compile(f"(?b)({FWD_RC}){ANCHOR_ERRORS}")


def find_complex_adapters_fast(sequence: str) -> Tuple[
    Optional[str], Optional[regex.Match], Optional[regex.Match], int]:
    """Ищет адаптеры только по краям рида."""
    seq_len = len(sequence)
    if seq_len < SEARCH_WINDOW * 2:
        return None, None, None, 0

    end_zone_start = seq_len - SEARCH_WINDOW

    # 1. PLUS-ориентация
    # Ищем только от 0 до SEARCH_WINDOW (начало рида)
    m_fwd = REG_FWD.search(sequence, 0, SEARCH_WINDOW)
    if m_fwd:
        # Ищем только от end_zone_start до конца рида
        m_rev_rc = REG_REV_RC.search(sequence, end_zone_start)
        if m_rev_rc:
            return "PLUS", m_fwd, m_rev_rc, sum(m_fwd.fuzzy_counts) + sum(m_rev_rc.fuzzy_counts)

    # 2. MINUS-ориентация
    m_rev = REG_REV.search(sequence, 0, SEARCH_WINDOW)
    if m_rev:
        m_fwd_rc = REG_FWD_RC.search(sequence, end_zone_start)
        if m_fwd_rc:
            return "MINUS", m_rev, m_fwd_rc, sum(m_rev.fuzzy_counts) + sum(m_fwd_rc.fuzzy_counts)

    return None, None, None, 0


def find_adapters_and_match_barcodes(sequence: str, fbs_reg: dict, rbs_reg: dict) -> Optional[dict]:
    orientation, m1, m2, anchor_errors = find_complex_adapters_fast(sequence)
    if not orientation:
        return None

    if orientation == "PLUS":
        fb_flank = sequence[max(0, m1.start() - 40): m1.start() + 5]
        rb_flank = sequence[m2.end() - 5: min(len(sequence), m2.end() + 40)]

        # Проверяем скомпилированные регулярки баркодов
        found_fb = next((id for id, reg in fbs_reg['fwd'].items() if reg.search(fb_flank)), None)
        found_rb = next((id for id, reg in rbs_reg['rc'].items() if reg.search(rb_flank)), None)
    else:
        rb_flank = sequence[max(0, m1.start() - 40): m1.start() + 5]
        fb_flank = sequence[m2.end() - 5: min(len(sequence), m2.end() + 40)]

        found_rb = next((id for id, reg in rbs_reg['fwd'].items() if reg.search(rb_flank)), None)
        found_fb = next((id for id, reg in fbs_reg['rc'].items() if reg.search(fb_flank)), None)

    if found_fb and found_rb:
        return {
            "sample_id": f"{found_fb}_{found_rb}",
            "orientation": orientation,
            "start": m1.end(),
            "end": m2.start(),
            "total_errors": anchor_errors
        }
    return None


def process_and_collect_stats(fastq_path: str, fbs: Dict[str, str], rbs: Dict[str, str], output_dir: str) -> dict:
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Очищаем баркоды и сразу КОМПИЛИРУЕМ их, чтобы не делать это в цикле
    print("Подготовка и компиляция баркодов...")
    fbs_reg = {'fwd': {}, 'rc': {}}
    rbs_reg = {'fwd': {}, 'rc': {}}

    for id, seq in fbs.items():
        clean_seq = seq.replace(UNIVERSAL_FWD_ADAPTER, "")
        fbs_reg['fwd'][id] = regex.compile(f"(?b)({clean_seq}){BARCODE_ERRORS}")
        fbs_reg['rc'][id] = regex.compile(f"(?b)({get_reverse_complement(clean_seq)}){BARCODE_ERRORS}")

    for id, seq in rbs.items():
        clean_seq = seq.replace(UNIVERSAL_REV_ADAPTER, "")
        rbs_reg['fwd'][id] = regex.compile(f"(?b)({clean_seq}){BARCODE_ERRORS}")
        rbs_reg['rc'][id] = regex.compile(f"(?b)({get_reverse_complement(clean_seq)}){BARCODE_ERRORS}")

    stats = {
        "total_reads": 0, "demuxed": 0, "unassigned": 0,
        "sample_counts": {}, "error_distribution": {i: 0 for i in range(16)}
    }
    handles = {}

    print("Обработка ридов (FastqGeneralIterator)...")

    # 3. FAST I/O: Читаем файл как простые строки, это в 5-10 раз быстрее
    with open(fastq_path, "r") as handle:
        for title, seq, qual in FastqGeneralIterator(handle):
            stats["total_reads"] += 1

            if stats["total_reads"] % 5000 == 0:
                print(f"Обработано {stats['total_reads']} ридов...", end="\r")

            res = find_adapters_and_match_barcodes(seq, fbs_reg, rbs_reg)

            if res:
                sample_id = res["sample_id"]
                stats["demuxed"] += 1
                stats["sample_counts"][sample_id] = stats["sample_counts"].get(sample_id, 0) + 1

                err_idx = min(res["total_errors"], 15)
                stats["error_distribution"][err_idx] += 1

                if sample_id not in handles:
                    handles[sample_id] = open(os.path.join(output_dir, f"{sample_id}.fastq"), "w")

                # Тримминг строк
                start, end = res["start"], res["end"]
                trimmed_seq = seq[start:end]
                trimmed_qual = qual[start:end]

                # Если минус цепь, нужно перевернуть и секвенс, и качество
                if res["orientation"] == "MINUS":
                    trimmed_seq = get_reverse_complement(trimmed_seq)
                    trimmed_qual = trimmed_qual[::-1]

                # Быстрая запись сырым текстом
                handles[sample_id].write(f"@{title}\n{trimmed_seq}\n+\n{trimmed_qual}\n")
            else:
                stats["unassigned"] += 1

    for h in handles.values():
        h.close()

    print(f"\nОбработка завершена. Всего: {stats['total_reads']}")
    return stats