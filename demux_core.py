import regex
import os
from Bio import Align
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from typing import Dict, Tuple, Optional
from utils import get_reverse_complement

# --- КОНСТАНТЫ И НАСТРОЙКИ ---
UNIVERSAL_FWD_ADAPTER = "TGTTGGTGCAGATATTGCGGG"
UNIVERSAL_REV_ADAPTER = "TTGCCTGTCGCTCTATCTTCAAA"

# 1. Возвращаем МЯГКИЕ лимиты на якоря (до 10 ошибок суммарно)
# Так как мы ищем оба якоря по краям, риск химер минимален, но мы спасаем грязные риды
ANCHOR_ERRORS = "{e<=10}"
SEARCH_WINDOW = 120  # Расширили окно поиска, на случай если фланки длинные

# 2. Настраиваем ту самую умную систему штрафов для баркодов
aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 2  # Награда за совпадение
aligner.mismatch_score = -1  # Штраф за замену
aligner.open_gap_score = -10  # Мощный штраф за индел (гэп)
aligner.extend_gap_score = -10

# Прекомпиляция адаптеров
FWD_RC = get_reverse_complement(UNIVERSAL_FWD_ADAPTER)
REV_RC = get_reverse_complement(UNIVERSAL_REV_ADAPTER)

REG_FWD = regex.compile(f"(?b)({UNIVERSAL_FWD_ADAPTER}){ANCHOR_ERRORS}")
REG_REV_RC = regex.compile(f"(?b)({REV_RC}){ANCHOR_ERRORS}")
REG_REV = regex.compile(f"(?b)({UNIVERSAL_REV_ADAPTER}){ANCHOR_ERRORS}")
REG_FWD_RC = regex.compile(f"(?b)({FWD_RC}){ANCHOR_ERRORS}")


def find_anchors(sequence: str) -> Tuple[Optional[str], Optional[regex.Match], Optional[regex.Match], int]:
    """Быстрый поиск грязных адаптеров по краям."""
    seq_len = len(sequence)
    if seq_len < SEARCH_WINDOW * 2:
        return None, None, None, 0

    end_zone = seq_len - SEARCH_WINDOW

    # PLUS-ориентация
    m_fwd = REG_FWD.search(sequence, 0, SEARCH_WINDOW)
    if m_fwd:
        m_rev_rc = REG_REV_RC.search(sequence, end_zone)
        if m_rev_rc:
            errs = sum(m_fwd.fuzzy_counts) + sum(m_rev_rc.fuzzy_counts)
            return "PLUS", m_fwd, m_rev_rc, errs

    # MINUS-ориентация
    m_rev = REG_REV.search(sequence, 0, SEARCH_WINDOW)
    if m_rev:
        m_fwd_rc = REG_FWD_RC.search(sequence, end_zone)
        if m_fwd_rc:
            errs = sum(m_rev.fuzzy_counts) + sum(m_fwd_rc.fuzzy_counts)
            return "MINUS", m_rev, m_fwd_rc, errs

    return None, None, None, 0


def get_best_barcode(flank: str, barcodes: Dict[str, str]) -> Tuple[Optional[str], int]:
    """Умный выбор баркода с учетом системы штрафов и динамического порога."""
    if not flank:
        return None, 0

    best_id = None
    best_score = -float('inf')
    best_errs_approx = 0

    for b_id, b_seq in barcodes.items():
        score = aligner.score(flank, b_seq)
        max_possible_score = len(b_seq) * aligner.match_score

        # Очень мягкий порог прохождения (достаточно набрать 30% от идеала),
        # но за счет best_score мы всегда выберем сильнейшего кандидата
        min_passing_score = max_possible_score * 0.3

        if score >= min_passing_score and score > best_score:
            best_score = score
            best_id = b_id
            # Имитируем кол-во ошибок для красивых графиков в seaborn
            best_errs_approx = max(0, int((max_possible_score - score) / 3))

    return best_id, best_errs_approx


def process_read(sequence: str, fbs: dict, rbs: dict, fbs_rc: dict, rbs_rc: dict) -> Optional[dict]:
    orientation, m1, m2, anchor_errs = find_anchors(sequence)
    if not orientation:
        return None

    if orientation == "PLUS":
        fb_flank = sequence[:m1.start()]
        rb_flank = sequence[m2.end():]

        found_fb, fb_err = get_best_barcode(fb_flank, fbs)
        found_rb, rb_err = get_best_barcode(rb_flank, rbs_rc)

    else:  # MINUS
        rb_flank = sequence[:m1.start()]
        fb_flank = sequence[m2.end():]

        found_rb, rb_err = get_best_barcode(rb_flank, rbs)
        found_fb, fb_err = get_best_barcode(fb_flank, fbs_rc)

    if found_fb and found_rb:
        return {
            "sample_id": f"{found_fb}_{found_rb}",
            "orientation": orientation,
            "trim_start": m1.end(),
            "trim_end": m2.start(),
            "total_errors": anchor_errs + fb_err + rb_err
        }
    return None


def run_demux(fastq_path: str, fbs: Dict[str, str], rbs: Dict[str, str], output_dir: str) -> dict:
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fbs_rc = {k: get_reverse_complement(v) for k, v in fbs.items()}
    rbs_rc = {k: get_reverse_complement(v) for k, v in rbs.items()}

    # Структура словаря четко под твою visualization.py
    stats = {
        "total": 0, "demuxed": 0, "unassigned": 0,
        "sample_counts": {},
        "error_distribution": {i: 0 for i in range(35)}  # расширил до 35 ошибок на всякий случай
    }
    handles = {}

    print("Анализ ридов запущен...")

    with open(fastq_path, "r") as handle:
        for title, seq, qual in FastqGeneralIterator(handle):
            stats["total"] += 1
            if stats["total"] % 5000 == 0:
                print(f"Обработано {stats['total']} ридов...", end="\r")

            res = process_read(seq, fbs, rbs, fbs_rc, rbs_rc)

            if res:
                s_id = res["sample_id"]
                stats["demuxed"] += 1
                stats["sample_counts"][s_id] = stats["sample_counts"].get(s_id, 0) + 1

                err_idx = min(res["total_errors"], 34)
                stats["error_distribution"][err_idx] += 1

                if s_id not in handles:
                    handles[s_id] = open(os.path.join(output_dir, f"{s_id}.fastq"), "w")

                t_start, t_end = res["trim_start"], res["trim_end"]
                trimmed_seq = seq[t_start:t_end]
                trimmed_qual = qual[t_start:t_end]

                if res["orientation"] == "MINUS":
                    trimmed_seq = get_reverse_complement(trimmed_seq)
                    trimmed_qual = trimmed_qual[::-1]

                handles[s_id].write(f"@{title}\n{trimmed_seq}\n+\n{trimmed_qual}\n")
            else:
                stats["unassigned"] += 1

    for h in handles.values():
        h.close()

    return stats