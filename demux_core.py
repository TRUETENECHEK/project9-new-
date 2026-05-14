import os
import sys
from Bio import Align
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from typing import Dict, Tuple, Optional
from utils import get_reverse_complement

# --- КОНСТАНТЫ И НАСТРОЙКИ ---
UNIVERSAL_FWD_ADAPTER = "TGTTGGTGCAGATATTGCGGG"
UNIVERSAL_REV_ADAPTER = "TTGCCTGTCGCTCTATCTTCAAA"

SEARCH_WINDOW = 120

FWD_RC = get_reverse_complement(UNIVERSAL_FWD_ADAPTER)
REV_RC = get_reverse_complement(UNIVERSAL_REV_ADAPTER)

# Настройка алайнера для адаптеров и баркодов
# Более мягкие штрафы, чтобы справляться с "грязными" данными и находить ЛУЧШЕЕ выравнивание
aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -3
aligner.extend_gap_score = -1


def find_best_alignment(sequence: str, query: str, min_score_ratio: float = 0.5) -> Optional[Tuple[int, int, int]]:
    """Находит лучшее локальное выравнивание и возвращает координаты и примерное количество ошибок."""
    alns = aligner.align(sequence, query)
    if not alns:
        return None
    best = alns[0]
    max_score = len(query) * aligner.match_score
    if best.score >= max_score * min_score_ratio:
        # Получаем координаты выравнивания в целевой последовательности
        target_coords = best.aligned[0]
        start = int(target_coords[0][0])
        end = int(target_coords[-1][1])
        # Примерный расчет ошибок: каждая ошибка уменьшает макс. скор примерно на 3
        errs = max(0, round((max_score - best.score) / 3.0))
        return start, end, errs
    return None


def find_anchors(sequence: str) -> Tuple[Optional[str], Optional[Tuple[int, int]], Optional[Tuple[int, int]], int]:
    """Ищет универсальные адаптеры с обоих концов рида, учитывая инверсию цепи."""
    seq_len = len(sequence)
    if seq_len < SEARCH_WINDOW * 2:
        return None, None, None, 0

    end_zone_start = seq_len - SEARCH_WINDOW
    start_seq = sequence[:SEARCH_WINDOW]
    end_seq = sequence[-SEARCH_WINDOW:]

    # Ищем адаптеры. Требуем хотя бы 55% сходства
    min_anchor_ratio = 0.55

    # Проверка PLUS цепи: FWD в начале, REV_RC в конце
    res_fwd = find_best_alignment(start_seq, UNIVERSAL_FWD_ADAPTER, min_anchor_ratio)
    if res_fwd:
        res_rev_rc = find_best_alignment(end_seq, REV_RC, min_anchor_ratio)
        if res_rev_rc:
            start_fwd, end_fwd, err_fwd = res_fwd
            start_rev_rc, end_rev_rc, err_rev_rc = res_rev_rc
            
            m1_coords = (start_fwd, end_fwd)
            m2_coords = (end_zone_start + start_rev_rc, end_zone_start + end_rev_rc)
            return "PLUS", m1_coords, m2_coords, err_fwd + err_rev_rc

    # Проверка MINUS цепи: REV в начале, FWD_RC в конце
    res_rev = find_best_alignment(start_seq, UNIVERSAL_REV_ADAPTER, min_anchor_ratio)
    if res_rev:
        res_fwd_rc = find_best_alignment(end_seq, FWD_RC, min_anchor_ratio)
        if res_fwd_rc:
            start_rev, end_rev, err_rev = res_rev
            start_fwd_rc, end_fwd_rc, err_fwd_rc = res_fwd_rc
            
            m1_coords = (start_rev, end_rev)
            m2_coords = (end_zone_start + start_fwd_rc, end_zone_start + end_fwd_rc)
            return "MINUS", m1_coords, m2_coords, err_rev + err_fwd_rc

    return None, None, None, 0


def get_best_barcode(flank: str, barcodes: Dict[str, str]) -> Tuple[Optional[str], int]:
    """Умный отбор баркода: перебирает все и выбирает лучший по штрафам."""
    if not flank:
        return None, 0

    best_id = None
    best_score_ratio = -float('inf')
    best_errs = 0

    # Проходим по всем баркодам и выбираем ТОТ, У КОТОРОГО ЛУЧШИЙ СКОР
    for b_id, b_seq in barcodes.items():
        max_score = len(b_seq) * aligner.match_score
        
        # Турбо-буст: если идеальное совпадение, проверяем, не перебьет ли это предыдущий лучший
        if b_seq in flank:
            if 1.0 > best_score_ratio:
                best_score_ratio = 1.0
                best_id = b_id
                best_errs = 0
            continue

        alns = aligner.align(flank, b_seq)
        if not alns:
            continue
        
        score = alns[0].score
        ratio = score / max_score
        
        # Минимальный порог 40% совпадения, но ищем максимум!
        if ratio >= 0.4 and ratio > best_score_ratio:
            best_score_ratio = ratio
            best_id = b_id
            best_errs = max(0, round((max_score - score) / 3.0))

    return best_id, best_errs


def process_read(sequence: str, fbs: dict, rbs: dict, fbs_rc: dict, rbs_rc: dict) -> Optional[dict]:
    """Полный цикл обработки одного рида."""
    orientation, m1, m2, anchor_errs = find_anchors(sequence)
    if not orientation:
        return None

    m1_start, m1_end = m1
    m2_start, m2_end = m2

    if orientation == "PLUS":
        fb_flank = sequence[:m1_start]
        rb_flank = sequence[m2_end:]

        found_fb, fb_err = get_best_barcode(fb_flank, fbs)
        found_rb, rb_err = get_best_barcode(rb_flank, rbs_rc)

    else:
        rb_flank = sequence[:m1_start]
        fb_flank = sequence[m2_end:]

        found_rb, rb_err = get_best_barcode(rb_flank, rbs)
        found_fb, fb_err = get_best_barcode(fb_flank, fbs_rc)

    if found_fb and found_rb:
        return {
            "sample_id": f"{found_fb}_{found_rb}",
            "orientation": orientation,
            "trim_start": m1_end,
            "trim_end": m2_start,
            "total_errors": anchor_errs + fb_err + rb_err
        }
    return None


def run_demux(fastq_path: str, fbs: Dict[str, str], rbs: Dict[str, str], output_dir: str) -> dict:
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Очищаем баркоды от универсальных адаптеров, так как они ищутся отдельно в anchors
    fbs_clean = {}
    for k, v in fbs.items():
        if v.endswith(UNIVERSAL_FWD_ADAPTER):
            fbs_clean[k] = v[:-len(UNIVERSAL_FWD_ADAPTER)]
        else:
            fbs_clean[k] = v

    rbs_clean = {}
    for k, v in rbs.items():
        if v.endswith(UNIVERSAL_REV_ADAPTER):
            rbs_clean[k] = v[:-len(UNIVERSAL_REV_ADAPTER)]
        else:
            rbs_clean[k] = v

    fbs_rc = {k: get_reverse_complement(v) for k, v in fbs_clean.items()}
    rbs_rc = {k: get_reverse_complement(v) for k, v in rbs_clean.items()}

    stats = {
        "total": 0, "demuxed": 0, "unassigned": 0,
        "sample_counts": {},
        "error_distribution": {i: 0 for i in range(35)}
    }
    handles = {}

    print("Анализ ридов запущен...\n")

    with open(fastq_path, "r") as handle:
        for title, seq, qual in FastqGeneralIterator(handle):
            stats["total"] += 1

            if stats["total"] % 1000 == 0:
                sys.stdout.write(
                    f"\rОбработано: {stats['total']} | Распознано: {stats['demuxed']} | Мусор: {stats['unassigned']}     ")
                sys.stdout.flush()

            res = process_read(seq, fbs_clean, rbs_clean, fbs_rc, rbs_rc)

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

    print("\n")
    return stats
