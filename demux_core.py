import regex
import os
from Bio import SeqIO
from Bio.Seq import Seq
from typing import Dict, Tuple, Optional

# --- КОНСТАНТЫ ---
# Универсальные части FB и RB, которые служат "якорями"
UNIVERSAL_FWD_ADAPTER = "TGTTGGTGCAGATATTGCGGG"
UNIVERSAL_REV_ADAPTER = "TTGCCTGTCGCTCTATCTTCAAA"

# Параметры ошибок для поиска якорей (очень лояльные)
# Настройка {s<=6,i<=2,d<=2} позволяет до 10 ошибок суммарно (6 замен, 2 вставки, 2 делеции)
ANCHOR_ERRORS = "{s<=6,i<=2,d<=2}" 

# Параметры ошибок для поиска баркодов во флангах (строгие)
BARCODE_ERRORS = "{e<=1}"

def get_reverse_complement(seq: str) -> str:
    """Возвращает reverse complement последовательности."""
    return str(Seq(seq).reverse_complement())

def find_complex_adapters(sequence: str) -> Tuple[Optional[str], Optional[regex.Match], Optional[regex.Match], int]:
    """
    Ищет универсальные адаптеры (якори) в риде с помощью regex.
    Возвращает: ориентацию ('PLUS'/'MINUS'), Match для первого якоря, Match для второго якоря и суммарное кол-во ошибок.
    """
    fwd_rc = get_reverse_complement(UNIVERSAL_FWD_ADAPTER)
    rev_rc = get_reverse_complement(UNIVERSAL_REV_ADAPTER)

    # 1. PLUS-ориентация: [FB][FWD_ANCHOR] ... [REV_RC_ANCHOR][RB_RC]
    m_fwd = regex.search(f"(?b)({UNIVERSAL_FWD_ADAPTER}){ANCHOR_ERRORS}", sequence)
    m_rev_rc = regex.search(f"(?b)({rev_rc}){ANCHOR_ERRORS}", sequence)

    if m_fwd and m_rev_rc and m_fwd.start() < m_rev_rc.start():
        total_errs = sum(m_fwd.fuzzy_counts) + sum(m_rev_rc.fuzzy_counts)
        return "PLUS", m_fwd, m_rev_rc, total_errs

    # 2. MINUS-ориентация: [RB][REV_ANCHOR] ... [FWD_RC_ANCHOR][FB_RC]
    m_rev = regex.search(f"(?b)({UNIVERSAL_REV_ADAPTER}){ANCHOR_ERRORS}", sequence)
    m_fwd_rc = regex.search(f"(?b)({fwd_rc}){ANCHOR_ERRORS}", sequence)

    if m_rev and m_fwd_rc and m_rev.start() < m_fwd_rc.start():
        total_errs = sum(m_rev.fuzzy_counts) + sum(m_fwd_rc.fuzzy_counts)
        return "MINUS", m_rev, m_fwd_rc, total_errs

    return None, None, None, 0

def find_adapters_and_match_barcodes(sequence: str, fbs_clean: Dict[str, str], rbs_clean: Dict[str, str]) -> Optional[dict]:
    """
    Реализует метод 'Якорей' шаг за шагом:
    1. Находит универсальные адаптеры для определения координат и ориентации.
    2. Выделяет короткие фланки (35-40 нуклеотидов) снаружи от найденных адаптеров.
    3. Ищет специфичные баркоды во флангах с минимальным допуском ошибок.
    """
    orientation, m1, m2, anchor_errors = find_complex_adapters(sequence)
    if not orientation:
        return None

    if orientation == "PLUS":
        # Фланки: берем область снаружи от якорей с небольшим перекрытием (5 нуклеотидов)
        fb_flank = sequence[max(0, m1.start() - 40) : m1.start() + 5]
        rb_flank = sequence[m2.end() - 5 : min(len(sequence), m2.end() + 40)]
        
        # Поиск баркодов во флангах
        found_fb = next((id for id, s in fbs_clean.items() if regex.search(f"(?b)({s}){BARCODE_ERRORS}", fb_flank)), None)
        found_rb = next((id for id, s in rbs_clean.items() if regex.search(f"(?b)({get_reverse_complement(s)}){BARCODE_ERRORS}", rb_flank)), None)
    else:
        # MINUS: [RB_clean][REV_ANCHOR] ... [FWD_RC_ANCHOR][FB_clean_RC]
        rb_flank = sequence[max(0, m1.start() - 40) : m1.start() + 5]
        fb_flank = sequence[m2.end() - 5 : min(len(sequence), m2.end() + 40)]
        
        found_rb = next((id for id, s in rbs_clean.items() if regex.search(f"(?b)({s}){BARCODE_ERRORS}", rb_flank)), None)
        found_fb = next((id for id, s in fbs_clean.items() if regex.search(f"(?b)({get_reverse_complement(s)}){BARCODE_ERRORS}", fb_flank)), None)

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
    """
    Основной пайплайн: читает FASTQ, демультиплексирует риды и собирает статистику для визуализации.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Очищаем баркоды от универсальных частей (якорей), так как ищем их во флангах
    fbs_clean = {id: seq.replace(UNIVERSAL_FWD_ADAPTER, "") for id, seq in fbs.items()}
    rbs_clean = {id: seq.replace(UNIVERSAL_REV_ADAPTER, "") for id, seq in rbs.items()}

    stats = {
        "total_reads": 0,
        "demuxed": 0,
        "unassigned": 0,
        "sample_counts": {},
        "error_distribution": {i: 0 for i in range(16)} # Статистика ошибок (от 0 до 15+)
    }

    # Словарь для хранения открытых файловых дескрипторов
    handles = {}

    print("Обработка ридов...")
    for record in SeqIO.parse(fastq_path, "fastq"):
        stats["total_reads"] += 1
        
        # Визуальный прогресс каждые 100 ридов
        if stats["total_reads"] % 100 == 0:
            print(f"Обработано {stats['total_reads']} ридов...", end="\r")

        res = find_adapters_and_match_barcodes(str(record.seq), fbs_clean, rbs_clean)

        if res:
            sample_id = res["sample_id"]
            stats["demuxed"] += 1
            stats["sample_counts"][sample_id] = stats["sample_counts"].get(sample_id, 0) + 1
            
            # Собираем распределение ошибок для графика
            err_idx = min(res["total_errors"], 15)
            stats["error_distribution"][err_idx] += 1

            # Сохраняем рид в соответствующий файл
            if sample_id not in handles:
                handles[sample_id] = open(os.path.join(output_dir, f"{sample_id}.fastq"), "w")
            SeqIO.write(record, handles[sample_id], "fastq")
        else:
            stats["unassigned"] += 1

    # Закрываем все открытые файлы
    for h in handles.values():
        h.close()
    
    print(f"\nОбработка завершена. Всего: {stats['total_reads']}")
    return stats
