from utils import parse_barcodes
from demux_core import process_and_collect_stats
from visualization import plot_statistics

IN_FASTQ = "read_file/test_reads.fastq"
IN_FASTA = "read_file/barcodes.fasta"
OUT_DIR = "demux_results"


def main():
    print("Начинаем анализ...")
    fbs, rbs = parse_barcodes(IN_FASTA)
    print(f"Загружено: {len(fbs)} Forward и {len(rbs)} Reverse баркодов.")

    # Запускаем процессинг и получаем словарь со статистикой
    stats = process_and_collect_stats(IN_FASTQ, fbs, rbs, OUT_DIR)

    print(f"\nГотово! Всего ридов: {stats['total_reads']}")
    print(f"Распознано: {stats['demuxed']} | Не распознано (мусор): {stats['unassigned']}")

    # Генерируем красивые графики
    plot_statistics(stats, OUT_DIR)


if __name__ == "__main__":
    main()