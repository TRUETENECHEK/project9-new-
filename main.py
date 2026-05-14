from utils import parse_barcodes
from demux_core import run_demux
from visualization import plot_statistics  # Возвращаем отрисовку

IN_FASTQ = "read_file/test_reads.fastq"
IN_FASTA = "read_file/barcodes.fasta"
OUT_DIR = "demux_results"


def main():
    print("Инициализация пайплайна...")
    fbs, rbs = parse_barcodes(IN_FASTA)
    print(f"Загружено: {len(fbs)} FWD и {len(rbs)} REV баркодов.")

    # Запуск движка
    stats = run_demux(IN_FASTQ, fbs, rbs, OUT_DIR)

    print(f"\n✅ Завершено! Всего ридов: {stats['total']}")
    print(f"🧬 Демультиплексировано: {stats['demuxed']}")
    print(f"🗑️ Мусор/Не распознано: {stats['unassigned']}")

    # Запуск твоей визулизации
    print("Генерация графиков...")
    plot_statistics(stats, OUT_DIR)


if __name__ == "__main__":
    main()