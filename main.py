import pandas as pd
import os
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

    # Генерация Excel отчета
    print("Генерация Excel отчета...")
    report_data = []
    for s_id, count in stats['sample_counts'].items():
        if '_' in s_id:
            fb, rb = s_id.split('_', 1)
        else:
            fb, rb = s_id, ""
        
        total_len = stats.get('sample_lengths', {}).get(s_id, 0)
        avg_len = round(total_len / count) if count > 0 else 0
        
        report_data.append({
            'Barcode1': fb,
            'Barcode2': rb,
            'SampleID': s_id,
            'Trimmed_Reads_Count': count,
            'Average_Trimmed_Length': avg_len
        })
    
    if report_data:
        df = pd.DataFrame(report_data)
        # Сортируем по количеству ридов по убыванию
        df = df.sort_values(by='Trimmed_Reads_Count', ascending=False)
        excel_path = os.path.join(OUT_DIR, 'demultiplexing_report.xlsx')

        df.to_excel(excel_path, index=False)
        print(f"Отчет сохранен в {excel_path}")

    if stats.get('benchmark_data'):
        print("Генерация Benchmark Excel отчета...")
        df_bench = pd.DataFrame(stats['benchmark_data'])
        bench_path = os.path.join(OUT_DIR, 'benchmark_trimmed_sequences.xlsx')
        df_bench.to_excel(bench_path, index=False)
        print(f"Бенчмарк отчет сохранен в {bench_path}")



if __name__ == "__main__":
    main()
