import matplotlib.pyplot as plt
import seaborn as sns
import os

def plot_statistics(stats: dict, output_dir: str):
    sns.set_theme(style="whitegrid")
    fig = plt.figure(figsize=(15, 10))

    # 1. Распределение по сэмплами (Bar chart)
    ax1 = plt.subplot(2, 2, (1, 2))
    samples = list(stats["sample_counts"].keys())
    counts = list(stats["sample_counts"].values())
    sns.barplot(x=samples, y=counts, ax=ax1, palette="viridis")
    ax1.set_title("Количество ридов по образцам")
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')

    # 2. Успешность демультиплексации (Pie chart)
    ax2 = plt.subplot(2, 2, 3)
    labels = ['Распознано', 'Не распознано']
    sizes = [stats["demuxed"], stats["unassigned"]]
    ax2.pie(sizes, labels=labels, autopct='%1.1f%%', colors=['#4CAF50', '#F44336'])
    ax2.set_title("Общая успешность пайплайна")

    # 3. Распределение ошибок в адаптерах
    ax3 = plt.subplot(2, 2, 4)
    err_labels = list(stats["error_distribution"].keys())
    err_counts = list(stats["error_distribution"].values())
    sns.barplot(x=err_labels, y=err_counts, ax=ax3, palette="mako")
    ax3.set_title("Количество допущенных ошибок при поиске")
    ax3.set_xlabel("Суммарно ошибок (mismatch/indel)")

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "demux_report.png"), dpi=300)
    print(f"Отчет сохранен в {output_dir}/demux_report.png")