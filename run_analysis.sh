#!/bin/bash

echo "=== ЗАПУСК ПОЛНОГО АНАЛИЗА SNP ==="
echo "Дата: $(date)"

# Переход в рабочую директорию
cd /home/esp/data_analyze/01.06.2025_v2

# Создание директории для результатов
mkdir -p results

echo "1. Запуск основного анализа..."
python complete_analysis.py

echo "2. Создание визуализаций..."
python create_visualizations.py

echo "3. Анализ .assoc файла..."
python analyze_assoc_file.py

echo "=== АНАЛИЗ ЗАВЕРШЕН ==="
echo "Результаты сохранены в директории: results/"
ls -la results/

echo "Основные файлы результатов:"
echo "- complete_analysis_report.json - полный JSON отчет"
echo "- analysis_summary.txt - краткий текстовый отчет"
echo "- found_alzheimer_snps_detailed.csv - детальные данные по найденным SNP"
echo "- gwas_analysis_plots.pdf - графики и визуализации"
