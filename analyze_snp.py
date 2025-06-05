import pandas as pd
import numpy as np
import os

def analyze_snp_data():
    """
    Анализ SNP данных: поиск SNP из Excel файла в GWAS результатах
    """
    
    # Пути к файлам
    excel_file = "/home/esp/data_analyze/01.06.2025_v2/data/init/alleli_alz.xlsx"
    gwas_file = "/home/esp/data_analyze/01.06.2025_v2/data/init/gwas_results.assoc"
    
    try:
        # Чтение Excel файла с аллелями
        print("Загрузка данных из Excel файла...")
        df_excel = pd.read_excel(excel_file)
        print(f"Структура Excel файла:")
        print(df_excel.head())
        print(f"Колонки: {df_excel.columns.tolist()}")
        
        # Попытка найти колонку с SNP ID
        snp_columns = [col for col in df_excel.columns if 'snp' in col.lower() or 'rs' in col.lower()]
        if not snp_columns:
            # Если не найдены специфичные колонки, показываем первые несколько строк для ручного выбора
            print("Автоматически не удалось определить колонку с SNP. Доступные колонки:")
            for i, col in enumerate(df_excel.columns):
                print(f"{i}: {col}")
            return
        
        snp_column = snp_columns[0]
        snp_list = df_excel[snp_column].dropna().unique().tolist()
        print(f"Найдено {len(snp_list)} уникальных SNP в колонке '{snp_column}'")
        
        # Чтение GWAS результатов
        print("\nЗагрузка GWAS результатов...")
        df_gwas = pd.read_csv(gwas_file, sep='\t', low_memory=False)
        print(f"Структура GWAS файла:")
        print(df_gwas.head())
        print(f"Колонки: {df_gwas.columns.tolist()}")
        
        # Поиск колонки с SNP ID в GWAS файле
        gwas_snp_columns = [col for col in df_gwas.columns if 'snp' in col.lower() or 'rs' in col.lower() or 'id' in col.lower()]
        if not gwas_snp_columns:
            print("Не удалось найти колонку с SNP в GWAS файле")
            return
        
        gwas_snp_column = gwas_snp_columns[0]
        print(f"Используется колонка '{gwas_snp_column}' для поиска SNP в GWAS данных")
        
        # Поиск пересечений
        print("\nПоиск пересечений...")
        found_snps = []
        not_found_snps = []
        
        for snp in snp_list:
            if pd.isna(snp):
                continue
            snp_str = str(snp).strip()
            if snp_str in df_gwas[gwas_snp_column].values:
                found_snps.append(snp_str)
            else:
                not_found_snps.append(snp_str)
        
        # Результаты поиска
        print(f"\n=== РЕЗУЛЬТАТЫ АНАЛИЗА ===")
        print(f"Всего SNP для поиска: {len(snp_list)}")
        print(f"Найдено в GWAS: {len(found_snps)}")
        print(f"Не найдено в GWAS: {len(not_found_snps)}")
        
        if found_snps:
            print(f"\nНайденные SNP:")
            for snp in found_snps[:10]:  # Показываем первые 10
                print(f"  {snp}")
            if len(found_snps) > 10:
                print(f"  ... и еще {len(found_snps) - 10}")
        
        if not_found_snps:
            print(f"\nНе найденные SNP (первые 10):")
            for snp in not_found_snps[:10]:
                print(f"  {snp}")
            if len(not_found_snps) > 10:
                print(f"  ... и еще {len(not_found_snps) - 10}")
        
        # Создание детального отчета
        if found_snps:
            print("\nСоздание детального отчета...")
            detailed_results = df_gwas[df_gwas[gwas_snp_column].isin(found_snps)].copy()
            
            # Сохранение результатов
            output_file = "/home/esp/data_analyze/01.06.2025_v2/found_alzheimer_snps.csv"
            detailed_results.to_csv(output_file, index=False)
            print(f"Детальные результаты сохранены в: {output_file}")
            
            # Статистика по найденным SNP
            if 'P' in detailed_results.columns:
                significant_snps = detailed_results[detailed_results['P'] < 0.05]
                print(f"SNP с p-value < 0.05: {len(significant_snps)}")
                
                highly_significant = detailed_results[detailed_results['P'] < 0.001]
                print(f"SNP с p-value < 0.001: {len(highly_significant)}")
        
        # Создание сводного отчета
        summary_report = {
            'total_snps_searched': len(snp_list),
            'found_in_gwas': len(found_snps),
            'not_found_in_gwas': len(not_found_snps),
            'found_snps_list': found_snps,
            'not_found_snps_list': not_found_snps
        }
        
        # Сохранение сводного отчета
        import json
        summary_file = "/home/esp/data_analyze/01.06.2025_v2/snp_analysis_summary.json"
        with open(summary_file, 'w', encoding='utf-8') as f:
            json.dump(summary_report, f, ensure_ascii=False, indent=2)
        print(f"\nСводный отчет сохранен в: {summary_file}")
        
    except Exception as e:
        print(f"Ошибка при анализе: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    analyze_snp_data()