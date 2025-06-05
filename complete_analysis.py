import pandas as pd
import numpy as np
import json
import os
from datetime import datetime

def complete_snp_analysis():
    """
    Полный анализ SNP данных с созданием итогового отчета
    """
    
    print("=== НАЧАЛО ПОЛНОГО АНАЛИЗА SNP ===")
    print(f"Время начала: {datetime.now()}")
    
    # Пути к файлам
    excel_file = "/home/esp/data_analyze/01.06.2025_v2/data/init/Аллели по болезни Альцгеймера .xlsx"
    gwas_file = "/home/esp/data_analyze/01.06.2025_v2/data/init/gwas_results.assoc"
    output_dir = "/home/esp/data_analyze/01.06.2025_v2/results"
    
    # Создание директории для результатов
    os.makedirs(output_dir, exist_ok=True)
    
    results = {
        'analysis_date': datetime.now().isoformat(),
        'input_files': {
            'excel_file': excel_file,
            'gwas_file': gwas_file
        },
        'summary': {},
        'detailed_results': {}
    }
    
    try:
        # 1. Анализ Excel файла
        print("\n1. Анализ Excel файла с аллелями Альцгеймера...")
        df_excel = pd.read_excel(excel_file)
        
        print(f"   Размер Excel файла: {df_excel.shape}")
        print(f"   Колонки: {df_excel.columns.tolist()}")
        
        # Поиск колонки с SNP
        snp_columns = [col for col in df_excel.columns 
                      if any(keyword in col.lower() for keyword in ['snp', 'rs', 'variant', 'id'])]
        
        if not snp_columns:
            print("   ВНИМАНИЕ: Не найдена колонка с SNP ID. Используется первая колонка.")
            snp_column = df_excel.columns[0]
        else:
            snp_column = snp_columns[0]
        
        print(f"   Используется колонка: '{snp_column}'")
        
        # Извлечение списка SNP
        snp_list = df_excel[snp_column].dropna().astype(str).str.strip().unique().tolist()
        snp_list = [snp for snp in snp_list if snp and snp != 'nan']
        
        results['summary']['total_snps_from_excel'] = len(snp_list)
        print(f"   Найдено уникальных SNP: {len(snp_list)}")
        
        # 2. Анализ GWAS файла
        print("\n2. Анализ GWAS результатов...")
        df_gwas = pd.read_csv(gwas_file, sep='\t', low_memory=False)
        
        print(f"   Размер GWAS файла: {df_gwas.shape}")
        print(f"   Колонки: {df_gwas.columns.tolist()}")
        
        # Поиск колонки с SNP в GWAS
        gwas_snp_columns = [col for col in df_gwas.columns 
                           if any(keyword in col.lower() for keyword in ['snp', 'rs', 'id', 'variant'])]
        
        if not gwas_snp_columns:
            print("   ВНИМАНИЕ: Не найдена колонка с SNP ID в GWAS. Используется первая колонка.")
            gwas_snp_column = df_gwas.columns[0]
        else:
            gwas_snp_column = gwas_snp_columns[0]
        
        print(f"   Используется колонка: '{gwas_snp_column}'")
        
        results['summary']['total_snps_in_gwas'] = len(df_gwas)
        
        # 3. Поиск пересечений
        print("\n3. Поиск пересечений...")
        
        # Создание множества для быстрого поиска
        gwas_snps_set = set(df_gwas[gwas_snp_column].astype(str).str.strip())
        
        found_snps = []
        not_found_snps = []
        
        for snp in snp_list:
            if snp in gwas_snps_set:
                found_snps.append(snp)
            else:
                not_found_snps.append(snp)
        
        results['summary']['found_snps_count'] = len(found_snps)
        results['summary']['not_found_snps_count'] = len(not_found_snps)
        results['summary']['match_percentage'] = (len(found_snps) / len(snp_list)) * 100 if snp_list else 0
        
        print(f"   Найдено в GWAS: {len(found_snps)} ({results['summary']['match_percentage']:.1f}%)")
        print(f"   Не найдено: {len(not_found_snps)}")
        
        # 4. Детальный анализ найденных SNP
        if found_snps:
            print("\n4. Детальный анализ найденных SNP...")
            
            found_gwas_data = df_gwas[df_gwas[gwas_snp_column].isin(found_snps)].copy()
            
            # Статистический анализ
            if 'P' in found_gwas_data.columns:
                p_values = found_gwas_data['P'].dropna()
                
                stats = {
                    'total_with_pvalue': len(p_values),
                    'significant_005': len(p_values[p_values < 0.05]),
                    'significant_001': len(p_values[p_values < 0.001]),
                    'genome_wide_significant': len(p_values[p_values < 5e-8]),
                    'min_pvalue': float(p_values.min()) if len(p_values) > 0 else None,
                    'median_pvalue': float(p_values.median()) if len(p_values) > 0 else None
                }
                
                results['detailed_results']['p_value_statistics'] = stats
                
                print(f"   SNP с p-value < 0.05: {stats['significant_005']}")
                print(f"   SNP с p-value < 0.001: {stats['significant_001']}")
                print(f"   Genome-wide significant (p < 5e-8): {stats['genome_wide_significant']}")
                
                # Топ значимые SNP
                if len(p_values) > 0:
                    top_significant = found_gwas_data.nsmallest(min(10, len(found_gwas_data)), 'P')
                    results['detailed_results']['top_10_significant'] = top_significant.to_dict('records')
            
            # Сохранение детальных результатов
            output_file = os.path.join(output_dir, 'found_alzheimer_snps_detailed.csv')
            found_gwas_data.to_csv(output_file, index=False)
            print(f"   Детальные результаты сохранены: {output_file}")
        
        # 5. Сохранение всех результатов
        print("\n5. Сохранение результатов...")
        
        # Списки SNP
        results['detailed_results']['found_snps'] = found_snps
        results['detailed_results']['not_found_snps'] = not_found_snps
        
        # JSON отчет
        json_file = os.path.join(output_dir, 'complete_analysis_report.json')
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(results, f, ensure_ascii=False, indent=2, default=str)
        
        # Текстовый отчет
        txt_file = os.path.join(output_dir, 'analysis_summary.txt')
        with open(txt_file, 'w', encoding='utf-8') as f:
            f.write("=== ОТЧЕТ ПО АНАЛИЗУ SNP АЛЬЦГЕЙМЕРА ===\n\n")
            f.write(f"Дата анализа: {results['analysis_date']}\n\n")
            f.write("ИСХОДНЫЕ ДАННЫЕ:\n")
            f.write(f"- Excel файл: {excel_file}\n")
            f.write(f"- GWAS файл: {gwas_file}\n\n")
            f.write("РЕЗУЛЬТАТЫ:\n")
            f.write(f"- Всего SNP из Excel: {results['summary']['total_snps_from_excel']}\n")
            f.write(f"- Найдено в GWAS: {results['summary']['found_snps_count']}\n")
            f.write(f"- Не найдено: {results['summary']['not_found_snps_count']}\n")
            f.write(f"- Процент совпадений: {results['summary']['match_percentage']:.1f}%\n\n")
            

            if 'p_value_statistics' in results['detailed_results']:
                stats = results['detailed_results']['p_value_statistics']
                f.write("СТАТИСТИКА P-VALUES:\n")
                f.write(f"- Значимые (p < 0.05): {stats['significant_005']}\n")
                f.write(f"- Высоко значимые (p < 0.001): {stats['significant_001']}\n")
                f.write(f"- Genome-wide significant (p < 5e-8): {stats['genome_wide_significant']}\n")
                f.write(f"- Минимальное p-value: {stats['min_pvalue']}\n")
                f.write(f"- Медианное p-value: {stats['median_pvalue']}\n\n")
            
            f.write("НАЙДЕННЫЕ SNP:\n")
            for snp in found_snps[:20]:  # Первые 20
                f.write(f"- {snp}\n")
            if len(found_snps) > 20:
                f.write(f"... и еще {len(found_snps) - 20} SNP\n")
            
            f.write("\nНЕ НАЙДЕННЫЕ SNP:\n")
            for snp in not_found_snps[:20]:  # Первые 20
                f.write(f"- {snp}\n")
            if len(not_found_snps) > 20:
                f.write(f"... и еще {len(not_found_snps) - 20} SNP\n")
        
        print(f"   JSON отчет: {json_file}")
        print(f"   Текстовый отчет: {txt_file}")
        
        print(f"\n=== АНАЛИЗ ЗАВЕРШЕН ===")
        print(f"Время завершения: {datetime.now()}")
        print(f"Все результаты сохранены в: {output_dir}")
        
        return results
        
    except Exception as e:
        print(f"ОШИБКА при анализе: {str(e)}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    results = complete_snp_analysis()