import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def analyze_assoc_file():
    """
    Детальный анализ .assoc файла
    """
    
    assoc_file = "/home/esp/data_analyze/01.06.2025_v2/data/init/gwas_results.assoc"
    
    try:
        # Загрузка данных
        print("Загрузка .assoc файла...")
        df = pd.read_csv(assoc_file, sep='\t', low_memory=False)
        
        print(f"=== ОБЩАЯ ИНФОРМАЦИЯ ===")
        print(f"Размер данных: {df.shape}")
        print(f"Колонки: {df.columns.tolist()}")
        print(f"\nПервые 5 строк:")
        print(df.head())
        
        print(f"\n=== СТАТИСТИЧЕСКАЯ СВОДКА ===")
        print(df.describe())
        
        # Анализ p-values если есть
        if 'P' in df.columns:
            print(f"\n=== АНАЛИЗ P-VALUES ===")
            p_values = df['P'].dropna()
            print(f"Всего p-values: {len(p_values)}")
            print(f"Значимые (p < 0.05): {len(p_values[p_values < 0.05])}")
            print(f"Высоко значимые (p < 0.001): {len(p_values[p_values < 0.001])}")
            print(f"Genome-wide significant (p < 5e-8): {len(p_values[p_values < 5e-8])}")
            
            # Топ значимые SNP
            top_snps = df.nsmallest(10, 'P')
            print(f"\nТоп 10 наиболее значимых SNP:")
            print(top_snps[['SNP', 'P'] if 'SNP' in df.columns else top_snps.iloc[:, [0, df.columns.get_loc('P')]]])
        
        # Анализ хромосом если есть
        if 'CHR' in df.columns:
            print(f"\n=== АНАЛИЗ ПО ХРОМОСОМАМ ===")
            chr_counts = df['CHR'].value_counts().sort_index()
            print("Количество SNP по хромосомам:")
            print(chr_counts)
        
        # Проверка на отсутствующие значения
        print(f"\n=== ОТСУТСТВУЮЩИЕ ЗНАЧЕНИЯ ===")
        missing_data = df.isnull().sum()
        print(missing_data[missing_data > 0])
        
        return df
        
    except Exception as e:
        print(f"Ошибка при анализе .assoc файла: {str(e)}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    df = analyze_assoc_file()