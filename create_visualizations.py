import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

def create_gwas_visualizations():
    """
    Создание визуализаций для GWAS анализа
    """
    
    try:
        # Загрузка данных
        gwas_file = "/home/esp/data_analyze/01.06.2025_v2/data/init/gwas_results.assoc"
        df = pd.read_csv(gwas_file, sep='\t', low_memory=False)
        
        # Создание PDF с графиками
        with PdfPages('/home/esp/data_analyze/01.06.2025_v2/gwas_analysis_plots.pdf') as pdf:
            
            # Manhattan plot (если есть нужные колонки)
            if 'P' in df.columns and 'CHR' in df.columns and 'BP' in df.columns:
                plt.figure(figsize=(15, 8))
                
                # Подготовка данных для Manhattan plot
                df_plot = df.dropna(subset=['P', 'CHR', 'BP']).copy()
                df_plot['-log10P'] = -np.log10(df_plot['P'])
                
                # Цвета для хромосом
                colors = ['#1f77b4', '#ff7f0e'] * 12  # Чередующиеся цвета
                
                x_pos = 0
                x_ticks = []
                x_labels = []
                
                for i, chr_num in enumerate(sorted(df_plot['CHR'].unique())):
                    chr_data = df_plot[df_plot['CHR'] == chr_num].copy()
                    chr_data['x_pos'] = chr_data['BP'] + x_pos
                    
                    plt.scatter(chr_data['x_pos'], chr_data['-log10P'], 
                              c=colors[i % len(colors)], alpha=0.6, s=10)
                    
                    x_ticks.append(x_pos + chr_data['BP'].median())
                    x_labels.append(str(chr_num))
                    x_pos += chr_data['BP'].max() + 1000000
                
                plt.axhline(y=-np.log10(5e-8), color='red', linestyle='--', 
                           label='Genome-wide significance (5e-8)')
                plt.axhline(y=-np.log10(0.05), color='blue', linestyle='--', 
                           label='Nominal significance (0.05)')
                
                plt.xlabel('Chromosome')
                plt.ylabel('-log10(P-value)')
                plt.title('Manhattan Plot')
                plt.xticks(x_ticks, x_labels)
                plt.legend()
                plt.tight_layout()
                pdf.savefig()
                plt.close()
            
            # QQ plot
            if 'P' in df.columns:
                plt.figure(figsize=(8, 8))
                
                observed_p = df['P'].dropna().sort_values()
                expected_p = np.linspace(1/len(observed_p), 1, len(observed_p))
                
                plt.scatter(-np.log10(expected_p), -np.log10(observed_p), alpha=0.6)
                plt.plot([0, max(-np.log10(expected_p))], [0, max(-np.log10(expected_p))], 
                        'r--', label='Expected')
                
                plt.xlabel('Expected -log10(P)')
                plt.ylabel('Observed -log10(P)')
                plt.title('QQ Plot')
                plt.legend()
                plt.tight_layout()
                pdf.savefig()
                plt.close()
            
            # Гистограмма p-values
            if 'P' in df.columns:
                plt.figure(figsize=(10, 6))
                
                p_values = df['P'].dropna()
                plt.hist(p_values, bins=50, alpha=0.7, edgecolor='black')
                plt.axvline(x=0.05, color='red', linestyle='--', label='p = 0.05')
                plt.axvline(x=5e-8, color='orange', linestyle='--', label='p = 5e-8')
                
                plt.xlabel('P-value')
                plt.ylabel('Frequency')
                plt.title('Distribution of P-values')
                plt.legend()
                plt.tight_layout()
                pdf.savefig()
                plt.close()
            
            # Распределение по хромосомам
            if 'CHR' in df.columns:
                plt.figure(figsize=(12, 6))
                
                chr_counts = df['CHR'].value_counts().sort_index()
                plt.bar(chr_counts.index, chr_counts.values, alpha=0.7)
                
                plt.xlabel('Chromosome')
                plt.ylabel('Number of SNPs')
                plt.title('SNP Distribution by Chromosome')
                plt.xticks(rotation=45)
                plt.tight_layout()
                pdf.savefig()
                plt.close()
        
        print("Визуализации сохранены в: /home/esp/data_analyze/01.06.2025_v2/gwas_analysis_plots.pdf")
        
    except Exception as e:
        print(f"Ошибка при создании визуализаций: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    create_gwas_visualizations()