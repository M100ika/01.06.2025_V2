import pandas as pd

# Путь к исходному файлу
input_path = "/home/esp/data_analyze/01.06.2025_v2/data/init/gwas_results.assoc"

# Путь для сохранения CSV
output_path = "/home/esp/data_analyze/01.06.2025_v2/data/output/gwas_results.csv"

# Загрузка файла
# Предположим, что разделитель — пробел или табуляция, и нужно обработать заголовки
df = pd.read_csv(input_path, delim_whitespace=True)

# Сохранение в CSV с разделителем ';'
df.to_csv(output_path, sep=';', index=False)

print(f"Файл успешно сохранён: {output_path}")