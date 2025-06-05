import pandas as pd
import glob
import os
import subprocess

# Пути
base_dir = os.path.abspath(os.path.dirname(__file__))
data_dir = os.path.join(base_dir, '../data', 'init')
out_dir = os.path.join(base_dir, '../data', 'output', 'first_step')
os.makedirs(out_dir, exist_ok=True)

pheno_path = os.path.join(data_dir, 'merged_all.phenotype')
plink2_path = os.path.join(base_dir, '../plink2_linux_avx2_20250411', 'plink2')

# 1. Обновление .fam-файлов с фенотипами
print('Обновление .fam-файлов с фенотипами...')
pheno = pd.read_csv(pheno_path, delim_whitespace=True)
pheno_dict = dict(zip(pheno['IID'], pheno['PHENO']))

fam_files = glob.glob(os.path.join(data_dir, '*', '*.fam'))
for fam_path in fam_files:
    fam = pd.read_csv(fam_path, delim_whitespace=True, header=None)
    fam[5] = fam[1].map(pheno_dict).fillna(-9).astype(int)
    fam.to_csv(fam_path, sep=' ', header=False, index=False)
    print(f'Обновлен: {fam_path}')

# 2. Объединение .bed файлов через PLINK 1.9 (так как PLINK 2.0 не поддерживает merge)
print('Объединение .bed-файлов с помощью PLINK 1.9...')
bfiles = [
    os.path.join(data_dir, '1-zapusk', '1-zapusk_binary'),
    os.path.join(data_dir, '2-zapusk', '2-zapusk_binary'),
    os.path.join(data_dir, '3.1-zapusk', '3-zapusk_binary'),
    os.path.join(data_dir, '3.2-zapusk', '13-sample_binary'),
]
bmerge_list_path = os.path.join(base_dir, 'bmerge_list.txt')
with open(bmerge_list_path, 'w') as f:
    for b in bfiles[1:]:
        f.write(b + '\n')
merged_prefix_bed = os.path.join(out_dir, 'merged_all_bed')

plink1_9_path = '../plink_1.9/plink'
subprocess.run([
    plink1_9_path,
    '--bfile', bfiles[0],
    '--merge-list', bmerge_list_path,
    '--make-bed',
    '--allow-no-sex',
    '--out', merged_prefix_bed
], check=True)

# 3. Конвертация merged_all_bed в форматы PLINK 2.0
print('Конвертация в pgen формат (PLINK 2.0)...')
merged_prefix_pgen = os.path.join(out_dir, 'merged_all')
subprocess.run([
    plink2_path,
    '--bfile', merged_prefix_bed,
    '--make-pgen',
    '--out', merged_prefix_pgen
], check=True)

# 4. GWAS-анализ через PLINK 2.0 (--glm)
print('Запуск GWAS-анализа (PLINK 2.0)...')
plink2_out = os.path.join(out_dir, 'gwas_results')
subprocess.run([
    plink2_path,
    '--pfile', merged_prefix_pgen,
    '--pheno', pheno_path,
    '--pheno-name', 'PHENO',
    '--glm', 'allow-no-covars',
    '--out', plink2_out
], check=True)

print('✅ Готово! Все результаты сохранены в', out_dir)
