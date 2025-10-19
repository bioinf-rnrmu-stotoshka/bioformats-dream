import pandas as pd
from pathlib import Path
from typing import Iterator, Dict, List
from abstract import GenomicDataReader
from record import VariantRecord


class VcfReader(GenomicDataReader):
    """
    Ридер для VCF файлов
    """

    def __init__(self, filepath: str | Path):
        super().__init__(filepath)
        self.headers = []
        self.samples = []
        self._variants_count = 0

    def _parse_header(self):
        """
        Парсит заголовок VCF файла
        """
        self.headers = []
        self.samples = []
        
        if not self.file:
            return

        # Сохраняем начальную позицию
        current_pos = self.file.tell()
        self.file.seek(0)
        
        line = self.file.readline()
        while line and line.startswith('#'):
            self.headers.append(line.strip())
            
            if line.startswith('#CHROM'):
                # Парсим строку с названиями колонок
                parts = line.strip().split('\t')
                if len(parts) >= 10:  # Минимум: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT + samples
                    self.samples = parts[9:]
            
            line = self.file.readline()
        
        # Возвращаемся к началу данных
        self.file.seek(current_pos)
        self._header_parsed = True

    def read(self) -> Iterator[VariantRecord]:
        """
        Читает варианты из VCF файла
        """
        if not self.file:
            raise ValueError("Файл не открыт")
        
        # Переходим к началу данных (после заголовков)
        self.file.seek(0)
        line = self.file.readline()
        while line and line.startswith('#'):
            line = self.file.readline()
        
        # Читаем данные
        self._variants_count = 0
        while line:
            line = line.strip()
            if line and not line.startswith('#'):
                try:
                    variant = self._parse_variant_line(line)
                    if variant:
                        self._variants_count += 1
                        yield variant
                except Exception as e:
                    print(f"Ошибка парсинга строки: {line[:50]}... - {e}")
            
            line = self.file.readline()

    def _parse_variant_line(self, line: str) -> VariantRecord | None:
        """
        Парсит одну строку с вариантом
        """
        parts = line.split('\t')
        if len(parts) < 8:
            return None
        
        try:
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            
            # Парсим INFO поле
            info_dict = {}
            if parts[7] != '.':
                for info_item in parts[7].split(';'):
                    if '=' in info_item:
                        key, value = info_item.split('=', 1)
                        info_dict[key] = value
                    else:
                        info_dict[info_item] = True
            
            return VariantRecord(chrom, pos, ref, alt, info_dict)
            
        except (ValueError, IndexError) as e:
            print(f"Ошибка парсинга варианта: {e}")
            return None

    def get_statistics(self) -> Dict:
        """
        Возвращает статистику по файлу
        """
        # Пересчитываем количество вариантов
        count = self._count_variants()
        return {
            'total_variants': count,
            'samples_count': len(self.samples),
            'header_lines': len(self.headers)
        }

    def _count_variants(self) -> int:
        """Считает общее количество вариантов в файле"""
        if not self.file:
            return 0
        
        current_pos = self.file.tell()
        self.file.seek(0)
        
        # Пропускаем заголовки
        line = self.file.readline()
        while line and line.startswith('#'):
            line = self.file.readline()
        
        # Считаем варианты
        count = 0
        while line:
            line = line.strip()
            if line and not line.startswith('#'): 
                count += 1
            line = self.file.readline()
        
        # Возвращаемся к исходной позиции
        self.file.seek(current_pos)
        return count

    def get_chromosomes(self) -> List[str]:
        """
        Возвращает список хромосом в файле
        """
        chromosomes = set()
        
        if not self.file:
            return list(chromosomes)
        
        current_pos = self.file.tell()
        self.file.seek(0)
        
        # Пропускаем заголовки
        line = self.file.readline()
        while line and line.startswith('#'):
            line = self.file.readline()
        
        # Собираем хромосомы
        while line:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split('\t')
                if parts:
                    chromosomes.add(parts[0])
            line = self.file.readline()
        
        # Возвращаемся к исходной позиции
        self.file.seek(current_pos)
        
        return sorted(list(chromosomes))

    def get_region_stats(self) -> pd.DataFrame:
        """
        Возвращает статистику по регионам (хромосомам)
        """
        region_counts = {}
        
        if not self.file:
            return pd.DataFrame(columns=['chromosome', 'variant_count'])
        
        current_pos = self.file.tell()
        self.file.seek(0)
        
        # Пропускаем заголовки
        line = self.file.readline()
        while line and line.startswith('#'):
            line = self.file.readline()
        
        # Считаем варианты по хромосомам
        while line:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split('\t')
                if parts:
                    chrom = parts[0]
                    region_counts[chrom] = region_counts.get(chrom, 0) + 1
            line = self.file.readline()
        
        # Возвращаемся к исходной позиции
        self.file.seek(current_pos)
        
        # Создаем DataFrame
        data = []
        for chrom, count in region_counts.items():
            data.append({'chromosome': chrom, 'variant_count': count})
        
        return pd.DataFrame(data)

    def get_variant_type_stats(self) -> pd.DataFrame:
        """
        Возвращает статистику по типам вариантов
        """
        type_counts = {
            'SNV': 0,
            'Insertion': 0,
            'Deletion': 0,
            'Complex': 0
        }
        
        if not self.file:
            return pd.DataFrame(columns=['variant_type', 'count'])
        
        current_pos = self.file.tell()
        self.file.seek(0)
        
        # Пропускаем заголовки
        line = self.file.readline()
        while line and line.startswith('#'):
            line = self.file.readline()
        
        # Анализируем типы вариантов
        while line:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 5:
                    ref = parts[3]
                    alt = parts[4]
                    
                    # Определяем тип варианта
                    if len(ref) == 1 and len(alt) == 1:
                        type_counts['SNV'] += 1
                    elif len(ref) < len(alt):
                        type_counts['Insertion'] += 1
                    elif len(ref) > len(alt):
                        type_counts['Deletion'] += 1
                    else:
                        type_counts['Complex'] += 1
            
            line = self.file.readline()
        
        # Возвращаемся к исходной позиции
        self.file.seek(current_pos)
        
        # Создаем DataFrame
        data = []
        for var_type, count in type_counts.items():
            if count > 0:
                data.append({'variant_type': var_type, 'count': count})
        
        return pd.DataFrame(data)

    def reset(self):
        """
        Сбрасывает позицию чтения файла
        """
        if self.file:
            self.file.seek(0)
            self._parse_header()

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Закрывает файл при выходе из контекста
        """
        super().__exit__(exc_type, exc_val, exc_tb)
        self._header_parsed = False
        self.headers = []
        self.samples = []
        self._variants_count = 0
