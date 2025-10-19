import pandas as pd
from pathlib import Path
from typing import Iterator, Dict, List
from abstract_vcf import GenomicDataReader
from record_vcf import VariantRecordVcf


class VcfReader(GenomicDataReader):
    """
    Реализация ридера для файлов в формате VCF (Variant Call Format).
    
    VCF - стандартный формат для хранения геномных вариаций, содержащий 
    информацию о хромосоме, позиции, референсном и альтернативном аллелях,
    качестве вызова и дополнительных аннотациях.
    
    Наследует от абстрактного класса GenomicDataReader и реализует
    специфичную для VCF логику парсинга и анализа данных.
    
    Атрибуты:
        headers (List[str]): Список строк заголовка VCF файла
        samples (List[str]): Список образцов (индивидуумов) из VCF файла
        _variants_count (int): Счетчик вариантов (для внутреннего использования)
        _header_parsed (bool): Флаг, указывающий был ли распарсен заголовок
    """

    def __init__(self, filepath: str | Path):
        """
        Инициализирует VCF ридер.
        
        Args:
            filepath (str | Path): Путь к VCF файлу. Может быть строкой или объектом Path.
            
        Инициализирует:
            - headers: пустой список для хранения строк заголовка
            - samples: пустой список для имен образцов
            - _variants_count: счетчик вариантов (0)
            - _header_parsed: флаг парсинга заголовка (False)
        """
        super().__init__(filepath)
        self.headers = []
        self.samples = []
        self._variants_count = 0

    def _parse_header(self):
        """
        Парсит заголовочную часть VCF файла.
        
        VCF заголовок состоит из строк, начинающихся с '#' и содержит:
        - Мета-информацию (##-строки)
        - Строку с названиями колонок (#CHROM...)
        - Информацию о формате и образцах
        
        Алгоритм:
        1. Сохраняет текущую позицию в файле
        2. Читает файл с начала до первой не-заголовочной строки
        3. Извлекает имена образцов из строки #CHROM
        4. Восстанавливает исходную позицию чтения
        
        Устанавливает:
            self.headers: все строки заголовка
            self.samples: список образцов (если есть)
            self._header_parsed: True после завершения парсинга
        """
        self.headers = []
        self.samples = []
        
        if not self.file:
            return

        # Сохраняем начальную позицию для восстановления после парсинга
        current_pos = self.file.tell()
        self.file.seek(0)
        
        line = self.file.readline()
        while line and line.startswith('#'):
            self.headers.append(line.strip())
            
            if line.startswith('#CHROM'):
                # Парсим строку с названиями колонок
                # Формат: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT [SAMPLE1 ...]
                parts = line.strip().split('\t')
                if len(parts) >= 10:  # Минимум: обязательные поля + samples
                    self.samples = parts[9:]  # Образцы начинаются с 9-й колонки
            
            line = self.file.readline()
        
        # Возвращаемся к начальной позиции чтения данных
        self.file.seek(current_pos)
        self._header_parsed = True

    def read(self) -> Iterator[VariantRecordVcf]:
        """
        Основной метод чтения вариантов из VCF файла.
        
        Генератор, который последовательно читает и парсит варианты из файла,
        возвращая объекты VariantRecord для каждой валидной строки с вариантом.
        
        Yields:
            VariantRecord: Объект, содержащий информацию о геномном варианте
            
        Raises:
            ValueError: Если файл не был открыт перед чтением
            
        Процесс:
        1. Пропускает все заголовочные строки
        2. Читает файл построчно
        3. Для каждой не-пустой и не-заголовочной строки вызывает _parse_variant_line()
        4. Увеличивает счетчик вариантов для каждого успешно распарсенного варианта
        """
        if not self.file:
            raise ValueError("Файл не открыт")
        
        # Переходим к началу данных (после заголовков)
        self.file.seek(0)
        line = self.file.readline()
        while line and line.startswith('#'):
            line = self.file.readline()
        
        # Читаем и парсим данные о вариантах
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

    def _parse_variant_line(self, line: str) -> VariantRecordVcf | None:
        """
        Парсит одну строку с вариантом из VCF файла.
        
        VCF строки вариантов содержат минимум 8 обязательных полей:
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
        
        Args:
            line (str): Строка из VCF файла, содержащая описание варианта
            
        Returns:
            VariantRecord | None: Объект варианта или None при ошибке парсинга
            
        Обрабатываемые поля:
            - CHROM: хромосома (str)
            - POS: позиция (int)
            - REF: референсный аллель (str)
            - ALT: альтернативный аллель (str)
            - INFO: поле с дополнительной информацией (парсится в dict)
        """
        parts = line.split('\t')
        if len(parts) < 8:
            return None
        
        try:
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            
            # Парсим INFO поле в словарь
            info_dict = {}
            if parts[7] != '.':
                for info_item in parts[7].split(';'):
                    if '=' in info_item:
                        key, value = info_item.split('=', 1)
                        info_dict[key] = value
                    else:
                        info_dict[info_item] = True
            
            return VariantRecordVcf(chrom, pos, ref, alt, info_dict)
            
        except (ValueError, IndexError) as e:
            print(f"Ошибка парсинга варианта: {e}")
            return None

    def get_statistics(self) -> Dict:
        """
        Возвращает базовую статистику по VCF файлу.
        
        Returns:
            Dict: Словарь со статистикой:
                - total_variants: общее количество вариантов
                - samples_count: количество образцов
                - header_lines: количество строк в заголовке
        """
        # Пересчитываем количество вариантов для актуальности
        count = self._count_variants()
        return {
            'total_variants': count,
            'samples_count': len(self.samples),
            'header_lines': len(self.headers)
        }

    def _count_variants(self) -> int:
        """
        Подсчитывает общее количество вариантов в файле.
        
        Returns:
            int: Количество строк с вариантами (не-заголовочных строк)
            
        Примечание:
            Метод независимо подсчитывает варианты, временно перемещая
            позицию чтения файла и восстанавливая её после подсчета.
        """
        if not self.file:
            return 0
        
        current_pos = self.file.tell()
        self.file.seek(0)
        
        # Пропускаем заголовки
        line = self.file.readline()
        while line and line.startswith('#'):
            line = self.file.readline()
        
        # Считаем варианты (все не-пустые не-заголовочные строки)
        count = 0
        while line:
            line = line.strip()
            if line and not line.startswith('#'): 
                count += 1
            line = self.file.readline()
        
        # Восстанавливаем исходную позицию чтения
        self.file.seek(current_pos)
        return count

    def get_chromosomes(self) -> List[str]:
        """
        Возвращает список уникальных хромосом, представленных в файле.
        
        Returns:
            List[str]: Отсортированный список названий хромосом
            
        Примечание:
            Для больших файлов может быть затратной операцией,
            так требует полного прохода по файлу.
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
        
        # Собираем уникальные хромосомы
        while line:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split('\t')
                if parts:
                    chromosomes.add(parts[0])  # CHROM - первая колонка
            line = self.file.readline()
        
        # Восстанавливаем позицию чтения
        self.file.seek(current_pos)
        
        return sorted(list(chromosomes))

    def get_region_stats(self) -> pd.DataFrame:
        """
        Возвращает распределение вариантов по хромосомам.
        
        Returns:
            pd.DataFrame: DataFrame с колонками:
                - chromosome: название хромосомы
                - variant_count: количество вариантов на хромосоме
                
        Использование:
            stats_df = reader.get_region_stats()
            print(stats_df)  # Показывает таблицу с распределением
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
        
        # Восстанавливаем позицию чтения
        self.file.seek(current_pos)
        
        # Создаем DataFrame из собранных данных
        data = []
        for chrom, count in region_counts.items():
            data.append({'chromosome': chrom, 'variant_count': count})
        
        return pd.DataFrame(data)

    def get_variant_type_stats(self) -> pd.DataFrame:
        """
        Анализирует и возвращает статистику по типам геномных вариантов.
        
        Классификация вариантов:
            - SNV (Single Nucleotide Variant): одинарные нуклеотидные замены
            - Insertion: вставки (REF короче ALT)
            - Deletion: делеции (REF длиннее ALT)  
            - Complex: сложные варианты (равная длина, но разные последовательности)
            
        Returns:
            pd.DataFrame: DataFrame с колонками:
                - variant_type: тип варианта
                - count: количество вариантов этого типа
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
        
        # Анализируем типы вариантов по длине REF и ALT
        while line:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 5:
                    ref = parts[3]
                    alt = parts[4]
                    
                    # Определяем тип варианта по длине последовательностей
                    if len(ref) == 1 and len(alt) == 1:
                        type_counts['SNV'] += 1
                    elif len(ref) < len(alt):
                        type_counts['Insertion'] += 1
                    elif len(ref) > len(alt):
                        type_counts['Deletion'] += 1
                    else:
                        type_counts['Complex'] += 1
            
            line = self.file.readline()
        
        # Восстанавливаем позицию чтения
        self.file.seek(current_pos)
        
        # Создаем DataFrame только для ненулевых типов
        data = []
        for var_type, count in type_counts.items():
            if count > 0:
                data.append({'variant_type': var_type, 'count': count})
        
        return pd.DataFrame(data)

    def reset(self):
        """
        Сбрасывает состояние ридера к начальному.
        
        Действия:
            - Перемещает позицию чтения файла в начало
            - Перепарсивает заголовок файла
            - Сбрасывает внутренние счетчики и состояния
        """
        if self.file:
            self.file.seek(0)
            self._parse_header()

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Обеспечивает корректное завершение работы с файлом при выходе из контекста.
        
        Args:
            exc_type: Тип исключения (если было)
            exc_val: Значение исключения (если было) 
            exc_tb: Traceback исключения (если было)
            
        Действия:
            - Вызывает родительский __exit__ для закрытия файла
            - Сбрасывает все внутренние состояния ридера
        """
        super().__exit__(exc_type, exc_val, exc_tb)
        self._header_parsed = False
        self.headers = []
        self.samples = []
        self._variants_count = 0
