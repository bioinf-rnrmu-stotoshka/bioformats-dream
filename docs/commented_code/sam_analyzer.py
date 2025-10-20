import pandas as pd
from pathlib import Path
from typing import Iterator, Dict, List, Any
import re
from abstract_sam import GenomicDataReader
from record_sam import AlignmentRecordSam

class SamReader(GenomicDataReader):
    """
    Парсер SAM-файлов для работы с данными геномных выравниваний.
    
    Наследуется от GenomicDataReader и предоставляет методы для чтения,
    анализа и фильтрации данных выравнивания из SAM-файлов.
    
    Attributes:
        filepath (Path): Путь к SAM-файлу
        file (IO): Файловый объект для чтения данных
        header (Dict[str, List[str]]): Заголовок SAM-файла, разобранный по тегам
    """
    
    def __init__(self, filepath: str | Path):
        """
        Инициализация парсера SAM-файла.
        
        Args:
            filepath: Путь к SAM-файлу (строка или Path-объект)
        """
        super().__init__(filepath)
        self.header: Dict[str, List[str]] = {}

    def _parse_header(self):
        """
        Приватный метод для парсинга заголовка SAM-файла.
        
        Читает все строки, начинающиеся с '@', и сохраняет их в атрибуте header.
        После парсинга возвращает позицию чтения файла в начало.
        """
        if not self.file:
            self.file = open(self.filepath, "r", encoding="utf-8", errors="replace")

        pos = self.file.tell()
        for line in self.file:
            if not line.startswith("@"):
                break
            self._parse_header_line(line)
        self.file.seek(pos)

    def _parse_header_line(self, line: str):
        """
        Парсит отдельную строку заголовка SAM-файла.
        
        Args:
            line: Строка заголовка, начинающаяся с '@'
        """
        parts = line.strip().split("\t")
        tag = parts[0]
        self.header.setdefault(tag, []).append("\t".join(parts[1:]))

    def get_header(self) -> Dict[str, List[str]]:
        """
        Возвращает заголовок SAM-файла.
        
        Returns:
            Словарь с ключами - тегами заголовка (например, '@SQ', '@RG')
            и значениями - списками соответствующих данных
        """
        return self.header

    def read(self) -> Iterator[AlignmentRecordSam]:
        """
        Основной метод для чтения записей выравнивания из SAM-файла.
        
        Пропускает строки заголовка и записи с некорректными данными.
        Генерирует объекты AlignmentRecord для каждой валидной записи.
        
        Yields:
            AlignmentRecord: Объект с данными выравнивания
            
        Notes:
            - Пропускает записи где RNAME или CIGAR равны '*'
            - Автоматически вычисляет конечную позицию выравнивания
        """
        if not self.file:
            self.file = open(self.filepath, "r", encoding="utf-8", errors="replace")

        for line in self.file:
            if line.startswith("@"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 11:
                continue

            if fields[2] == "*" or fields[5] == "*":
                continue

            qname, flag, rname, pos, mapq, cigar = (
                fields[0],
                fields[1],
                fields[2],
                fields[3],
                fields[4],
                fields[5],
            )

            aligned_len = self._calc_aligned_length(cigar)
            end_pos = int(pos) + aligned_len - 1 if aligned_len > 0 else int(pos)

            rec = AlignmentRecordSam(
                id=qname, chrom=rname, start=int(pos), cigar=cigar, mapq=int(mapq)
            )
            rec.flag = int(flag)
            rec.end = end_pos
            yield rec

    @staticmethod
    def _calc_aligned_length(cigar: str) -> int:
        """
        Вычисляет длину выравнивания на основе CIGAR-строки.
        
        Args:
            cigar: CIGAR-строка в стандартном формате
            
        Returns:
            Суммарная длина выравнивания с учетом операций:
            M (match), D (deletion), N (skip), = (match), X (mismatch)
            
        Example:
            "50M10D20N" -> 80
        """
        if not cigar or cigar == "*":
            return 0

        total = 0
        for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
            if op in ("M", "D", "N", "=", "X"):
                total += int(length)
        return total

    def get_chromosomes(self) -> List[str]:
        """
        Возвращает список всех хромосом, присутствующих в выравниваниях.
        
        Returns:
            Отсортированный список уникальных названий хромосом
        """
        chromosomes = set()
        # Сохраняем текущую позицию в файле
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for rec in self.read():
            if rec.chrom != "*" and rec.chrom:  # проверяем что не пустая и не *
                chromosomes.add(rec.chrom)

        # Восстанавливаем позицию в файле
        if self.file:
            self.file.seek(current_pos)
        return sorted(chromosomes)

    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        """
        Проверяет валидность геномной координаты.
        
        Args:
            chrom: Название хромосомы
            pos: Позиция на хромосоме
            
        Returns:
            True если хромосома существует и позиция положительна
        """
        return chrom in self.get_chromosomes() and pos > 0

    def calculate_coverage(self, chrom: str) -> Dict[int, int]:
        """
        Вычисляет покрытие для указанной хромосомы.
        
        Args:
            chrom: Название хромосомы для анализа
            
        Returns:
            Словарь где ключи - позиции в хромосоме,
            значения - количество прочитанных оснований в этой позиции
        """
        coverage = {}
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for rec in self.read():
            if rec.chrom == chrom and rec.cigar != "*":
                aligned_len = self._calc_aligned_length(rec.cigar)
                for i in range(rec.start, rec.start + aligned_len):
                    coverage[i] = coverage.get(i, 0) + 1

        if self.file:
            self.file.seek(current_pos)
        return coverage

    def filter_alignments(self, flag: int) -> Iterator[AlignmentRecordSam]:
        """
        Фильтрует выравнивания по флагу SAM.
        
        Args:
            flag: SAM-флаг для фильтрации (битовая маска)
            
        Yields:
            Записи выравнивания, удовлетворяющие критерию флага
        """
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for rec in self.read():
            if hasattr(rec, "flag") and rec.flag & flag:
                yield rec

        if self.file:
            self.file.seek(current_pos)

    def count_alignments(self) -> int:
        """
        Подсчитывает общее количество выравниваний в файле.
        
        Returns:
            Количество валидных записей выравнивания
        """
        cnt = 0
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for _ in self.read():
            cnt += 1

        if self.file:
            self.file.seek(current_pos)
        return cnt

    def stats_by_chromosome(self) -> pd.DataFrame:
        """
        Статистика по количеству выравниваний на хромосому.
        
        Returns:
            DataFrame с колонками:
            - chrom: название хромосомы
            - count: количество выравниваний
        """
        from collections import Counter

        cnt = Counter()
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for rec in self.read():
            cnt[rec.chrom] += 1

        if self.file:
            self.file.seek(current_pos)

        df = pd.DataFrame(list(cnt.items()), columns=["chrom", "count"])
        return df

    def filter_by_region(
        self, chrom: str, start: int, end: int
    ) -> Iterator[AlignmentRecordSam]:
        """
        Фильтрует выравнивания по геномному региону.
        
        Args:
            chrom: Хромосома
            start: Начальная позиция региона (1-based)
            end: Конечная позиция региона (1-based)
            
        Yields:
            Записи выравнивания, пересекающие указанный регион
        """
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for rec in self.read():
            if rec.chrom == chrom and rec.end >= start and rec.start <= end:
                yield rec

        if self.file:
            self.file.seek(current_pos)

    def get_records_in_region(
        self, chrom: str, start: int, end: int
    ) -> List[AlignmentRecordSam]:
        """
        Получает все записи выравнивания в указанном регионе.
        
        Args:
            chrom: Хромосома
            start: Начальная позиция региона
            end: Конечная позиция региона
            
        Returns:
            Список объектов AlignmentRecord в регионе
        """
        return list(self.filter_by_region(chrom, start, end))

    def filter_records(self, **filters) -> List[AlignmentRecordSam]:
        """
        Универсальный метод фильтрации записей по различным критериям.
        
        Args:
            **filters: Ключевые аргументы для фильтрации:
                - chrom: фильтр по хромосоме
                - min_mapq: минимальное качество маппинга
                - flag: SAM-флаг для фильтрации
                
        Returns:
            Список отфильтрованных записей AlignmentRecord
        """
        results = []
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)

        for rec in self.read():
            match = True

            if "chrom" in filters and rec.chrom != filters["chrom"]:
                match = False
            if "min_mapq" in filters and rec.mapq < filters["min_mapq"]:
                match = False
            if (
                "flag" in filters
                and hasattr(rec, "flag")
                and not (rec.flag & filters["flag"])
            ):
                match = False

            if match:
                results.append(rec)

        if self.file:
            self.file.seek(current_pos)
        return results