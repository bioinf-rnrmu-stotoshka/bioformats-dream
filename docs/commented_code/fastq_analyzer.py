from pathlib import Path
from typing import Iterator, List, Dict
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict
from abstract_fastq import SequenceReader
from record_fastq import SequenceRecordFastq


class FastqAnalyzer(SequenceReader):
    """
    Анализатор FASTQ файлов для биоинформатического анализа.
    
    Наследуется от абстрактного класса SequenceReader и предоставляет функциональность
    для чтения и анализа данных из файлов в формате FASTQ.
    
    FASTQ формат содержит последовательности ДНК/РНК с оценками качества:
    - Строка 1: Заголовок (начинается с @)
    - Строка 2: Нуклеотидная последовательность
    - Строка 3: Разделитель (обычно +)
    - Строка 4: Строка качества в формате Phred
    
    Атрибуты:
        filepath (Path): Путь к анализируемому FASTQ файлу
        _seq_count (int): Счетчик прочитанных последовательностей
        _total_length (int): Общая длина всех последовательностей
    """

    def __init__(self, filepath: str | Path):
        """
        Инициализирует анализатор FASTQ файлов.
        
        Args:
            filepath: Путь к FASTQ файлу в виде строки или объекта Path
        """
        super().__init__(filepath)
        self._seq_count = 0
        self._total_length = 0

    def __enter__(self):
        """Реализация контекстного менеджера для использования с 'with'."""
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Заглушка для контекстного менеджера (очистка ресурсов не требуется)."""
        pass

    def read(self) -> Iterator[SequenceRecordFastq]:
        """
        Читает последовательности из FASTQ файла и возвращает итератор.
        
        Генератор, который читает файл построчно, группирует строки по 4 (полная запись FASTQ)
        и создает объекты SequenceRecord для каждой последовательности.
        
        Yields:
            SequenceRecord: Объект с информацией о последовательности (ID и последовательность)
            
        Пример FASTQ записи:
            @seq1
            ATGCGATCGATCG
            +
            IIIIIIIIIIIII
        """
        with open(self.filepath, "r") as file:
            lines = []
            for line in file:
                lines.append(line.strip())
                if len(lines) == 4:
                    header = lines[0][1:]  # Убираем символ @ из заголовка
                    sequence = lines[1]
                    record = SequenceRecordFastq(id=header, sequence=sequence)
                    
                    # Обновляем статистику
                    self._seq_count += 1
                    self._total_length += len(sequence)
                    
                    yield record
                    lines = []  # Сбрасываем буфер для следующей записи

    def get_seq_count(self) -> int:
        """
        Возвращает общее количество прочитанных последовательностей.
        
        Returns:
            int: Количество последовательностей в файле
        """
        return self._seq_count

    def get_mean_seq_length(self) -> float:
        """
        Вычисляет среднюю длину последовательностей.
        
        Returns:
            float: Средняя длина последовательности в нуклеотидах
                   Возвращает 0.0 если последовательностей нет
        """
        if self._seq_count == 0:
            return 0.0
        return self._total_length / self._seq_count

    def phred_to_quality(self, phred_char: str, offset: int = 33) -> int:
        """
        Конвертирует символ Phred в числовое значение качества.
        
        Phred качество кодируется одним символом ASCII. По умолчанию используется
        кодировка Sanger (offset=33), где качество = ASCII код - 33.
        
        Args:
            phred_char: Один символ из строки качества
            offset: Смещение для декодирования (33 для Sanger, 64 для Illumina 1.3+)
            
        Returns:
            int: Численное значение качества от 0 до 40+
            
        Пример:
            'I' (ASCII 73) -> 73 - 33 = 40
        """
        return ord(phred_char) - offset

    def get_sequences_with_quality(self) -> Iterator[tuple]:
        """
        Генератор для чтения последовательностей вместе с информацией о качестве.
        
        Читает полные записи FASTQ и возвращает кортежи с заголовком, 
        последовательностью и строкой качества.
        
        Yields:
            tuple: (header, sequence, quality_string) для каждой записи FASTQ
            
        Используется для анализа качества последовательностей.
        """
        with open(self.filepath, "r") as file:
            lines = []
            for line in file:
                lines.append(line.strip())
                if len(lines) == 4:
                    header = lines[0][1:]  # Без символа @
                    sequence = lines[1]
                    quality = lines[3]  # Строка качества
                    yield header, sequence, quality
                    lines = []  # Сброс для следующей записи

    def per_base_sequence_quality(self, output_file: str = "fastq_quality_plot.png"):
        """
        Строит график качества оснований по позициям в ридах.
        
        Анализирует качество каждого основания в каждой позиции рида и
        визуализирует среднее качество по всем последовательностям.
        
        Args:
            output_file: Имя файла для сохранения графика
            
        Returns:
            str: Путь к сохраненному файлу с графиком
        """
        # Словарь для хранения качества по позициям: {позиция: [качества]}
        quality_by_position = defaultdict(list)
        
        # Сбор данных о качестве для каждой позиции
        for _, sequence, quality_str in self.get_sequences_with_quality():
            for pos, qual_char in enumerate(quality_str):
                quality_score = self.phred_to_quality(qual_char)
                quality_by_position[pos].append(quality_score)
        
        # Проверка наличия данных
        if not quality_by_position:
            print("Нет данных для графика качества")
            return output_file
            
        # Подготовка данных для графика
        positions = sorted(quality_by_position.keys())
        mean_qualities = [np.mean(quality_by_position[pos]) for pos in positions]
        
        # Построение графика
        plt.figure(figsize=(12, 6))
        plt.plot(positions, mean_qualities, linewidth=2)
        plt.xlabel('Position in read (bp)')
        plt.ylabel('Quality score')
        plt.title('Per Base Sequence Quality')
        plt.grid(True, alpha=0.3)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"График качества создан: {output_file}")
        return output_file

    def per_base_sequence_content(self, output_file: str = "fastq_content_plot.png"):
        """
        Анализирует и визуализирует содержание нуклеотидов по позициям.
        
        Вычисляет процентное содержание A, T, G, C оснований для каждой позиции
        во всех последовательностях и строит линейный график.
        
        Args:
            output_file: Имя файла для сохранения графика
            
        Returns:
            str: Путь к сохраненному файлу с графиком
        """
        base_content = {'A': [], 'T': [], 'G': [], 'C': []}
        max_length = 0
        
        # Первый проход: сбор всех последовательностей и определение максимальной длины
        sequences_data = []
        for _, sequence, _ in self.get_sequences_with_quality():
            sequences_data.append(sequence)
            max_length = max(max_length, len(sequence))
        
        if max_length == 0:
            print("Нет данных для графика содержания")
            return output_file
        
        # Инициализация счетчиков для каждой позиции
        for base in base_content:
            base_content[base] = [0] * max_length
        total_bases = [0] * max_length  # Общее количество оснований на каждой позиции
        
        # Второй проход: подсчет оснований по позициям
        for sequence in sequences_data:
            for pos, base in enumerate(sequence):
                if base in base_content and pos < max_length:
                    base_content[base][pos] += 1
                    total_bases[pos] += 1
        
        # Расчет процентного содержания для каждого основания
        percentages = {}
        for base in base_content:
            percentages[base] = [100 * count / total if total > 0 else 0 
                               for count, total in zip(base_content[base], total_bases)]
        
        # Построение графика
        plt.figure(figsize=(12, 6))
        positions = list(range(max_length))
        
        # Цвета для разных оснований
        for base, color in zip(['A', 'T', 'G', 'C'], ['green', 'red', 'black', 'blue']):
            plt.plot(positions, percentages[base], label=base, color=color, linewidth=2)
        
        plt.xlabel('Position in read (bp)')
        plt.ylabel('Percentage (%)')
        plt.title('Per Base Sequence Content')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"График содержания создан: {output_file}")
        return output_file

    def sequence_length_distribution(self, output_file: str = "fastq_length_plot.png"):
        """
        Анализирует и визуализирует распределение длин последовательностей.
        
        Строит гистограмму, показывающую частоту встречаемости ридов разной длины.
        
        Args:
            output_file: Имя файла для сохранения гистограммы
            
        Returns:
            str: Путь к сохраненному файлу с графиком
        """
        # Сбор данных о длинах последовательностей
        lengths = []
        for record in self.read():
            lengths.append(len(record.sequence))
        
        if not lengths:
            print("Нет данных для графика распределения длин")
            return output_file
            
        # Построение гистограммы
        plt.figure(figsize=(10, 6))
        plt.hist(lengths, bins=20, alpha=0.7, edgecolor='black')
        plt.xlabel('Sequence Length (bp)')
        plt.ylabel('Frequency')
        plt.title('Sequence Length Distribution')
        plt.grid(True, alpha=0.3)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"График распределения длин создан: {output_file}")
        return output_file


def demo_fastq_analysis(filepath: str):
    """
    Демонстрационная функция для комплексного анализа FASTQ файла.
    
    Выполняет полный анализ FASTQ файла включая:
    - Проверку читаемости файла
    - Базовую статистику
    - Создание диагностических графиков
    
    Args:
        filepath: Путь к FASTQ файлу для анализа
    """
    print("=== FASTQ Analysis ===")
    print(f"File: {filepath}")
    
    # Предварительная проверка файла
    try:
        with open(filepath, 'r') as f:
            first_lines = [next(f).strip() for _ in range(4)]
        print("Файл читается:")
        for i, line in enumerate(first_lines):
            print(f"  {i+1}: {line[:50]}...")  # Показываем первые 50 символов каждой строки
    except Exception as e:
        print(f"❌ Ошибка чтения файла: {e}")
        return
    
    # Создание анализатора и выполнение анализа
    analyzer = FastqAnalyzer(filepath)
    
    # Базовая статистика
    sequences = list(analyzer.read())
    print(f"1. Количество последовательностей: {analyzer.get_seq_count()}")
    print(f"2. Средняя длина последовательностей: {analyzer.get_mean_seq_length():.2f}")
    
    if analyzer.get_seq_count() == 0:
        print("Нет последовательностей для анализа!")
        return
    
    # Создание диагностических графиков
    print("3. Создание графиков...")
    quality_plot = analyzer.per_base_sequence_quality()
    content_plot = analyzer.per_base_sequence_content()
    length_plot = analyzer.sequence_length_distribution()
    
    print(f"4. График качества: {quality_plot}")
    print(f"5. График содержания: {content_plot}")
    print(f"6. Распределение длин: {length_plot}")


if __name__ == "__main__":
    """
    Точка входа для запуска анализа из командной строки.
    
    Использование:
        python fastq_analyzer.py <путь_к_fastq_файлу>
        
    Если файл не указан, запрашивает путь у пользователя.
    """
    import sys
    if len(sys.argv) > 1:
        demo_fastq_analysis(sys.argv[1])
    else:
        print("Usage: python fastq_analyzer.py <fastq_file>")
        print("Using example file...")
        demo_fastq_analysis(input("Введите путь к FASTQ файлу: "))