from pathlib import Path
from typing import Iterator, List, Dict
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict
from .abstract import SequenceReader
from .record import SequenceRecord


class FastqAnalyzer(SequenceReader):
    """
    Анализатор FASTQ файлов.
    
    Поддерживает:
    - Подсчет количества последовательностей
    - Расчет средней длины
    - Графики качества: per base sequence quality, per base sequence content, sequence length distribution
    """

    def __init__(self, filepath: str | Path):
        super().__init__(filepath)
        self.file = None
        self._seq_count = 0
        self._total_length = 0

    def __enter__(self):
        self.file = open(self.filepath, "r")
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        if self.file:
            self.file.close()
            self.file = None

    def read(self) -> Iterator[SequenceRecord]:
        """Читает последовательности из FASTQ файла."""
        if not self.file:
            self.file = open(self.filepath, "r")

        lines = []
        for line in self.file:
            lines.append(line.strip())
            if len(lines) == 4:  # FASTQ блок: header, sequence, +, quality
                header = lines[0][1:]  # Убираем @
                sequence = lines[1]
                quality = lines[3]
                
                record = SequenceRecord(id=header, sequence=sequence)
                
                # Обновляем статистику
                self._seq_count += 1
                self._total_length += len(sequence)
                
                yield record
                lines = []

    def get_seq_count(self) -> int:
        """Получить количество последовательностей."""
        return self._seq_count

    def get_mean_seq_length(self) -> float:
        """Получить среднюю длину последовательностей."""
        if self._seq_count == 0:
            return 0.0
        return self._total_length / self._seq_count

    # ДОПОЛНИТЕЛЬНЫЕ МЕТОДЫ ДЛЯ FASTQ

    def phred_to_quality(self, phred_char: str, offset: int = 33) -> int:
        """Конвертировать символ Phred в числовое значение качества."""
        return ord(phred_char) - offset

    def get_sequences_with_quality(self) -> Iterator[tuple]:
        """Генератор для чтения последовательностей с качеством."""
        if not self.file:
            self.file = open(self.filepath, "r")

        lines = []
        for line in self.file:
            lines.append(line.strip())
            if len(lines) == 4:
                header = lines[0][1:]
                sequence = lines[1]
                quality = lines[3]
                yield header, sequence, quality
                lines = []

    def per_base_sequence_quality(self, output_file: str = "fastq_quality_plot.png"):
        """Построить график качества по позициям."""
        quality_by_position = defaultdict(list)
        
        for _, sequence, quality_str in self.get_sequences_with_quality():
            for pos, qual_char in enumerate(quality_str):
                quality_score = self.phred_to_quality(qual_char)
                quality_by_position[pos].append(quality_score)
        
        # Расчет статистики
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
        
        return output_file

    def per_base_sequence_content(self, output_file: str = "fastq_content_plot.png"):
        """Построить график содержания нуклеотидов по позициям."""
        base_content = {'A': [], 'T': [], 'G': [], 'C': []}
        max_length = 0
        
        # Первый проход: найти максимальную длину
        for _, sequence, _ in self.get_sequences_with_quality():
            max_length = max(max_length, len(sequence))
        
        # Инициализация списков
        for base in base_content:
            base_content[base] = [0] * max_length
        total_bases = [0] * max_length
        
        # Второй проход: подсчет оснований
        for _, sequence, _ in self.get_sequences_with_quality():
            for pos, base in enumerate(sequence):
                if base in base_content and pos < max_length:
                    base_content[base][pos] += 1
                    total_bases[pos] += 1
        
        # Расчет процентов
        percentages = {}
        for base in base_content:
            percentages[base] = [100 * count / total if total > 0 else 0 
                               for count, total in zip(base_content[base], total_bases)]
        
        # Построение графика
        plt.figure(figsize=(12, 6))
        positions = list(range(max_length))
        
        for base, color in zip(['A', 'T', 'G', 'C'], ['green', 'red', 'black', 'blue']):
            plt.plot(positions, percentages[base], label=base, color=color, linewidth=2)
        
        plt.xlabel('Position in read (bp)')
        plt.ylabel('Percentage (%)')
        plt.title('Per Base Sequence Content')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        return output_file

    def sequence_length_distribution(self, output_file: str = "fastq_length_plot.png"):
        """Построить распределение длин последовательностей."""
        lengths = []
        for record in self.read():
            lengths.append(len(record.sequence))
        
        plt.figure(figsize=(10, 6))
        plt.hist(lengths, bins=50, alpha=0.7, edgecolor='black')
        plt.xlabel('Sequence Length (bp)')
        plt.ylabel('Frequency')
        plt.title('Sequence Length Distribution')
        plt.grid(True, alpha=0.3)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        return output_file


def demo_fastq_analysis(filepath: str):
    """Демонстрационная функция для FASTQ анализа."""
    print("=== FASTQ Analysis ===")
    print(f"File: {filepath}")
    
    with FastqAnalyzer(filepath) as analyzer:
        # Базовая статистика
        sequences = list(analyzer.read())
        
        print(f"1. Количество последовательностей: {analyzer.get_seq_count()}")
        print(f"2. Средняя длина последовательностей: {analyzer.get_mean_seq_length():.2f}")
        
        # Графики качества
        quality_plot = analyzer.per_base_sequence_quality()
        content_plot = analyzer.per_base_sequence_content()
        length_plot = analyzer.sequence_length_distribution()
        
        print(f"3. График качества сохранен как: {quality_plot}")
        print(f"4. График содержания нуклеотидов: {content_plot}")
        print(f"5. Распределение длин: {length_plot}")


if __name__ == "__main__":
    """Запуск демо при прямом выполнении файла."""
    import sys
    if len(sys.argv) > 1:
        demo_fastq_analysis(sys.argv[1])
    else:
        print("Usage: python fastq_analyzer.py <fastq_file>")
        print("Using example file...")
        demo_fastq_analysis("examples/example.fastq")
