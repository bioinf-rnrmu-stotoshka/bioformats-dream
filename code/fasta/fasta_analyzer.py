from pathlib import Path
from typing import Iterator
from abstract import SequenceReader
from record import SequenceRecord

class FastaAnalyzer(SequenceReader):
    """
    Анализатор FASTA файлов.
    Реализует получение количества последовательностей и средней длины.
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
        """Читает последовательности из FASTA файла."""
        if not self.file:
            self.file = open(self.filepath, "r")

        current_header = None
        current_sequence = []
        
        for line in self.file:
            line = line.strip()
            
            if not line:
                continue
                
            if line.startswith(">"):
                if current_header is not None:
                    sequence = "".join(current_sequence)
                    record = SequenceRecord(id=current_header, sequence=sequence)
                    
                    # ОБНОВЛЯЕМ СТАТИСТИКУ
                    self._seq_count += 1
                    self._total_length += len(sequence)
                    
                    yield record
                
                current_header = line[1:].strip()
                current_sequence = []
            else:
                current_sequence.append(line)
        
        if current_header is not None:
            sequence = "".join(current_sequence)
            record = SequenceRecord(id=current_header, sequence=sequence)
            
            self._seq_count += 1
            self._total_length += len(sequence)
            
            yield record

    def get_seq_count(self) -> int:
        """ПОЛУЧЕНИЕ КОЛИЧЕСТВА ПОСЛЕДОВАТЕЛЬНОСТЕЙ - ТРЕБОВАНИЕ 1"""
        return self._seq_count

    def get_mean_seq_length(self) -> float:
        """ПОЛУЧЕНИЕ СРЕДНЕЙ ДЛИНЫ - ТРЕБОВАНИЕ 2"""
        if self._seq_count == 0:
            return 0.0
        return self._total_length / self._seq_count


def demo_fasta_analysis(filepath: str):
    """Демонстрационная функция для задания."""
    print("=== FASTA Analysis ===")
    print(f"File: {filepath}")
    
    with FastaAnalyzer(filepath) as analyzer:
        # Читаем все последовательности (обновляет статистику)
        sequences = list(analyzer.read())
        
        # ВЫВОДИМ ТО, ЧТО ТРЕБУЕТСЯ В ЗАДАНИИ:
        print(f"1. Количество последовательностей: {analyzer.get_seq_count()}")
        print(f"2. Средняя длина последовательностей: {analyzer.get_mean_seq_length():.2f}")

if __name__ == "__main__":
    """Запуск демо при прямом выполнении файла."""
    import sys
    if len(sys.argv) > 1:
        demo_fasta_analysis(sys.argv[1])
    else:
        print("Usage: python fastq_analyzer.py <fastq_file>")
        print("Using example file...")
        demo_fasta_analysis(input("Введите путь к FASTA файлу: "))
