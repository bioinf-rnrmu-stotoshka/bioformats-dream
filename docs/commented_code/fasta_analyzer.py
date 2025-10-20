from pathlib import Path
from typing import Iterator
from abstract_fasta import SequenceReader
from record_fasta import SequenceRecordFasta

class FastaAnalyzer(SequenceReader):
    """
    Анализатор FASTA файлов для чтения и анализа биологических последовательностей.
    
    Реализует:
    1. Подсчет количества последовательностей в файле
    2. Расчет средней длины последовательностей
    3. Итеративное чтение записей с поддержкой контекстного менеджера
    
    Атрибуты:
        filepath (Path): Путь к анализируемому FASTA файлу
        file (file object): Открытый файловый объект
        _seq_count (int): Внутренний счетчик последовательностей
        _total_length (int): Суммарная длина всех последовательностей
    """

    def __init__(self, filepath: str | Path):
        """
        Инициализация анализатора FASTA файлов.
        
        Args:
            filepath: Путь к FASTA файлу в виде строки или объекта Path
        """
        super().__init__(filepath)
        self.file = None  # Файловый объект будет открыт в контекстном менеджере
        self._seq_count = 0  # Счетчик последовательностей
        self._total_length = 0  # Накопитель общей длины последовательностей

    def __enter__(self):
        """Реализация протокола контекстного менеджера для автоматического открытия файла."""
        self.file = open(self.filepath, "r")
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Гарантированное закрытие файла при выходе из контекста."""
        self.close()

    def close(self):
        """Явное закрытие файлового дескриптора."""
        if self.file:
            self.file.close()
            self.file = None

    def read(self) -> Iterator[SequenceRecordFasta]:
        """
        Основной метод чтения FASTA файла с подсчетом статистики.
        
        Генератор, который последовательно возвращает записи из FASTA файла.
        Одновременно обновляет внутреннюю статистику для последующих запросов.
        
        Логика обработки:
        - Пропускает пустые строки
        - Строки, начинающиеся с '>' - заголовки последовательностей
        - Последующие строки до следующего заголовка - данные последовательности
        
        Yields:
            SequenceRecord: Объект записи с идентификатором и последовательностью
            
        Обновляет:
            _seq_count: увеличивает на 1 для каждой найденной последовательности
            _total_length: добавляет длину каждой последовательности
        """
        # Автоматическое открытие файла если он не был открыт через контекстный менеджер
        if not self.file:
            self.file = open(self.filepath, "r")

        current_header = None  # Текущий обрабатываемый заголовок
        current_sequence = []  # Накопитель для строк последовательности
        
        for line in self.file:
            line = line.strip()
            
            if not line:  # Пропуск пустых строк
                continue
                
            if line.startswith(">"):  # Обнаружен новый заголовок
                # Если уже обрабатывали предыдущую последовательность - отдаем ее
                if current_header is not None:
                    sequence = "".join(current_sequence)
                    record = SequenceRecordFasta(id=current_header, sequence=sequence)
                    
                    # ОБНОВЛЕНИЕ СТАТИСТИКИ
                    self._seq_count += 1
                    self._total_length += len(sequence)
                    
                    yield record
                
                # Начинаем новую последовательность
                current_header = line[1:].strip()  # Убираем '>' и лишние пробелы
                current_sequence = []
            else:
                # Накопление данных последовательности
                current_sequence.append(line)
        
        # Обработка последней последовательности в файле
        if current_header is not None:
            sequence = "".join(current_sequence)
            record = SequenceRecordFasta(id=current_header, sequence=sequence)
            
            self._seq_count += 1
            self._total_length += len(sequence)
            
            yield record

    def get_seq_count(self) -> int:
        """
        ПОЛУЧЕНИЕ КОЛИЧЕСТВА ПОСЛЕДОВАТЕЛЬНОСТЕЙ - ТРЕБОВАНИЕ 1
        
        Returns:
            int: Общее количество последовательностей в файле
            
        Note:
            Для корректного подсчета требуется предварительный вызов read()
        """
        return self._seq_count

    def get_mean_seq_length(self) -> float:
        """
        ПОЛУЧЕНИЕ СРЕДНЕЙ ДЛИНЫ - ТРЕБОВАНИЕ 2
        
        Returns:
            float: Средняя длина последовательности с плавающей точкой
            
        Note:
            Возвращает 0.0 если файл не содержит последовательностей
        """
        if self._seq_count == 0:
            return 0.0
        return self._total_length / self._seq_count


def demo_fasta_analysis(filepath: str):
    """
    Демонстрационная функция анализа FASTA файла.
    
    Показывает использование основных возможностей класса FastaAnalyzer:
    1. Чтение всех последовательностей
    2. Подсчет общего количества
    3. Расчет средней длины
    
    Args:
        filepath: Путь к FASTA файлу для анализа
    """
    print("=== FASTA Analysis ===")
    print(f"File: {filepath}")
    
    # Использование контекстного менеджера для гарантированного закрытия файла
    with FastaAnalyzer(filepath) as analyzer:
        # Чтение всех последовательностей (одновременно обновляет статистику)
        sequences = list(analyzer.read())
        
        # ВЫВОД РЕЗУЛЬТАТОВ СОГЛАСНО ТРЕБОВАНИЯМ ЗАДАНИЯ:
        print(f"1. Количество последовательностей: {analyzer.get_seq_count()}")
        print(f"2. Средняя длина последовательностей: {analyzer.get_mean_seq_length():.2f}")

if __name__ == "__main__":
    """
    Точка входа при прямом запуске скрипта.
    
    Обрабатывает аргументы командной строки:
    - Если передан путь к файлу - анализирует указанный файл
    - Иначе запрашивает путь у пользователя
    """
    import sys
    if len(sys.argv) > 1:
        demo_fasta_analysis(sys.argv[1])
    else:
        print("Usage: python fastq_analyzer.py <fastq_file>")
        print("Using example file...")
        demo_fasta_analysis(input("Введите путь к FASTA файлу: "))