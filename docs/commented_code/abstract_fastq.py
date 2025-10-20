from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterator
from record_fastq import SequenceRecordFastq


class SequenceReader(ABC):
    """
    Абстрактный базовый класс для чтения файлов с биологическими последовательностями.

    Этот класс определяет общий интерфейс для парсеров файлов различных форматов
    (FASTA, FASTQ, GenBank и т.д.). Конкретные реализации должны наследовать от этого класса
    и реализовать все абстрактные методы.

    Attributes:
        filepath (Path): Путь к файлу с последовательностями, преобразованный в объект Path.

    Example:
        class FastaReader(SequenceReader):
            def read(self) -> Iterator[SequenceRecord]:
                # Реализация чтения FASTA
                ...

        with FastaReader("sequences.fasta") as reader:
            for record in reader.read():
               print(record.header)
    """

    def __init__(self, filepath: str | Path) -> None:
        """
        Инициализирует парсер последовательностей.

        Аргументы:
            filepath (str | Path): Путь к входному файлу. Может быть строкой или объектом Path.
                Внутренне преобразуется в Path для единообразия работы с путями.

        Исключения:
            FileNotFoundError: Если указанный файл не существует.
            PermissionError: Если нет прав на чтение файла.
        """
        self.filepath = Path(filepath)
        # Проверка существования файла может быть добавлена здесь
        # if not self.filepath.exists():
        #     raise FileNotFoundError(f"Файл {self.filepath} не найден")

    @abstractmethod
    def read(self) -> Iterator[SequenceRecordFastq]:
        """
        Генератор для последовательного чтения записей из файла.

        Возвращает:
            Iterator[SequenceRecord]: Итератор объектов SequenceRecord. Каждый объект содержит:
                - header (str): Заголовок последовательности
                - sequence (str): Строка с биологической последовательностью (ДНК, РНК, белок)
                - quality (str, опционально): Строка качества для форматов с поддержкой качества

        Особенности:
            - Метод должен работать лениво (generator), чтобы эффективно обрабатывать большие файлы
            - Реализация должна корректно обрабатывать многострочные последовательности
            - Пустые файлы должны возвращать пустой итератор
            - Кодировка файла: рекомендуется использовать UTF-8 или ASCII

        Пример:
            >>> reader = FastaReader("test.fasta")
            >>> for record in reader.read():
            ...    print(f"{record.header}: {len(record.sequence)} нт")
        """
        pass

    @abstractmethod
    def get_seq_count(self) -> int:
        """
        Подсчитывает общее количество последовательностей в файле.

        Возвращает:
            int: Число уникальных последовательностей в файле

        Особенности:
            - Метод должен поддерживать повторный вызов без повторной загрузки файла
            - Реализация может требовать предварительного полного чтения файла
            - Для больших файлов рекомендуется использовать эффективные алгоритмы подсчета

        Пример:
            reader = FastqReader("reads.fastq")
            total_reads = reader.get_seq_count()
            print(f"Всего прочтений: {total_reads}")
        """
        pass

    @abstractmethod
    def get_mean_seq_length(self) -> float:
        """
        Вычисляет среднюю длину последовательностей в файле.

        Возвращает:
            float: Средняя длина последовательности с плавающей точкой

        Особенности:
            - Пустые последовательности должны игнорироваться в расчетах
            - Для пустого файла должно возвращаться 0.0
            - Реализация должна учитывать все символы последовательности (включая пробелы и переводы строк)

        Пример:
            reader = GenBankReader("sequences.gb")
            mean_len = reader.get_mean_seq_length()
            print(f"Средняя длина: {mean_len:.2f} нт")

        Примечание:
            Для точных расчетов рекомендуется использовать формулу:
            total_length / total_sequences
        """
        pass