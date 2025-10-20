from abc import ABC, abstractmethod
from typing import Iterator, List, Dict
from record import Record, SequenceRecord
from pathlib import Path


class Reader(ABC):
    """
    Абстрактный класс, который показывает логику работу ридеров
    """

    def __init__(self, filepath: str | Path):
        self.filepath = Path(filepath)
        self.file = None

    @abstractmethod
    def read(self) -> Iterator[Record]:
        """
        Абстрактный метод
        Должен возвращать итератор по объектам типа Record
        """
        pass

    def close(self):
        if self.file and not self.file.closed:
            self.file.close()
            self.file = None

    def __enter__(self):
        self.file = open(self.filepath, "r")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class SequenceReader(Reader):
    """
    Абстрактный класс, который показывает логику работы ридеров FASTA FASTQ
    """

    def __init__(self, filepath: str | Path):
        super().__init__(filepath)

    @abstractmethod
    def read(self) -> Iterator[SequenceRecord]:
        """
        Абстрактный метод
        Должен возвращать итератор объектов SequenceRecord
        """
        pass


    @abstractmethod
    def _get_sequence(self, seq_id: str, seq: str) -> SequenceRecord:
        """
        Создаёт объект SequenceRecord из идентификатора и последовательности.
        """
        pass

    @abstractmethod
    def _validate_sequence(self, seq: str) -> bool:
        """
        Проверяет корректность последовательности
        """
        pass


class GenomicDataReader(Reader):
    """
    Абстрактный класс для ридеров геномных данных
    """

    def __init__(self, filepath: str | Path):
        super().__init__(filepath)
        self._header_parsed = False

    def __enter__(self):
        """
        Открывает файл и парсит заголовок при входе в контекст
        """
        try:
            self.file = open(self.filepath, "r")
            self._parse_header()
            return self
        except Exception as e:
            if self.file and not self.file.closed:
                self.file.close()
            raise RuntimeError(f"Ошибка при открытии или парсинге файла {self.filepath}: {e}")

    @abstractmethod
    def _parse_header(self):
        """
        Абстрактный метод для парсинга заголовка файла
        Должен быть реализован в дочерних классах
        """
        pass

    @abstractmethod
    def read(self) -> Iterator[Record]:
        """
        Чтение записей с пропуском заголовочных строк
        """
        pass

    def close(self):
        """Закрывает файл"""
        super().close()
        self._header_parsed = False
