from abc import ABC, abstractmethod
from typing import Iterator, List, Dict, Tuple
from record import Record, SequenceRecord
from pathlib import Path


class Reader(ABC):

    def __init__(self, filepath: str | Path):
        self.filepath = Path(filepath)
        self.file = None

    @abstractmethod
    def read(self) -> Iterator[Record]:
        pass

    def close(self):
        if self.file and not self.file.closed:
            self.file.close()
            self.file = None

    def __enter__(self):
        self.file = open(self.filepath, "r", encoding="utf-8")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class SequenceReader(Reader):

    def __init__(self, filepath: str | Path):
        super().__init__(filepath)

    @abstractmethod
    def read(self) -> Iterator[SequenceRecord]:
        pass

    def __enter__(self):
        self.file = open(self.filepath, "r", encoding="utf-8")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        if self.file and not self.file.closed:
            self.file.close()
            self.file = None

    @abstractmethod
    def get_sequence(self, seq: str, id: str) -> SequenceRecord:
        pass

    @abstractmethod
    def validate_sequence(self, seq: str) -> bool:
        pass


class GenomicDataReader(Reader):

    def __init__(self, filepath: str | Path):
        super().__init__(filepath)
        self._header_parsed = False
        self._chromosomes = []

    def __enter__(self):
        try:
            self.file = open(self.filepath, "r", encoding="utf-8")
            self._parse_header()
            return self
        except Exception as e:
            if self.file and not self.file.closed:
                self.file.close()
            raise RuntimeError(f"Ошибка при открытии или парсинге файла {self.filepath}: {e}")

    @abstractmethod
    def _parse_header(self):
        pass

    @abstractmethod
    def read(self) -> Iterator[Record]:
        pass

    @abstractmethod
    def get_chromosomes(self) -> List[str]:
        pass

    @abstractmethod
    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        pass

    @abstractmethod
    def get_records_in_region(self, chrom: str, start: int, end: int) -> Iterator[Record]:
        pass

    @abstractmethod
    def filter_records(self, **filters) -> Iterator[Record]:
        pass

    def get_chromosome_length(self, chrom: str) -> int:
        pass

    def get_statistics(self) -> Dict[str, any]:
        return {
            "file_path": str(self.filepath),
            "file_size": self.filepath.stat().st_size if self.filepath.exists() else 0,
            "chromosomes": self.get_chromosomes(),
            "chromosome_count": len(self.get_chromosomes()),
        }

    def close(self):
        """Закрывает файл"""
        super().close()
        self._header_parsed = False
