from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterator
from record import SequenceRecord

class SequenceReader(ABC):
    """Абстрактный базовый класс для чтения последовательностей."""
    
    def __init__(self, filepath: str | Path):
        self.filepath = Path(filepath)
    
    @abstractmethod
    def read(self) -> Iterator[SequenceRecord]:
        """Чтение последовательностей."""
        pass
    
    @abstractmethod
    def get_seq_count(self) -> int:
        """Получить количество последовательностей."""
        pass
    
    @abstractmethod
    def get_mean_seq_length(self) -> float:
        """Получить среднюю длину последовательностей."""
        pass
