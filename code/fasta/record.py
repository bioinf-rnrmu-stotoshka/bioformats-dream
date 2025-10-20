from dataclasses import dataclass

@dataclass
class SequenceRecord:
    """Класс для хранения информации о последовательности."""
    id: str
    sequence: str
    
    def __len__(self):
        return len(self.sequence)
