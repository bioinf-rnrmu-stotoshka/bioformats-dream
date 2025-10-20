from abc import ABC, abstractmethod
from typing import Iterator, List, Dict, Tuple
from record_vcf import RecordVcf, SequenceRecordVcf
from pathlib import Path


class Reader(ABC):
    """
    Абстрактный базовый класс для всех ридеров геномных данных.
    
    Предоставляет общий интерфейс и базовую функциональность для чтения 
    файлов различных форматов. Реализует контекстный менеджер для безопасной
    работы с файлами.
    
    Attributes:
        filepath (Path): Путь к файлу для чтения
        file (FileIO | None): Файловый объект, открытый для чтения
    """
    
    def __init__(self, filepath: str | Path):
        """
        Инициализирует ридер с указанным путем к файлу.
        
        Args:
            filepath: Путь к файлу в виде строки или объекта Path
            
        Raises:
            FileNotFoundError: Если файл не существует
        """
        self.filepath = Path(filepath)
        self.file = None

    @abstractmethod
    def read(self) -> Iterator[RecordVcf]:
        """
        Абстрактный метод для чтения записей из файла.
        
        Returns:
            Iterator[Record]: Итератор, возвращающий объекты Record
            
        Raises:
            IOError: При ошибках чтения файла
            ValueError: При ошибках парсинга данных
        """
        pass

    def close(self):
        """
        Безопасно закрывает файл, если он открыт.
        
        Метод идемпотентен - многократный вызов не вызывает ошибок.
        """
        if self.file and not self.file.closed:
            self.file.close()
            self.file = None

    def __enter__(self):
        """
        Вход в контекстный менеджер.
        
        Returns:
            Reader: Экземпляр ридера с открытым файлом
            
        Raises:
            FileNotFoundError: Если файл не существует
            PermissionError: Если нет прав на чтение файла
        """
        self.file = open(self.filepath, "r", encoding="utf-8")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Выход из контекстного менеджера.
        
        Args:
            exc_type: Тип исключения (если было)
            exc_val: Значение исключения (если было)
            exc_tb: Traceback исключения (если было)
        """
        self.close()


class SequenceReader(Reader):
    """
    Абстрактный класс для ридеров последовательностей (FASTA, FASTQ).
    
    Предназначен для работы с файлами, содержащими биологические 
    последовательности (нуклеотидные или аминокислотные).
    
    Inherits:
        Reader: Базовый класс ридера
    """
    
    def __init__(self, filepath: str | Path):
        """
        Инициализирует ридер последовательностей.
        
        Args:
            filepath: Путь к файлу в виде строки или объекта Path
        """
        super().__init__(filepath)

    @abstractmethod
    def read(self) -> Iterator[SequenceRecordVcf]:
        """
        Абстрактный метод для чтения последовательностей из файла.
        
        Returns:
            Iterator[SequenceRecord]: Итератор объектов SequenceRecord
            
        Yields:
            SequenceRecord: Объекты, содержащие последовательности и метаданные
            
        Raises:
            ValueError: При нарушении формата файла
            IOError: При ошибках чтения файла
        """
        pass

    def __enter__(self):
        """
        Открывает файл при входе в блок with.
        
        Returns:
            SequenceReader: Экземпляр ридера с открытым файлом
        """
        self.file = open(self.filepath, "r", encoding="utf-8")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Автоматически закрывает файл при выходе из блока 'with'.
        
        Args:
            exc_type: Тип исключения (если было)
            exc_val: Значение исключения (если было) 
            exc_tb: Traceback исключения (если было)
        """
        self.close()

    def close(self):
        """
        Закрывает открытый файл.
        
        Рекомендуется использовать контекстный менеджер (with) 
        для автоматического управления ресурсами.
        """
        if self.file and not self.file.closed:
            self.file.close()
            self.file = None

    @abstractmethod
    def get_sequence(self, seq: str, id: str) -> SequenceRecordVcf:
        """
        Создаёт объект SequenceRecord из идентификатора и последовательности.
        
        В зависимости от формата (FASTA/FASTQ) может добавлять дополнительную
        информацию, такую как качество последовательности.
        
        Args:
            seq: Строка с последовательностью (нуклеотиды или аминокислоты)
            id: Идентификатор последовательности
            
        Returns:
            SequenceRecord: Объект с последовательностью и метаданными
            
        Raises:
            ValueError: Если последовательность содержит некорректные символы
        """
        pass

    @abstractmethod
    def validate_sequence(self, seq: str) -> bool:
        """
        Проверяет корректность последовательности.
        
        Валидация зависит от типа последовательности (ДНК, РНК, белок).
        Проверяет допустимые символы и возможные дополнительные правила.
        
        Args:
            seq: Последовательность для валидации
            
        Returns:
            bool: True если последовательность корректна, иначе False
        """
        pass


class GenomicDataReader(Reader):
    """
    Абстрактный класс для ридеров структурированных геномных данных.
    
    Предназначен для форматов SAM, VCF и других, содержащих аннотированные
    геномные варианты или выравнивания. Поддерживает работу с хромосомами,
    регионами и фильтрацию записей.
    
    Attributes:
        _header_parsed (bool): Флаг, указывающий был ли распарсен заголовок
        _chromosomes (List[str]): Список хромосом, обнаруженных в файле
        
    Inherits:
        Reader: Базовый класс ридера
    """
    
    def __init__(self, filepath: str | Path):
        """
        Инициализирует ридер геномных данных.
        
        Args:
            filepath: Путь к файлу в виде строки или объекта Path
        """
        super().__init__(filepath)
        self._header_parsed = False
        self._chromosomes = []

    def __enter__(self):
        """
        Открывает файл и парсит заголовок при входе в контекст.
        
        Returns:
            GenomicDataReader: Экземпляр ридера с распарсенным заголовком
            
        Raises:
            RuntimeError: При ошибках открытия файла или парсинга заголовка
            FileNotFoundError: Если файл не существует
        """
        try:
            self.file = open(self.filepath, "r", encoding="utf-8")
            self._parse_header()
            return self
        except Exception as e:
            # если ошибка — аккуратно закрываем файл, чтобы не остался висеть
            if self.file and not self.file.closed:
                self.file.close()
            raise RuntimeError(f"Ошибка при открытии или парсинге файла {self.filepath}: {e}")

    @abstractmethod
    def _parse_header(self):
        """
        Абстрактный метод для парсинга заголовка файла.
        
        Должен быть реализован в дочерних классах для извлечения информации
        о хромосомах, версиях формата, референсных последовательностях и т.д.
        
        Raises:
            ValueError: При нарушении формата заголовка
        """
        pass

    @abstractmethod
    def read(self) -> Iterator[RecordVcf]:
        """
        Чтение записей с пропуском заголовочных строк.
        
        Returns:
            Iterator[Record]: Итератор по объектам Record
            
        Yields:
            Record: Объекты, представляющие геномные записи (варианты, выравнивания)
            
        Raises:
            ValueError: При нарушении формата данных
        """
        pass

    @abstractmethod
    def get_chromosomes(self) -> List[str]:
        """
        Возвращает список хромосом, представленных в файле.
        
        Returns:
            List[str]: Список идентификаторов хромосом (например, ['chr1', 'chr2', ...])
        """
        pass

    @abstractmethod
    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        """
        Проверяет корректность координат.
        
        Args:
            chrom: Идентификатор хромосомы
            pos: Позиция на хромосоме (1-based)
            
        Returns:
            bool: True если координаты корректны, иначе False
        """
        pass

    @abstractmethod
    def get_records_in_region(self, chrom: str, start: int, end: int) -> Iterator[RecordVcf]:
        """
        Возвращает записи, попадающие в указанный геномный регион.
        
        Args:
            chrom: Идентификатор хромосомы
            start: Начальная позиция региона (включительно, 1-based)
            end: Конечная позиция региона (включительно, 1-based)
            
        Returns:
            Iterator[Record]: Итератор записей в указанном регионе
            
        Raises:
            ValueError: При некорректных координатах или неизвестной хромосоме
        """
        pass

    @abstractmethod
    def filter_records(self, **filters) -> Iterator[RecordVcf]:
        """
        Фильтрация записей по различным критериям.
        
        Конкретные критерии фильтрации зависят от формата файла и должны
        быть реализованы в дочерних классах.
        
        Args:
            **filters: Ключевые аргументы с критериями фильтрации.
                     Пример: quality=20, type='SNV', chromosome='chr1'
                     
        Returns:
            Iterator[Record]: Итератор отфильтрованных записей
            
        Raises:
            ValueError: При некорректных критериях фильтрации
        """
        pass

    def get_chromosome_length(self, chrom: str) -> int:
        """
        Получить длину хромосомы (если доступно из заголовка).
        
        Базовая реализация возвращает None. Должна быть переопределена 
        в дочерних классах, если информация о длинах хромосом доступна.
        
        Args:
            chrom: Идентификатор хромосомы
            
        Returns:
            int: Длина хромосомы в парах оснований или 0 если неизвестно
        """
        # Базовая реализация, может быть переопределена в дочерних классах
        pass

    def get_statistics(self) -> Dict[str, any]:
        """
        Возвращает базовую статистику по файлу.
        
        Должна быть расширена в дочерних классах для предоставления
        специфичной для формата статистики.
        
        Returns:
            Dict[str, any]: Словарь со статистикой, содержащий:
                - file_path (str): Путь к файлу
                - file_size (int): Размер файла в байтах
                - chromosomes (List[str]): Список хромосом
                - chromosome_count (int): Количество хромосом
        """
        return {
            "file_path": str(self.filepath),
            "file_size": self.filepath.stat().st_size if self.filepath.exists() else 0,
            "chromosomes": self.get_chromosomes(),
            "chromosome_count": len(self.get_chromosomes()),
        }

    def close(self):
        """
        Закрывает файл и сбрасывает состояние парсера.
        """
        super().close()
        self._header_parsed = False