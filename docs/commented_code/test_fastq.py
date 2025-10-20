import pytest
import tempfile
import os
import sys

# Добавляем родительскую директорию в путь для импорта модулей
# Это позволяет импортировать модули из соседних пакетов
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fastq_analyzer import FastqAnalyzer


class TestFastqAnalyzer:
    """
    Набор тестов для класса FastqAnalyzer.
    
    Этот тестовый класс проверяет функциональность анализа FASTQ файлов,
    включая подсчет последовательностей, расчет средней длины и конвертацию 
    качества Phred.
    
    Attributes:
        Нет атрибутов экземпляра. Все тестовые данные создаются временно.
    """
    
    def create_test_fastq(self, content: str) -> str:
        """
        Создает временный FASTQ файл с заданным содержимым.
        
        Args:
            content (str): Содержимое для записи в FASTQ файл. Должно соответствовать
                         стандартному формату FASTQ с четырьмя строками на запись.
        
        Returns:
            str: Путь к созданному временному файлу.
            
        Note:
            Файл создается с суффиксом .fastq и не удаляется автоматически.
            Удаление файла должно выполняться в вызывающем коде.
        """
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fastq') as f:
            f.write(content)
            return f.name
    
    def test_fastq_count_sequences(self):
        """
        Тестирует корректность подсчета последовательностей в FASTQ файле.
        
        Создает тестовый FASTQ файл с двумя последовательностями и проверяет,
        что анализатор правильно определяет количество последовательностей.
        
        Test Case:
            - Входные данные: файл с двумя последовательностями (read1 и read2)
            - Ожидаемый результат: get_seq_count() возвращает 2
            
        Raises:
            AssertionError: если количество последовательностей не равно 2
        """
        # Тестовые данные в формате FASTQ
        fastq_content = """@read1
ATCG
+
IIII
@read2
GGCC
+
IIII
"""
        # Создаем временный файл
        test_file = self.create_test_fastq(fastq_content)
        
        try:
            # Используем контекстный менеджер для гарантированного закрытия файла
            with FastqAnalyzer(test_file) as analyzer:
                # Читаем все последовательности для инициализации подсчета
                sequences = list(analyzer.read())
                # Проверяем корректность подсчета
                assert analyzer.get_seq_count() == 2
        finally:
            # Всегда удаляем временный файл
            os.unlink(test_file)
    
    def test_fastq_mean_length(self):
        """
        Тестирует расчет средней длины последовательностей в FASTQ файле.
        
        Создает тестовый FASTQ файл с последовательностями разной длины
        и проверяет корректность расчета средней длины.
        
        Test Case:
            - Входные данные: две последовательности длиной 4 и 6 нуклеотидов
            - Ожидаемый результат: (4 + 6) / 2 = 5.0
            
        Raises:
            AssertionError: если средняя длина не соответствует расчетной
        """
        fastq_content = """@read1
ATCG
+
IIII
@read2
GGCCAA
+
IIIIII
"""
        test_file = self.create_test_fastq(fastq_content)
        
        try:
            with FastqAnalyzer(test_file) as analyzer:
                sequences = list(analyzer.read())
                # Проверяем расчет средней длины
                expected_mean = (4 + 6) / 2
                assert analyzer.get_mean_seq_length() == expected_mean
        finally:
            os.unlink(test_file)
    
    def test_phred_quality_conversion(self):
        """
        Тестирует конвертацию символов Phred в числовые значения качества.
        
        Проверяет корректность преобразования ASCII символов в соответствующие
        значения качества с использованием offset 33 (стандарт Sanger).
        
        Test Case:
            - Символ 'I' (ASCII 73) с offset 33: 73 - 33 = 40
            - Ожидаемый результат: 40
            
        Raises:
            AssertionError: если преобразование возвращает неверное значение
        """
        # Создаем экземпляр анализатора с заглушкой имени файла
        analyzer = FastqAnalyzer("dummy.fastq")
        # Проверяем преобразование Phred качества
        # Символ 'I' имеет ASCII код 73, при offset 33: 73-33=40
        assert analyzer.phred_to_quality('I') == 40


def test_demo_fastq_function():
    """
    Интеграционный тест для демонстрационной функции анализа FASTQ.
    
    Проверяет что демонстрационная функция может быть импортирована и
    выполняется без ошибок на корректном FASTQ файле.
    
    Test Case:
        - Создается валидный FASTQ файл с одной последовательностью
        - Вызывается demo_fastq_analysis()
        - Ожидается выполнение без исключений
        
    Raises:
        Любое исключение будет проваливать тест
    """
    from fastq_analyzer import demo_fastq_analysis
    
    # Создаем тестовый FASTQ файл
    fastq_content = """@test_read
ATCG
+
IIII
"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
        f.write(fastq_content)
        temp_file = f.name
    
    try:
        # Проверяем, что функция выполняется без ошибок
        demo_fastq_analysis(temp_file)
    finally:
        # Удаляем временный файл
        os.unlink(temp_file)

if __name__ == "__main__":
    # Запуск демонстрационного теста при прямом выполнении файла
    test_demo_fastq_function()
    print("Все тесты FASTQ пройдены!")