import sys
import os
import tempfile

"""
МОДУЛЬ ТЕСТИРОВАНИЯ FASTAAnalyzer

Этот модуль содержит тесты для проверки функциональности класса FastaAnalyzer.
Использует временные файлы для изолированного тестирования парсера FASTA.
"""

# Добавляем путь к bio_analyzer в Python path
# Это необходимо для импорта модулей из родительского каталога
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Импортируем тестируемый класс после корректировки sys.path
from fasta_analyzer import FastaAnalyzer


def test_basic_functionality():
    """
    Тест базовой функциональности FastaAnalyzer.
    
    Проверяет:
    1. Корректность парсинга последовательностей
    2. Подсчет количества последовательностей
    3. Расчет средней длины последовательностей
    
    Процесс тестирования:
    1. Создает временный FASTA-файл с тестовыми данными
    2. Использует FastaAnalyzer для обработки файла
    3. Проверяет соответствие ожидаемым результатам
    4. Автоматически очищает временные ресурсы
    """
    # Тестовые данные в формате FASTA
    fasta_content = """>seq1
ATCG
>seq2
GGCC
"""
    # Создание временного файла
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as f:
        f.write(fasta_content)
        temp_file = f.name  # Сохраняем путь к временному файлу
    
    try:
        # Инициализация анализатора с использованием менеджера контекста
        with FastaAnalyzer(temp_file) as analyzer:
            # Чтение и материализация всех последовательностей
            sequences = list(analyzer.read())
            
            # Проверка количества последовательностей
            assert analyzer.get_seq_count() == 2, "Ожидается 2 последовательности"
            
            # Проверка средней длины (4+4)/2=4.0
            assert analyzer.get_mean_seq_length() == 4.0, "Ожидается средняя длина 4.0"
            
            print("Тест пройден! Базовая функциональность работает корректно.")
            
    finally:
        # Обязательная очистка: удаление временного файла
        os.unlink(temp_file)


if __name__ == "__main__":
    """
    Основная точка входа при прямом запуске скрипта.
    Запускает тест и автоматически сообщает о результате.
    """
    test_basic_functionality()