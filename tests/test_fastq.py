import pytest
import tempfile
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from bio_analyzer.fastq_analyzer import FastqAnalyzer


class TestFastqAnalyzer:
    """Тесты для FASTQ анализатора."""
    
    def create_test_fastq(self, content: str) -> str:
        """Создает временный FASTQ файл."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fastq') as f:
            f.write(content)
            return f.name
    
    def test_fastq_count_sequences(self):
        """Тест подсчета последовательностей в FASTQ."""
        fastq_content = """@read1
ATCG
+
IIII
@read2
GGCC
+
IIII
"""
        test_file = self.create_test_fastq(fastq_content)
        
        try:
            with FastqAnalyzer(test_file) as analyzer:
                sequences = list(analyzer.read())
                assert analyzer.get_seq_count() == 2
        finally:
            os.unlink(test_file)
    
    def test_fastq_mean_length(self):
        """Тест расчета средней длины в FASTQ."""
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
                assert analyzer.get_mean_seq_length() == (4 + 6) / 2
        finally:
            os.unlink(test_file)
    
    def test_phred_quality_conversion(self):
        """Тест конвертации Phred качества."""
        analyzer = FastqAnalyzer("dummy.fastq")
        # Символ 'I' имеет ASCII код 73, при offset 33: 73-33=40
        assert analyzer.phred_to_quality('I') == 40


def test_demo_fastq_function():
    """Тест демонстрационной функции."""
    from bio_analyzer.fastq_analyzer import demo_fastq_analysis
    
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
        os.unlink(temp_file)

if __name__ == "__main__":
    test_demo_fastq_function()
    print("Все тесты FASTQ пройдены!")
