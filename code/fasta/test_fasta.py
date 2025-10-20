import sys
import os
import tempfile

# Добавляем путь к bio_analyzer в Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fasta.fasta_analyzer import FastaAnalyzer

def test_basic_functionality():
    """Простой тест базовой функциональности."""
    fasta_content = """>seq1
ATCG
>seq2
GGCC
"""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.fasta') as f:
        f.write(fasta_content)
        temp_file = f.name
    
    try:
        with FastaAnalyzer(temp_file) as analyzer:
            sequences = list(analyzer.read())
            assert analyzer.get_seq_count() == 2
            assert analyzer.get_mean_seq_length() == 4.0
            print("Тест пройден!")
    finally:
        os.unlink(temp_file)

if __name__ == "__main__":
    test_basic_functionality()
