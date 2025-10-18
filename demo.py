# test_run.py
from sam_reader import SamReader

def test_sam(file_path):
    """Быстрая проверка SAM ридера"""
    print("Testing SAM reader...")
    try:
        with SamReader(file_path) as reader:
            # Базовые методы
            header = reader.get_header()
            count = reader.count_alignments()
            chroms = reader.get_chromosomes()
            
            print(f"Header: {len(header)} lines")
            print(f"Alignments: {count}")
            print(f"Chromosomes: {chroms[:3]}")  # первые 3
            
            # Первая запись
            for rec in reader.read():
                print(f"First record: {rec.chrom}:{rec.start}")
                break
                
    except Exception as e:
        print(f"SAM Error: {e}")

if __name__ == "__main__":
    # Укажи путь к своему SAM файлу
    sam_file = "Col0_C1.100k.sam"  # поменяй на свой SAM
    
    test_sam(sam_file)
