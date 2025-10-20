import os
from vcf_reader import VcfReader


def demo_vcf():
    """
    Демонстрация работы с VCF ридером на файле homo_sapiens-chrMT.vcf
    """
    vcf_file = input("Введите путь к VCF файлу: ")
    
    if not os.path.exists(vcf_file):
        print(f"Файл {vcf_file} не найден!")
        print("Убедитесь, что файл находится в той же директории")
        return
    
    try:
        print(f"Демонстрация VCF ридера")
        print(f"Файл: {vcf_file}")
        print("=" * 60)
        
        with VcfReader(vcf_file) as reader:
            # Получаем всю статистику
            print("1. Статистика файла:")
            print("-" * 30)
            
            stats = reader.get_statistics()
            chromosomes = reader.get_chromosomes()
            region_stats = reader.get_region_stats()
            type_stats = reader.get_variant_type_stats()
            
            print(f"   Всего вариантов: {stats.get('total_variants', 0):,}")
            print(f"   Образцы: {stats.get('samples_count', 0)}")
            print(f"   Строк заголовка: {stats.get('header_lines', 0)}")
            print(f"   Хромосомы: {chromosomes}")
            print()

            # Показываем первые несколько вариантов
            print("2. Первые 5 вариантов:")
            print("-" * 30)
            
            variant_count = 0
            for variant in reader.read():
                if variant_count >= 5:
                    break
                
                # Форматируем информацию для красивого вывода
                info_preview = format_info_preview(variant.info)
                print(f"   {variant.id}")
                print(f"     REF: {variant.ref} -> ALT: {variant.alt}")
                print(f"     INFO: {info_preview}")
                print()
                variant_count += 1
            
            if variant_count == 0:
                print("   Варианты не найдены")
            print()

            # Статистика по регионам
            print("3. Распределение по хромосомам:")
            print("-" * 30)
            
            if not region_stats.empty:
                for _, row in region_stats.iterrows():
                    print(f"   {row['chromosome']}: {row['variant_count']:,} вариантов")
            else:
                print("   Данные по регионам отсутствуют")
            print()

            # Статистика по типам вариантов
            print("4. Типы вариантов:")
            print("-" * 30)
            
            if not type_stats.empty:
                total_variants = type_stats['count'].sum()
                for _, row in type_stats.iterrows():
                    percentage = (row['count'] / total_variants * 100) if total_variants > 0 else 0
                    print(f"   {row['variant_type']}: {row['count']:,} ({percentage:.1f}%)")
            else:
                print("   Данные по типам вариантов отсутствуют")

    except Exception as e:
        print(f"Ошибка при чтении VCF файла: {e}")
        import traceback
        traceback.print_exc()


def format_info_preview(info_dict):
    """
    Форматирует INFO поле для красивого вывода
    """
    if not info_dict:
        return "нет информации"
    
    preview_items = []
    for key, value in list(info_dict.items())[:3]:  # Показываем первые 3 поля
        if value is True:
            preview_items.append(key)
        else:
            preview_items.append(f"{key}={value}")
    
    preview = "; ".join(preview_items)
    if len(info_dict) > 3:
        preview += f" ... (еще {len(info_dict) - 3} полей)"
    
    return preview


if __name__ == "__main__":
    demo_vcf()
