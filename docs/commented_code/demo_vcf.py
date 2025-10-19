import os
from vcf_analyzer import VcfReader


def demo_vcf():
    """
    Демонстрационная функция для работы с VCF файлами.
    
    Основной поток выполнения:
    1. Запрашивает путь к VCF файлу у пользователя
    2. Проверяет существование файла
    3. Читает и анализирует файл с помощью VcfReader
    4. Выводит статистику и примеры вариантов
    
    Исключения:
    - FileNotFoundError: если файл не существует
    - Exception: любые ошибки при чтении/парсинге VCF файла
    """
    # Запрос пути к файлу с консоли
    vcf_file = input("Введите путь к VCF файлу: ")
    
    # Валидация существования файла
    if not os.path.exists(vcf_file):
        print(f"Файл {vcf_file} не найден!")
        print("Убедитесь, что файл находится в той же директории")
        return  # Завершение функции при отсутствии файла
    
    try:
        # Инициализация анализатора VCF
        print(f"Демонстрация VCF ридера")
        print(f"Файл: {vcf_file}")
        print("=" * 60)  # Разделитель для улучшения читаемости вывода
        
        # Использование контекстного менеджера для гарантированного закрытия файла
        with VcfReader(vcf_file) as reader:
            # ========== БЛОК СБОРКИ СТАТИСТИКИ ==========
            print("1. Статистика файла:")
            print("-" * 30)
            
            # Получение различных видов статистики через API VcfReader
            stats = reader.get_statistics()  # Основная статистика
            chromosomes = reader.get_chromosomes()  # Список хромосом
            region_stats = reader.get_region_stats()  # Распределение по регионам
            type_stats = reader.get_variant_type_stats()  # Статистика по типам вариантов
            
            # Форматированный вывод базовой статистики
            print(f"   Всего вариантов: {stats.get('total_variants', 0):,}")  # Разделители тысяч для чисел
            print(f"   Образцы: {stats.get('samples_count', 0)}")  # Количество образцов в файле
            print(f"   Строк заголовка: {stats.get('header_lines', 0)}")  # Мета-информация VCF
            print(f"   Хромосомы: {chromosomes}")  # Список присутствующих хромосом
            print()

            # ========== БЛОК ВЫВОДА ПРИМЕРОВ ВАРИАНТОВ ==========
            print("2. Первые 5 вариантов:")
            print("-" * 30)
            
            variant_count = 0
            # Итерация по вариантам с использованием генератора read()
            for variant in reader.read():
                # Ограничение вывода первыми 5 вариантами
                if variant_count >= 5:
                    break
                
                # Форматирование вывода для читаемости
                info_preview = format_info_preview(variant.info)  # Специальное форматирование INFO
                print(f"   {variant.id}")  # Идентификатор варианта (rs-номер)
                print(f"     REF: {variant.ref} -> ALT: {variant.alt}")  # Референсный и альтернативный аллели
                print(f"     INFO: {info_preview}")  # Дополнительная информация
                print()  # Разделитель между вариантами
                variant_count += 1
            
            # Обработка случая пустого файла
            if variant_count == 0:
                print("   Варианты не найдены")
            print()

            # ========== БЛОК СТАТИСТИКИ ПО ХРОМОСОМАМ ==========
            print("3. Распределение по хромосомам:")
            print("-" * 30)
            
            # Проверка наличия данных перед итерацией
            if not region_stats.empty:  # region_stats - DataFrame pandas
                for _, row in region_stats.iterrows():
                    print(f"   {row['chromosome']}: {row['variant_count']:,} вариантов")
            else:
                print("   Данные по регионам отсутствуют")
            print()

            # ========== БЛОК СТАТИСТИКИ ПО ТИПАМ ВАРИАНТОВ ==========
            print("4. Типы вариантов:")
            print("-" * 30)
            
            if not type_stats.empty:  # type_stats - DataFrame pandas
                total_variants = type_stats['count'].sum()  # Общее количество для расчета процентов
                for _, row in type_stats.iterrows():
                    # Расчет процентного соотношения
                    percentage = (row['count'] / total_variants * 100) if total_variants > 0 else 0
                    print(f"   {row['variant_type']}: {row['count']:,} ({percentage:.1f}%)")
            else:
                print("   Данные по типам вариантов отсутствуют")

    except Exception as e:
        # Обработка ошибок с выводом трассировки стека
        print(f"Ошибка при чтении VCF файла: {e}")
        import traceback
        traceback.print_exc()  # Для диагностики сложных ошибок


def format_info_preview(info_dict):
    """
    Форматирует поля INFO для компактного отображения.
    
    Аргументы:
        info_dict (dict): Словарь с информацией о генетическом варианте
    
    Возвращает:
        str: Отформатированная строка с ключевой информацией
    
    Логика форматирования:
    - Для флаговых значений (True) выводится только ключ
    - Для пар ключ-значение выводится в формате "ключ=значение"
    - Ограничивается вывод первыми 3 полями
    - Добавляется указатель на наличие дополнительных полей
    """
    # Обработка пустого словаря
    if not info_dict:
        return "нет информации"
    
    preview_items = []
    # Ограничение вывода первыми тремя полями для компактности
    for key, value in list(info_dict.items())[:3]:
        if value is True:
            # Случай флаговых параметров (например, "CONFIDENT")
            preview_items.append(key)
        else:
            # Стандартные пары ключ-значение
            preview_items.append(f"{key}={value}")
    
    # Сборка основной строки
    preview = "; ".join(preview_items)
    
    # Добавление индикатора оставшихся полей
    if len(info_dict) > 3:
        preview += f" ... (еще {len(info_dict) - 3} полей)"
    
    return preview


if __name__ == "__main__":
    # Точка входа при прямом запуске скрипта
    demo_vcf()