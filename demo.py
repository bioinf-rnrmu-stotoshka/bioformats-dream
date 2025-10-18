
import os
from vcf_reader import VcfReader


def demo():
    vcf_file = "homo_sapiens-chrMT.vcf"
    
    if not os.path.exists(vcf_file):
        print(f"File not found: {vcf_file}")
        return
    
    try:
        with VcfReader(vcf_file) as reader:
            stats = reader.get_statistics()
            print("File Information:")
            print(f"Variants: {stats.get('total_variants', 0):,}")
            print(f"Samples: {stats.get('samples_count', 0)}")
            print(f"Chromosomes: {reader.get_chromosomes()}")
            print()

            print("First 3 variants:")
            for i, variant in enumerate(reader.read()):
                if i >= 3:
                    break
                print(f"{variant.id}: {variant.ref}→{variant.alt}")
            print()

            region_stats = reader.get_region_stats()
            print("Variants by chromosome:")
            for _, row in region_stats.iterrows():
                print(f"{row['chromosome']}: {row['variant_count']:,}")

            type_stats = reader.get_variant_type_stats()
            print("\nVariant types:")
            for _, row in type_stats.iterrows():
                print(f"{row['variant_type']}: {row['count']:,}")

    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    demo()
