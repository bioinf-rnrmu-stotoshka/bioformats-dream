

import pandas as pd
from pathlib import Path
from typing import Iterator, Dict, List, Union
from collections import defaultdict, Counter
import re

from abstract import GenomicDataReader
from record import VariantRecord


class VcfReader(GenomicDataReader):
    
    def __init__(self, filepath: Union[str, Path]):
        super().__init__(filepath)
        self.header_lines = []
        self.metadata = defaultdict(list)
        self.column_names = []
        self.sample_names = []
        self._chromosome_lengths = {}
    
    def _parse_header(self):
        if self._header_parsed:
            return
            
        try:
            if self.file:
                self.file.seek(0)
            
            for line in self.file:
                line = line.strip()
                if line.startswith('##'):
                    self.header_lines.append(line)
                    if '=' in line:
                        key, value = line[2:].split('=', 1)
                        self.metadata[key].append(value)
                        if key == 'contig' and 'length=' in value:
                            chrom_match = re.search(r'ID=([^,>]+)', value)
                            length_match = re.search(r'length=(\d+)', value)
                            if chrom_match and length_match:
                                chrom = chrom_match.group(1)
                                length = int(length_match.group(1))
                                self._chromosome_lengths[chrom] = length
                elif line.startswith('#CHROM'):
                    self.column_names = line[1:].split('\t')
                    if len(self.column_names) > 9:
                        self.sample_names = self.column_names[9:]
                    break
            
            self._header_parsed = True
            self._chromosomes = list(self._chromosome_lengths.keys())
            
        except Exception as e:
            raise RuntimeError(f"Error parsing header: {e}")
    
    def read(self) -> Iterator[VariantRecord]:
        if not self._header_parsed:
            self._parse_header()
        
        if self.file:
            self.file.seek(0)
            for line in self.file:
                if not line.startswith('#'):
                    break
        
        for line in self.file:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 8:
                continue
            
            try:
                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[3]
                alt = parts[4]
                info_str = parts[7] if len(parts) > 7 else ""
                info = self._parse_info_field(info_str)
                
                variant = VariantRecord(chrom=chrom, pos=pos, ref=ref, alt=alt, info=info)
                yield variant
            except Exception:
                continue
    
    @staticmethod
    def _parse_info_field(info_str: str) -> Dict:
        if not info_str or info_str == '.':
            return {}
        
        info_dict = {}
        for item in info_str.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                try:
                    if '.' in value:
                        value = float(value)
                    else:
                        value = int(value)
                except (ValueError, TypeError):
                    pass
                info_dict[key] = value
            else:
                info_dict[item] = True
        
        return info_dict
    
    def get_chromosomes(self) -> List[str]:
        if not self._header_parsed:
            self._parse_header()
        
        if self._chromosomes:
            return self._chromosomes
        
        chromosomes = set()
        current_pos = self.file.tell() if self.file else 0
        
        if self.file:
            self.file.seek(0)
            for line in self.file:
                if not line.startswith('#'):
                    break
            
            for line in self.file:
                if line.startswith('#'):
                    continue
                parts = line.split('\t')
                if parts:
                    chromosomes.add(parts[0])
                if len(chromosomes) > 50:
                    break
        
        if self.file:
            self.file.seek(current_pos)
        
        self._chromosomes = sorted(chromosomes)
        return self._chromosomes
    
    def validate_coordinate(self, chrom: str, pos: int) -> bool:
        if chrom not in self.get_chromosomes():
            return False
        if pos <= 0:
            return False
        if chrom in self._chromosome_lengths and pos > self._chromosome_lengths[chrom]:
            return False
        return True
    
    def get_records_in_region(self, chrom: str, start: int, end: int) -> Iterator[VariantRecord]:
        if not self.validate_coordinate(chrom, start) or not self.validate_coordinate(chrom, end):
            raise ValueError(f"Invalid coordinates: {chrom}:{start}-{end}")
        
        for variant in self.read():
            if variant.chrom == chrom and start <= variant.pos <= end:
                yield variant
    
    def filter_records(self, **filters) -> Iterator[VariantRecord]:
        for variant in self.read():
            if self._variant_matches_filters(variant, filters):
                yield variant
    
    def _variant_matches_filters(self, variant: VariantRecord, filters: Dict) -> bool:
        if 'chrom' in filters and variant.chrom != filters['chrom']:
            return False
        
        if 'min_pos' in filters and variant.pos < filters['min_pos']:
            return False
        if 'max_pos' in filters and variant.pos > filters['max_pos']:
            return False
        
        for filter_key, filter_value in filters.items():
            if filter_key.startswith('info_'):
                info_field = filter_key[5:]
                if info_field not in variant.info:
                    return False
                if variant.info[info_field] != filter_value:
                    return False
        
        return True
    
    def get_statistics(self) -> Dict:
        base_stats = super().get_statistics()
        
        vcf_stats = {
            "file_format": self.metadata.get('fileformat', ['Unknown'])[0] if self.metadata.get('fileformat') else 'Unknown',
            "samples_count": len(self.sample_names),
            "header_lines_count": len(self.header_lines),
        }
        
        variant_count = 0
        chrom_counts = Counter()
        
        current_pos = self.file.tell() if self.file else 0
        if self.file:
            self.file.seek(0)
            for line in self.file:
                if not line.startswith('#'):
                    break
            
            for line in self.file:
                if line.startswith('#'):
                    continue
                parts = line.split('\t')
                if parts:
                    chrom_counts[parts[0]] += 1
                    variant_count += 1
        
        if self.file:
            self.file.seek(current_pos)
        
        vcf_stats.update({
            "total_variants": variant_count,
            "variants_by_chromosome": dict(chrom_counts),
        })
        
        return {**base_stats, **vcf_stats}
    
    def get_region_stats(self) -> pd.DataFrame:
        chrom_counts = Counter()
        
        for variant in self.read():
            chrom_counts[variant.chrom] += 1
        
        data = []
        for chrom, count in chrom_counts.items():
            data.append({
                'chromosome': chrom,
                'variant_count': count,
            })
        
        return pd.DataFrame(data)
    
    def get_variant_type_stats(self) -> pd.DataFrame:
        type_counts = Counter()
        
        for variant in self.read():
            ref_len = len(variant.ref)
            alt_len = len(variant.alt)
            
            if ref_len == 1 and alt_len == 1:
                type_counts['SNP'] += 1
            elif ref_len > alt_len:
                type_counts['Deletion'] += 1
            elif ref_len < alt_len:
                type_counts['Insertion'] += 1
            else:
                type_counts['Complex'] += 1
        
        return pd.DataFrame(list(type_counts.items()), columns=['variant_type', 'count'])


if __name__ == "__main__":
    try:
        with VcfReader("homo_sapiens-chrMT.vcf") as reader:
            variants = list(reader.read())
            print(f"Variants: {len(variants):,}")
            print(f"Chromosomes: {reader.get_chromosomes()}")
    except FileNotFoundError:
        print("File not found")
    except Exception as e:
        print(f"Error: {e}")
