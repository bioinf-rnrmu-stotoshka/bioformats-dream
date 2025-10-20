class Record:
    """
    Класс, хранящий биологическую последовательность и информацию
    """

    def __init__(self, id: str):
        self.id = id

    def __repr__(self):
        return f"<Record id={self.id}>"


class SequenceRecord(Record):
    """
    FASTA FASTQ
    """

    def __init__(self, id: str, sequence: str, quality: list[int] | None = None):
        super().__init__(id)
        self.sequence = sequence
        self.quality = quality


class AlignmentRecord(Record):
    """
    SAM
    """

    def __init__(self, id: str, chrom: str, start: int, cigar: str, mapq: int):
        super().__init__(id)
        self.chrom = chrom
        self.start = start
        self.cigar = cigar
        self.mapq = mapq
        self.end: int = start
        self.flag: int = 0

    def __repr__(self):
        return (
            f"<AlignmentRecord id={self.id}, {self.chrom}:{self.start}-{self.end}, "
            f"MAPQ={self.mapq}, FLAG={self.flag}>"
        )

class VariantRecord(Record):
    """
    VCF
    """

    def __init__(self, chrom: str, pos: int, ref: str, alt: str, info: dict):
        super().__init__(f"{chrom}:{pos}")
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.info = info

    def __repr__(self):
        return f"<VariantRecord {self.chrom}:{self.pos} {self.ref}>{self.alt}>"
