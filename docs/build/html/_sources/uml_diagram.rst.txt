UML Диаграмма отношения классов в программе
===========================================

.. mermaid::

   classDiagram
    class Reader {
        <<abstract>>
        #filename: string
        +__init__(filename: string)
        +read() Iterator[Record]
        +close()
        #parse_line(line: str) Record
    }

    class SequenceReader {
        <<abstract>>
        +get_sequence(seq_id: str) Sequence
        +validate_sequence(sequence: str) bool
    }

    class GenomicDataReader {
        <<abstract>>
        +get_chromosomes() List[str]
        +get_reference_genome() string
        +validate_coordinate(chrom: str, pos: int) bool
    }

    class FastaReader {
        +read_sequences() Dict[str, Sequence]
        +get_sequence_length(seq_id: str) int
    }

    class FastqReader {
        +get_quality_scores(seq_id: str) List[int]
        +get_average_quality(seq_id: str) float
    }

    class SamReader {
        +read_alignments() List[Alignment]
        +get_header() SamHeader
        +filter_alignments(flag: int) List[Alignment]
        +calculate_coverage(chrom: str) Dict[int, int]

    }

    class VcfReader {
        +read_variants() List[Variant]
        +get_header() VcfHeader
        +filter_by_quality(min_qual: float) List[Variant]
        +get_genotype(sample: str, variant: Variant) Genotype
    }

    Reader <|-- SequenceReader
    Reader <|-- GenomicDataReader
    SequenceReader <|-- FastaReader
    SequenceReader <|-- FastqReader
    GenomicDataReader <|-- SamReader
    GenomicDataReader <|-- VcfReader
