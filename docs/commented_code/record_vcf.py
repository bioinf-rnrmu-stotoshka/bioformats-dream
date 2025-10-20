class RecordVcf:
    """
    Базовый класс для хранения биологических данных с идентификатором.
    
    Используется как родительский класс для различных типов биологических записей.
    Содержит минимальную общую информацию - идентификатор записи.
    
    Attributes
    ----------
    id : str
        Уникальный идентификатор записи в исходном файле данных
    
    Examples
    --------
    >>> record = RecordVcf("seq1")
    >>> print(record)
    <RecordVcf id=seq1>
    """

    def __init__(self, id: str):
        """
        Инициализирует запись с идентификатором.
        
        Parameters
        ----------
        id : str
            Идентификатор записи (например, название последовательности 
            или координаты варианта)
        """
        self.id = id

    def __repr__(self):
        """
        Возвращает строковое представление объекта для отладки.
        
        Returns
        -------
        str
            Строка в формате: <RecordVcf id={идентификатор}>
        """
        return f"<RecordVcf id={self.id}>"


class SequenceRecordVcf(RecordVcf):
    """
    Класс для хранения последовательностей в форматах FASTA/FASTQ.
    
    Наследует от RecordVcf и добавляет информацию о биологической 
    последовательности и качестве оснований.
    
    Attributes
    ----------
    sequence : str
        Биологическая последовательность (нуклеотиды или аминокислоты)
    quality : list[int] | None
        Список значений качества для каждого символа последовательности
        в формате Phred (только для FASTQ). None для FASTA записей.
    
    Examples
    --------
    FASTA запись:
    >>> seq_rec = SequenceRecordVcf("seq1", "ATCG")
    >>> print(seq_rec.sequence)
    ATCG
    
    FASTQ запись:
    >>> qual = [40, 40, 39, 38]
    >>> fastq_rec = SequenceRecordVcf("read1", "ATCG", qual)
    >>> print(fastq_rec.quality)
    [40, 40, 39, 38]
    """

    def __init__(self, id: str, sequence: str, quality: list[int] | None = None):
        """
        Инициализирует запись последовательности.
        
        Parameters
        ----------
        id : str
            Идентификатор последовательности
        sequence : str
            Строка с биологической последовательностью
        quality : list[int] | None, optional
            Список значений качества, по умолчанию None
        """
        super().__init__(id)
        self.sequence = sequence
        self.quality = quality


class AlignmentRecordVcf(RecordVcf):
    """
    Класс для хранения записей выравнивания в формате SAM/BAM.
    
    Содержит информацию о позиции выравнивания последовательности
    на референсный геном и параметрах качества выравнивания.
    
    Attributes
    ----------
    chrom : str
        Название хромосомы/контига референсного генома
    start : int
        Позиция начала выравнивания (0-based индексация)
    cigar : str
        CIGAR строка, описывающая соответствие между последовательностью
        и референсом (M - match, I - insertion, D - deletion и т.д.)
    mapq : int
        Качество выравнивания (Mapping Quality)
    end : int
        Позиция конца выравнивания (вычисляется автоматически)
    flag : int
        SAM флаг выравнивания (по умолчанию 0)
    
    Examples
    --------
    >>> aln = AlignmentRecordVcf("read1", "chr1", 1000, "50M", 60)
    >>> print(aln)
    <AlignmentRecordVcf id=read1, chr1:1000-1050, MAPQ=60, FLAG=0>
    """

    def __init__(self, id: str, chrom: str, start: int, cigar: str, mapq: int):
        """
        Инициализирует запись выравнивания.
        
        Parameters
        ----------
        id : str
            ID прочитанной последовательности (QNAME в SAM)
        chrom : str
            Референсная хромосома (RNAME в SAM)
        start : int
            Левая позиция выравнивания (POS в SAM)
        cigar : str
            CIGAR строка выравнивания
        mapq : int
            Качество отображения (MAPQ в SAM)
        """
        super().__init__(id)
        self.chrom = chrom
        self.start = start
        self.cigar = cigar
        self.mapq = mapq
        self.end: int = start  # Должно вычисляться на основе CIGAR
        self.flag: int = 0

    def __repr__(self):
        """
        Возвращает детальное строковое представление выравнивания.
        
        Returns
        -------
        str
            Строка в формате: <AlignmentRecordVcf id={id}, {chrom}:{start}-{end}, 
            MAPQ={mapq}, FLAG={flag}>
        """
        return (
            f"<AlignmentRecordVcf id={self.id}, {self.chrom}:{self.start}-{self.end}, "
            f"MAPQ={self.mapq}, FLAG={self.flag}>"
        )


class VariantRecordVcf(RecordVcf):
    """
    Класс для хранения генетических вариантов в формате VCF.
    
    Содержит информацию о позиции варианта, референсном и альтернативном 
    аллелях, а также дополнительную информацию в формате словаря.
    
    Attributes
    ----------
    chrom : str
        Название хромосомы
    pos : int
        Позиция варианта (1-based)
    ref : str
        Референсный аллель
    alt : str
        Альтернативный аллель
    info : dict
        Словарь с дополнительной информацией о варианте
        (соответствует полю INFO в VCF)
    
    Examples
    --------
    >>> info = {"DP": 50, "AF": 0.25}
    >>> var = VariantRecordVcf("chr1", 123456, "A", "T", info)
    >>> print(var)
    <VariantRecordVcf chr1:123456 A>T>
    """

    def __init__(self, chrom: str, pos: int, ref: str, alt: str, info: dict):
        """
        Инициализирует запись варианта.
        
        Parameters
        ----------
        chrom : str
            Название хромосомы
        pos : int
            Позиция на хромосоме
        ref : str
            Референсная последовательность в данной позиции
        alt : str
            Альтернативная последовательность
        info : dict
            Дополнительная информация о варианте
        """
        super().__init__(f"{chrom}:{pos}")
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.info = info

    def __repr__(self):
        """
        Возвращает компактное представление варианта.
        
        Returns
        -------
        str
            Строка в формате: <VariantRecordVcf {chrom}:{pos} {ref}>{alt}>
        """
        return f"<VariantRecordVcf {self.chrom}:{self.pos} {self.ref}>{self.alt}>"