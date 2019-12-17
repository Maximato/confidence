from Clss.FileSys.Extractor import Extractor
from Clss.FileSys.RecordsWriter import RecordsWriter
from Clss.Model.Records import Records


class RecordsController:

    @staticmethod
    def grouping(filename, outdir, minsog, maxsog):
        """
        Grouping sequences by names. All big sequences (genomes) with minsog <= size <= maxsog form separate group 'cds'

        :param filename: filename with sequences in fasta
        :param outdir: output directory for saving groups
        :param minsog: int, min size of genome
        :param maxsog: int, max size of genome
        """
        records = Records(Extractor.extract_records(filename))

        groups = records.group(minsog, maxsog)
        for key in groups:
            rw = RecordsWriter(groups[key])
            rw.write_to_dir(key + ".fasta", outdir)

    @staticmethod
    def filtrating(filename, out_file, organism, minsize, maxsize):
        """
        Filtrating sequences by parameters: organizm, minsize and maxsize

        :param filename: filename with sequences in fasta
        :param out_file: out filename
        :param organism: full name of organism
        :param minsize: minimal size of sequences
        :param maxsize: maximal size of sequences
        """
        records = Records(Extractor.extract_records(filename))
        fseqs = records.filtr_organism_by_size(organism, minsize, maxsize)
        RecordsWriter(fseqs).write_to(out_file)
