from Clss.FileSys.Extractor import Extractor
from Clss.FileSys.RecordsWriter import RecordsWriter
from Clss.Model.Records import Records


class RecordsController:

    @staticmethod
    def grouping(filename, minsog, maxsog, outdir):
        # minsog - min size of genome
        # maxsog - max size of genome
        records = Records(Extractor.extract_records(filename))

        groups = records.group(minsog, maxsog)
        for key in groups:
            rw = RecordsWriter(groups[key])
            rw.write_to_dir(key + ".fasta", outdir)

    @staticmethod
    def filtrating(filename, out_file, organism, minsize, maxsize):
        records = Records(Extractor.extract_records(filename))
        fseqs = records.filtr_organism_by_size(organism, minsize, maxsize)
        RecordsWriter(fseqs).write_to(out_file)
