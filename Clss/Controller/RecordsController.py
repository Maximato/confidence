from Clss.Extractor.Extractor import Extractor
from Clss.FileSys.RecordsWriter import RecordsWriter
from Clss.Model.Records import Records


class RecordsController:

    @staticmethod
    def grouping(filename, size_of_genome, outdir):
        records = Records(Extractor.recs_extractor(filename))

        groups = records.group(size_of_genome)
        for key in groups:
            rw = RecordsWriter(groups[key])
            rw.write_to_dir(key + ".fasta", outdir)

    @staticmethod
    def filtrating(filename, out_file, organism, minsize, maxsize):
        records = Records(Extractor.recs_extractor(filename))
        fseqs = records.filtr_organism_by_size(organism, minsize, maxsize)
        RecordsWriter(fseqs).write_to(out_file)
