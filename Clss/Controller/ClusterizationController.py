from Clss.Model.Clusterization import Clusterization
from Clss.FileSys.Extractor import Extractor
from Clss.FileSys.RecordsWriter import RecordsWriter


class ClusterisationController:

    @staticmethod
    def clusterization(filename, outdir, eps=1, ms=2):
        clusterisation = Clusterization(Extractor.extract_records(filename))
        clusterisation.clusterize(eps=eps, ms=ms)
        clusters = clusterisation.get_clusters()

        for key in clusters:
            rw = RecordsWriter(clusters[key])
            rw.write_to_dir(f"cluster_{key}.fasta", outdir)


ClusterisationController.clusterization("../../Data/Groups/TBEV_sequences/HE.fasta", "df")
