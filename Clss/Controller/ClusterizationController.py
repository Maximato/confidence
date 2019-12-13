from Clss.Model.Clusterization import Clusterization
from Clss.FileSys.Extractor import Extractor
from Clss.FileSys.RecordsWriter import RecordsWriter
from Clss.FileSys.ClusterVisualisation import ClusterVisualisation


class ClusterisationController:

    @staticmethod
    def clusterization(filename, outdir, eps, ms, isdm=None):
        #if isdm:

        clusterisation = Clusterization(Extractor.extract_records(filename))
        clusterisation.clusterize(eps=eps, ms=ms)
        clusters = clusterisation.get_clusters()

        for key in clusters:
            rw = RecordsWriter(clusters[key])
            rw.write_to_dir(f"cluster_{key}.fasta", outdir)

        # visualisation and saving
        cv = ClusterVisualisation(clusterisation.db, clusterisation.dm)
        cv.write_dm(outdir)
        cv.visualize(outdir)
