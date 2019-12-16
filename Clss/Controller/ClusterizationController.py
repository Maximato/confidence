from Clss.Model.Clusterization import Clusterization
from Clss.FileSys.Extractor import Extractor
from Clss.FileSys.RecordsWriter import RecordsWriter
from Clss.FileSys.ClusterVisualisation import ClusterVisualisation
from Clss.Model.Records import Records


class ClusterisationController:

    @staticmethod
    def clusterization(filename, outdir, eps, ms, isdm=None):

        if isdm:
            dm = Extractor.extract_dist_matrix(filename)
        else:
            records = Records(Extractor.extract_records(filename))
            dm = records.create_dist_matrix()

        c = Clusterization(dm)
        ci = c.get_clusters_indexes(eps=eps, ms=ms)
        coordinates = c.get_coordinates()

        if not isdm:
            clusters = records.get_clusters(ci)
            for key in clusters:
                rw = RecordsWriter(clusters[key])
                rw.write_to_dir(f"cluster_{key}.fasta", outdir)

        # visualisation and saving
        cv = ClusterVisualisation(coordinates, dm, ci)
        cv.write_dm(outdir)
        cv.visualize(outdir)
