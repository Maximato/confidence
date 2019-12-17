from Clss.Model.Clusterization import Clusterization
from Clss.FileSys.Extractor import Extractor
from Clss.FileSys.RecordsWriter import RecordsWriter
from Clss.FileSys.ClusterVisualisation import ClusterVisualisation
from Clss.Model.Records import Records


class ClusterisationController:

    @staticmethod
    def clusterization(filename, outdir, eps, ms, dm=None):
        """
        Clusterization of sequences. Visualisation and saving of clusters into files.

        :param filename: filename of sequences in fasta format to for clusterisation
        :param outdir: out directory name for out files
        :param eps: The maximum distance between two samples
        :param ms: The number of samples (or total weight) in a neighborhood for a point to be considered
        as a core point. This includes the point itself.
        :param dm: path to distance matrix (optional). If is not None dm will be used as distance matrix.
        """
        records = Records(Extractor.extract_records(filename))
        if dm is None:
            dm = records.create_dist_matrix()
        else:
            dm = Extractor.extract_dist_matrix(dm)

        c = Clusterization(dm)
        ci = c.get_clusters_indexes(eps=eps, ms=ms)
        coordinates = c.get_coordinates()

        clusters = records.get_clusters(ci)
        for key in clusters:
            rw = RecordsWriter(clusters[key])
            rw.write_to_dir(f"cluster_{key}.fasta", outdir)

        # visualisation and saving
        cv = ClusterVisualisation(coordinates, dm, ci)
        cv.write_dm(outdir)
        cv.visualize(outdir)
