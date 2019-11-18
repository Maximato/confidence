from Bio import SeqIO

estimated_size = 10000


class Sequences:
    @staticmethod
    def extract_from(filename):
        return list(SeqIO.parse(filename, "fasta"))

    @staticmethod
    def group_seqs(lseqs):
        groups = {"cds": []}
        for record in lseqs:
            seq_id = record.id
            group_id = seq_id[0:5]
            length = len(record.seq)
            if length > estimated_size:
                groups["cds"].append(record)
            else:
                if group_id in groups.keys():
                    groups[group_id].append(record)
                else:
                    groups[group_id] = [record]
        return groups

    @staticmethod
    def write_groups(groups):
        for key in groups:
            SeqIO.write(groups[key], "./Groups/" + key + ".fasta", "fasta")


sequences = Sequences()
records = sequences.extract_from("./TBEV_all_nucleotides.fasta")
groups = sequences.group_seqs(records)
sequences.write_groups(groups)
