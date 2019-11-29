class Records:
    def __init__(self, records):
        self.records = records

    def get_seqs(self):
        seqs = []
        for record in self.records:
            seqs.append(record.seq)
        return seqs

    def group(self, genome_min=9000, genome_max=12000):
        groups = {"cds": []}
        for record in self.records:
            rec_id = record.id
            group_id = rec_id[0:2]
            length = len(record.seq)
            if genome_min < length < genome_max:
                groups["cds"].append(record)
            else:
                if group_id in groups.keys():
                    groups[group_id].append(record)
                else:
                    groups[group_id] = [record]
        print("Number of groups: ", len(groups))
        return groups

    def filtr_organism_by_size(self, organizm, minsize, maxsize):
        print("Initially number of sequences: ", len(self.records))
        fseqs = []
        for record in self.records:
            description = record.description.lower()
            if organizm in description and (minsize <= len(record.seq) <= maxsize):
                if ("chimeric" not in description) and ("chimera" not in description):
                    fseqs.append(record)
        print("Number of sequences after filtrating: ", len(fseqs))
        return fseqs
