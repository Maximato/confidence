class Records:
    def __init__(self, records):
        self.records = records

    def group(self, estimated_size=10000):
        groups = {"cds": []}
        for record in self.records:
            rec_id = record.id
            group_id = rec_id[0:2]
            length = len(record.seq)
            if length > estimated_size:
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
