import numpy as np

class Sequence_Node:
    seq = ""
    id = ""

    def __init__(self, seq, id):
        self.seq = seq
        self.id = id
        self.count = np.array(id.rstrip().split(' ')[1].split('-'), dtype=int)
        self.total_count = np.sum(self.count)

    def __hash__(self):
        return hash(self.seq)
