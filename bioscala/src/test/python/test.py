import numpy as np
namem = ['ZNF615_Imbeault2017_OM_RCADE', 'AIRE_HUMAN.H11MO.0.C']
namem = [x.split('_') for x in namem]

namem = ['_'.join([', '.join(list(x[0].split(', ')))] + x[1:]) for x in namem]
namem = np.array(namem)
# Remove motifs having no TF in current dataset
t1 = [not x.startswith('_') for x in namem]
# dw,dh=[x[:,t1] for x in [dw,dh]]
namem = namem[t1]
if len(namem) != len(set(namem)):
    from collections import Counter
    t1 = [x[0] for x in Counter(namem).items() if x[1] > 1][:3]
    raise ValueError(
        "Found non - unique motif name suffices. Each motif name is recommended to contain a unique suffix. First three non - unique motif names: {}".format(', '.join(t1)))
