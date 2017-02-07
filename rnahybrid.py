#### RNAhybrid module: compute the features of miR and its target,
#### by Zixing Wang

import subprocess
import os


class RNAHybrid:
    def __init__(self, _RNAhybrid):
        self._executable = _RNAhybrid

    def __call__(self, mrna_seq, mirna_seq,
                 energy_threshold=None, helix_constraint=None):

        # TODO: Change the dataset type (3utr_human). But do we need to change
        # if we are not interested in the p-value?
        args = "-c -s 3utr_human -m %d" % len(mrna_seq)
        if helix_constraint is not None:
            args += " -f %d,%d" % (helix_constraint[0], helix_constraint[1])
        if energy_threshold is not None:
            args += " -e %0.3f" % energy_threshold


         cmd = "%s %s %s %s" % (self._executable, args, mrna_seq, mirna_seq)
         output = subprocess.check_output(cmd, shell=True)

        return output
