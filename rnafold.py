#### RNAhybrid module: compute the features of miR and its target,
#### by Zixing Wang


import subprocess


class RNAFold:
    def __init__(self, _RNAFold):
        self._executable = _RNAFold

    def __call__(self, mrna_seq):
        cmd = "echo %s | %s --noPS" % (mrna_seq, self._executable)
        output = subprocess.check_output(cmd, shell=True)
        return output
