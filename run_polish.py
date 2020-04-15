import argparse
import polisher

parser = argparse.ArgumentParser(description="Assumes that corresponding assembly script was used for assembly of the genome. Polishing done with either short or long reads.")

parser.add_argument("-illumina ", metavar='PATH',
                    type=str, dest="shortdir", help='Directory with all short reads.')
parser.add_argument("-minion", metavar="PATH",
                    type=str, dest="longdir", help="Directory with output(s) from Assembly.py")
args = parser.parse_args()

assmb = polish_run.InputArg(args)
assmb.run_polishing()
