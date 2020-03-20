import argparse
import polish_run

parser = argparse.ArgumentParser(description="To be used for guppybase calling")

parser.add_argument("-illumina ", metavar='PATH',
                    type=str, dest="shortdir", help='Directory with all short reads, if available')
parser.add_argument("-minion", metavar="PATH",
                    type=str, dest="longdir", help="Directory with output(s) from assembly.py.")
args = parser.parse_args()

assmb = polish_run.InputArg(args)
assmb.run_polishing()
