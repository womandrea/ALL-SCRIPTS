import filtering
import argparse

parser = argparse.ArgumentParser(description="Providing long reads or short reads, filters reads for quality")
parser.add_argument("-longreads", metavar="DIR PATH", dest="LR", type=str,help="Default None")
parser.add_argument("-shortreads", metavar="DIR PATH", dest="SR", type=str, help="Default None")
parser.add_argument("-genomesize", metavar="NUM", dest="genomesize", help="Use integer values for size of genome. Default None.")
parser.add_argument("-filter", metavar="Y/N", dest="shortreads", help="Process short reads with bbduk for filtlong?")
parser.add_argument("-minlen", metavar="NUM", dest="minlen", help="filtlong only: what is the min length ")

args = parser.parse_args()

process = filtering.FilteringInput(args)
process.filtering()
