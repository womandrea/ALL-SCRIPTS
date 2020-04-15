import filtering
import argparse

parser = argparse.ArgumentParser(description="Providing long reads or short reads, filters reads for quality.")
parser.add_argument("-longreads", metavar="DIR PATH", dest="LR", type=str,help="Full file path to directory with long reads.", default=False)
parser.add_argument("-shortreads", metavar="DIR PATH", dest="SR", type=str, help="ull file path to directory with Illumina reads", default=False)
parser.add_argument("-genomesize", metavar="NUM", dest="genomesize", help="Use integer values for size of genome.", default=False)
parser.add_argument("-filter", action="store_true", default=False, dest="shortreads", help="Specify if short reads should be filtered with bbduk before acting as a reference for Filtlong.")
parser.add_argument("-minlen", metavar="NUM", dest="minlen", help="For Filtlong only: minimum length of reads to filter.")

args = parser.parse_args()

process = filtering.FilteringInput(args)
process.filtering()
