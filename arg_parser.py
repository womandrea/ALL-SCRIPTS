import argparse
import Assembler

parser = argparse.ArgumentParser(description="RIGHT NOW: Will take an input directory and convert all the \
fasta files either within the directory or within its subdirectories.")

parser.add_argument("-longreads",
                    metavar='input path', type=str,
                    dest="long", help="Path to directory.")
parser.add_argument("-shasta_out",
                    metavar="output directory", type=str,
                    dest="output", help="Output directory for the fasta files.")  # DO NOT USE THIS PARAM, remove
parser.add_argument("-minreads",
                    metavar="READ CUTOFF", type=str,
                    dest="minr", help= "Minimum required read length")
parser.add_argument("-shortreads",
                    metavar="FORWARD", type=str,
                    dest="sr", help="Directory containing all short reads in fastq.gz files.")
parser.add_argument("-type",
                    metavar="[contig, nanpore, short]",
                    type=str, dest="type", help="If we used SPAdes, what do we do with the long reads?")


args = parser.parse_args()
directory = Assembler.AssemblyCall(args)

directory.fasta_files()
directory.assembly()
