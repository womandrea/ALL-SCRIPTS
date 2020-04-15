import argparse
import Assembler

parser = argparse.ArgumentParser(description="RIGHT NOW: Will take an input directory and convert all the \
fasta files either within the directory or within its subdirectories.\n All outputs will be put in the directories where their reads belong.")

parser.add_argument("-longreads",
                    metavar='input_path', type=str,
                    dest="long", help="Path to directory with all reads.", required=False)
parser.add_argument("-minreads",
                    metavar="READ CUTOFF", type=str,
                    dest="minr", help= "Minimum required read length", default=False, required=True)
parser.add_argument("-shortreads",
                    metavar="illumina_path", type=str,
                    dest="sr", help="Directory containing all short reads in fastq.gz files.", required=Fase)
parser.add_argument("-assembly_type",
                    metavar="[contig, long, short]", required=True,
                    type=str, dest="type", help="Contig corresponds to hybrid assembly, long to long read assembly, and short to read read assembly.")


args = parser.parse_args()
directory = Assembler.AssemblyCall(args)

directory.fasta_files()
directory.assembly()
