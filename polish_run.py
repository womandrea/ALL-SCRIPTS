import subprocess
import os
from psutil import virtual_memory
from math import ceil

# Incorporate "__main__" maybe?


class InputArg:
    def __init__(self, args):
        self.samples = {}
        self.args = args

    def initialize(self):
        """
        :return: Dictionary of samples, including the path to : final_contigs.fasta from Spades/Skesa, the final
        assembly from Shasta and the original reads (long or short)
        """
        short_dir = self.args.shortdir  # Directory with short reads
        longdir = self.args.longdir  # Directory with subdirectories and assemblies

        for sd in os.listdir(longdir):  # will always have a subdir, because of the output of assembly.
            location = os.path.join(longdir, sd)
            f = sd.split("_")[1]
            self.samples[f] = Sample(f)
            self.samples[f].dir = "{}".format(location)

            for files in os.listdir(os.path.join(longdir, sd)):
                if ".fastq" in files:
                    self.samples[f].lr = "{}/{}".format(location, files)

            for r in os.listdir(short_dir):
                if f in r and "R1" in r:
                    self.samples[f].sr1 = "{}{}".format(short_dir, r)
                elif f in r and "R2" in r:
                    self.samples[f].sr2 = "{}{}".format(short_dir, r)

            location = "{}{}".format(longdir, sd)
            contents = os.listdir(location)

            if "spades" in contents:
                max_kmer = max(list(map(lambda x: int(x.split("K")[-1]),
                                        filter(lambda x: "K" in x, os.listdir("{}/spades".format(location))))))

                self.samples[f].short = "{}/spades/K{}/final_contigs.fasta".format(location, max_kmer)
            elif "medaka_consensus" in contents:
                self.samples[f].long = "{}/medaka_consensus/consensus.fasta".format(location)
            elif "ShastaRun" in contents:
                self.samples[f].long = "{}/ShastaRun/Assembly.fasta".format(location)  #
            elif "skesa" in contents:
                self.samples[f].short = "{}/skesa.fasta".format(location)
        return self.samples

    def run_polishing(self):
        samples_dict = self.initialize()

        for f in samples_dict.values():
            f.polishing()


class Sample:
    def __init__(self, name):
        self.short = None  # skesa or Spades
        self.long = None  # Shasta
        self.sr1 = None  # file path to short reads
        self.sr2 = None  # file path to short reads
        self.lr = None  # file path to long reads
        self.dir = None  # The directory where the polished assemblies will go!
        self.name = name
        self.enum = 0

    def polishing(self):
        """
        Calls the appropriate method for filtering.

        :return: None, but calls methods to allow for read polishing
        """
        # if self.short and self.lr:
        # right now only one of them is set
        # use medaka first, and then pilon
        # must switch btw environments
        # is it worth it to do Shasta -> Medaka, Spades -> Pilon

        if "medaka" in self.long or self.short:
            os.mkdir("{}/pilon/".format(self.dir))
            self.iterate_pilon()
            # Use Pilon on skesa/spades assembly

        elif self.long:
            self.medaka()
            # Use Medaka on shasta assembly

    def medaka(self):
        """
        Calls Medaka on Sample class object.

        :return: None
        """
        cmd = ["medaka_consensus", "-i", self.lr, "-d", self.long, "-o", "{}/medaka_consensus".format(self.dir),
               "-m", "r941_min_high_g344"]
        subprocess.run(cmd)

        os.remove("{}/medaka_consensus/calls_to_draft.bam".format(self.dir))
        os.remove("{}/medaka_consensus/calls_to_draft.bam.bai".format(self.dir))
        os.remove("{}/medaka_consensus/consensus_probs.hdf".format(self.dir))

    def iterate_pilon(self):
        """
        Calls pilon x times or until the changes file is empty.
        TO DO: incorporate putting everything in a pilon directory
        :return: None.
        """
        # Be more elegant
        if "medaka" in self.long:
            assembly_file = self.long
        else:
            assembly_file = self.short

        while self.enum < 6:

            if self.enum == 0:
                Sample.pilon(assembly_file, self.dir, self.name, self.sr1, self.sr2, self.enum)
                self.enum += 1
            elif os.stat("{}/pilon/{}_polished_{}.changes".format(self.dir, self.name, (int(self.enum) - 1))).st_size == 0:
                break
            else:
                variable = "{}/pilon/{}_polished_{}.fasta".format(self.dir, self.name, int(self.enum) - 1)
                Sample.pilon(variable, self.dir, self.name, self.sr1, self.sr2, self.enum)
                self.enum += 1

    @staticmethod
    def pilon(assembly, out, name, sr1, sr2, enum):
        """
        Requires that reads are pre-processed before use. Utilizes bwa, samtools, and pilon. Will produce the
        required bam file and index the assemblies. Requires being run in a virtual environment with bwa, samtools,
        and pilon installed.

        :return: None.
        """
        limit = ceil(virtual_memory().total / (10 ** 9) * .65)

        if not os.path.isdir(out):
            os.mkdir(out)
        # add pilon directory!
        rmdup_name = "{}/{}_nd_{}.bam".format(out, name, enum)
        final_bam = "{}/{}_final_{}.bam".format(out, name, enum)
        polished_file = "{}/pilon/{}_polished_{}".format(out, name, enum)
        print(assembly)
        # index assembly
        bw_index_cmd = ["bwa", "index", "{}".format(assembly)]
        # align short reads to indexed assembly
        bwa_mem_cmd = ["bwa", "mem", "-M", "-t", "{}".format(len(os.sched_getaffinity(0))), assembly, sr1, sr2]
        # remove unmapped reads and return a bam file with the appropriate header
        st_view_cmd = ["samtools", "view", "-h", "-b", "-F", "4"]
        # sort .bam file based on position
        st_sort_cmd = ["samtools", "sort", "-o", rmdup_name]
        # remove duplicates from bam file
        rm_dup_cmd = ["samtools", "rmdup", rmdup_name, final_bam]
        # index final bam file
        st_index = ["samtools", "index", final_bam]
        # command for Pilon
        pilon_cmd = ["java", "-Xmx{}G".format(limit),
                     "-jar", "pilon-1.23.jar",
                     "--genome", assembly,
                     "--frags", final_bam,
                     "--vcf",
                     "--output", polished_file,
                     "--threads", "{}".format(len(os.sched_getaffinity(0))),
                     "--changes"]

        # runs listed commands
        subprocess.run(bw_index_cmd)
        bwa = subprocess.Popen(bwa_mem_cmd,
                               stdout=subprocess.PIPE)  # pipe out bwa alignment with short reads to assembly
        st_view = subprocess.Popen(st_view_cmd, stdin=bwa.stdout,
                                   stdout=subprocess.PIPE)  # pipe out unmapped reads removed bam file
        bwa.stdout.close()
        st_sort = subprocess.Popen(st_sort_cmd, stdin=st_view.stdout, stdout=subprocess.PIPE)
        st_view.stdout.close()
        st_sort.communicate()[0]
        subprocess.run(rm_dup_cmd)
        subprocess.run(st_index)
        subprocess.run(pilon_cmd)

        # remove temporary files produced
        os.remove(final_bam)
        os.remove(rmdup_name)
        os.remove("{}.amb".format(assembly))
        os.remove("{}.ann".format(assembly))
        os.remove("{}.pac".format(assembly))
        os.remove("{}.bwt".format(assembly))
        os.remove("{}.sa".format(assembly))
        os.remove("{}.bai".format(final_bam))
