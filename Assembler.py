import os
import subprocess
import gzip
import shutil


class AssemblySample:
    """
    Sample class possessing features indicating relevant information for skesa, shasta, or spades assembly.
    """
    def __init__(self, args):
        self.longread = False  # path file to the fastq.gz long read file
        self.fasta = False  # updated only for shasta, the path file for the fasta version of the longread samples
        self.r1 = False # path file to the for the forward illumina reads
        self.r2 = False  # path file to for the reverse illumina reads
        self.args = args  # arguments from arg_parse
        self.type = args.type  # type of assembly [long, contig, short]
        if args.minr:
            self.minread = args.minr  # Sets only if set by the user
        else:
            self.minread = 0
        self.contigs = False  # updated for spades assembly, the output of the shasta long read
        self.name = None  # based on the sample name, used for file naming

    @staticmethod
    def shasta_run(reads, args):
        """
        Will take a file and call shasta on it.

        :param reads: type Str: File path for the fasta used in assembly.
        :param args: type args: From argparser
        :returns: None. Output for all shasta runs. Shasta will be run in the output directory. Requires that shasta is downloaded and available.
        """
        sample_path = "{}/{}".format("/".join(reads.split("/")[:-1]), "ShastaRun")
        # IMPROVEMENT: Identify what shasta version is used.
        cmd = ["shasta-Linux-0.1.0",
               "--input", reads,
               "--memoryMode", "filesystem",
               "--memoryBacking", "2M",
               "--Reads.minReadLength", args.minr,
               "--output", sample_path
               ]

        subprocess.run(cmd)
        os.remove(reads)

    @staticmethod
    def spades_run(r1, r2, long, assmb):
        """
        :param r1: type Str: File path to the forward reads
        :param r2: type Str: File path to the reverse reads
        :param long: type Str: the output directory path
        :param assmb: type Str/Bool: the path to the scaffold assembly, if produced by Shasta. If not, a boolean (False). 
        :return:
        """

        output_dir = "{}/spades/".format("/".join(long.split("/")[:-1]))
        if not os.path.isdir(output_dir):   
            os.mkdir(output_dir)
        cmd = ["python3", "/home/bioinfo/SPAdes-3.12.0-Linux/bin/spades.py",
               "-1", r1,
               "-2", r2,
               "-o", output_dir,
               "--only-assembler"]
        if assmb:
            cmd.extend(["--trusted-contigs", assmb])

        subprocess.run(cmd)

    @staticmethod
    def skesa(r1, r2, output):
        """
        :param r1: type Str: File path to the forward reads.
        :param r2: type Str: File path to the reverse reads.
        :param output: type Str: Output file path
        :return: None. 
        """
        cmd = ["skesa", "--fastq", "{},{}".format(r1, r2), "--contigs_out", output]
        subprocess.run(cmd)

        sample_name = output.split("/")[-1].split("_")[0]
        path = '/'.join(output.split("/")[:-1])
        for sd in [f for f in os.listdir(path) if os.path.isdir("{}/{}".format(path, f))]:
            if sample_name in sd:
                os.rename(output, "{}/{}/{}".format(path, sd, output.split("/")[-1]))

    def assembly_run(self):
        """
        Decides which series of assembly stuff to run based on args, using long, r1, r2 or assembly as an input
        :return: None.
        """
        # RUN THIS TO MAKE SURE
        if self.type == "long":
            AssemblySample.shasta_run(self.fasta, self.args)
        elif self.type.lower() == "contig":
            AssemblySample.shasta_run(self.fasta, self.args)
            self.contigs = "{}filtered_{}/ShastaRun/Assembly.fasta".format(self.longread, self.name)  # might not be the right file path
            AssemblySample.spades_run(self.r1, self.r2, self.fasta, self.contigs)
        elif self.type.lower() == "short":
            file_name = "{}/{}_skesa_assembly".format(os.getcwd(),self.name)
            AssemblySample.skesa(self.r1, self.r2, file_name)
        else:
            print("Please input long for shasta assembly, contig for hybrid, and short for skesa assembly")


class AssemblyCall:
    """
    The class with all the information for the assembly call. Does not distinguish between the samples in the input
    directories.
    """

    def __init__(self, args):
        # self.assembly = False
        if args.long:
            self.dirlong = args.long  # directory with all the long read data
        else: 
            self.dirlong = False
        if args.sr:
            self.dirshort = args.sr  # default None, directory with all the short reads (forward and reverse)
        else:
            self.dirshort = False
        self.args = args

    @staticmethod
    def samples_to_run(directory, sample_dir):

        """
        Sorts reads in the directories with short reads and long reads,
        using the sample name and features to create a dictionary of Samples
        :param directory: directory with short reads
        :param sample_dir: dictionary to be added to.
        :return: Dictionary
        """
        illum_samples = set(map(lambda x: x.split('_')[0], os.listdir(directory)))
        for s in illum_samples:
            if s in sample_dir:
                for r in os.listdir(directory):
                    if s in r and "R1" in r:
                        sample_dir[s].r1 = "{}{}".format(directory, r)
                    elif s in r and "R2" in r:
                        sample_dir[s].r2 = "{}{}".format(directory, r)
        return sample_dir

    def fasta_finder(self):

        """
        Go through subdirectories with copies of both fasta and fastq files. Will add fasta pathfile
        to the sample class for shasta input.

        :return: Dict

        FIX: Relies on short reads being provided. The issue is the initialization of
        """

        current_loc = self.dirlong
        subdir = [f for f in os.listdir(current_loc) if os.path.isdir("{}/{}".format(current_loc, f))]
        sample_dir = {}

        if len(subdir) != 0:
            for sd in subdir:
                for d in os.listdir("{}/{}".format(current_loc, sd)):
                    if "fasta" in d:
                        d_name = d.split("/")[-1].split("_")[1]
                        sample_dir[d_name] = AssemblySample(self.args)
                        sample_dir[d_name].fasta = "{}/{}/{}".format(current_loc, sd, d)
                        sample_dir[d_name].longread = "{}/{}".format(current_loc, sd)

        else:
            for f in os.listdir(current_loc):
                f_name = f.split("/")[-1].split("_")[-1]
                sample_dir[f_name] = AssemblySample(self.args)
                sample_dir[f_name].fasta = "{}{}/{}.fasta".format(current_loc, f, f)
                sample_dir[f_name].longread = current_loc

        if self.dirshort:
            sample_dir = AssemblyCall.samples_to_run(self.dirshort, sample_dir)

        return sample_dir

    def assembly(self):

        """
        Calls for the assemble of each sample object stored in the dictionary
        :return:
        """
        dictionary = AssemblyCall.fasta_finder(self)
        for d in dictionary.values():
            d.assembly_run()

    def fasta_files(self):

        """
        Converts long read file (fastq.gz) to a fasta file
        :return: None
        """
        fastq_dir = FastqFasta(self.args)
        fastq_dir.fastq_convert()


class FastqFasta:
    def __init__(self, arg):
        self.path = arg.long  # location of long reads, can be in subdir or in dir directly
        self.args = arg

    @staticmethod
    def fastq_to_a(input_path):

        """
        Results in a fastq file being converted to a fasta file. Will create a folder based on the sample ID name,
        and place both the fastq and fasta file in it. The output of shasta and spades will be directed to this dir

        :param input_path: (str) is the path to the long read file
        :returns: None.
        """
        acting_dir = '/'.join(input_path.split('/')[:-1])
        fastq_name = (input_path.split(".fastq")[0]).split('/')[-1]
        fasta_name = "{}.fasta".format(fastq_name)

        if not os.path.isdir("{}/{}".format(acting_dir, fastq_name)):
            os.mkdir("{}/{}".format(acting_dir, fastq_name))

        counter = 0
        with open(fasta_name, "w") as ftw:
            with gzip.open(input_path, "rb") as ftr:
                for line in ftr:
                    if counter == 0:
                        ftw.write(line.decode().replace("@", ">"))
                        counter += 1
                    elif counter == 1:
                        ftw.write(line.decode())
                        counter += 1
                    elif counter == 2:
                        counter += 1
                    elif counter == 3:
                        counter = 0

        shutil.move("{}/{}.fastq.gz".format(acting_dir, fastq_name), "{}/{}/".format(acting_dir, fastq_name)) # should not happen
        shutil.move("{}".format(fasta_name), "{}/{}".format(acting_dir, fastq_name)) # should be in samples

    def fastq_convert(self):

        """
        Calls fastq_to_a on fastq.gz files in the directory, or within the subdirectories of the input directory
        :return: None
        """
        parent_path = self.path
        subdir = [sd for sd in os.listdir(parent_path) if os.path.isdir("{}/{}".format(os.getcwd(), sd))]

        if len(subdir) != 0:
            for sd in subdir:
                # print(sd)
                files = [f for f in os.listdir("{}/{}".format(parent_path, sd)) if
                         os.path.isfile("{}/{}/{}".format(parent_path, sd, f))]
                for f in files:
                    if "fastq" in f or "fq" in f:
                        if "filtering" in f:
                            FastqFasta.fastq_to_a("{}/{}/{}".format(parent_path, sd, f))
                        else:
                            FastqFasta.fastq_to_a("{}/{}/{}".format(parent_path, sd, f))
        else:
            files = [f for f in os.listdir(parent_path) if os.path.isfile("{}/{}".format(parent_path, f))]
            for f in files:

                if "fastq" in f or "fq" in f:
                    if "filtering" in f:
                        FastqFasta.fastq_to_a("{}/{}".format(parent_path, f))
                    else:
                        FastqFasta.fastq_to_a("{}/{}".format(parent_path, f))


