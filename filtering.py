import subprocess
import os


class FilteringInput:
    """
    Creates samples based of off input directories
    """

    def __init__(self, args):
        self.samples = {}  # A dictionary of samples
        self.args = args

        if args.LR:  # Directory with longreads
            self.londir = args.LR
        else:
            self.londir = None

        if args.SR:  # Directory with shortreads
            self.shortdir = args.SR
        else:
            self.shortdir = None

        if args.genomesize:  # For filtlong
            self.genomesize = args.genomesize
        else:
            self.genomesize = None

        if args.minlen:  # for filtlong
            self.minlen = args.minlen
        else:
            self.minlen = None

    def samples_sorting(self):
        """
        Using the pathfiles for the longread directory, the shortread directory, or both, initiates a Sample class
        object, and updates parameters. All samples are stored in a dictionary as values, with the corresponding
        key being their sample name
        :return: Dict {Sample Name: <Sample Object>}
        """

        if self.londir and self.shortdir:
            # Filtering will only be done with short reads that match samples for long reads
            longdir = [f for f in os.listdir(self.londir) if os.path.isfile("{}/{}".format(self.londir, f))]
            shortdir = [f for f in os.listdir(self.shortdir) if os.path.isfile("{}/{}".format(self.shortdir,
                                                                                              f))]

            # if '/' not in the end of the shortdir, add it!
            samples_long = list(set(map(lambda x: x.split("_")[0], longdir)))
            for f in samples_long:
                self.samples[f] = FilteringSample(f)
                self.samples[f].minlen = self.minlen
                for r in shortdir:
                    if f in r and "R1" in r:
                        self.samples[f].sr1 = "{}/{}".format(self.shortdir, r)
                    if f in r and "R2" in r:
                        self.samples[f].sr2 = "{}/{}".format(self.shortdir, r)

                for r in longdir:
                    if f in r:
                        self.samples[f].longread = "{}/{}".format(self.londir, r)
                        # add minlen

        elif self.londir:

            longdir = [f for f in os.listdir(self.londir) if os.path.isfile("{}/{}".format(self.londir, f))]
            # if '/' not in the end of the shortdir, add it!
            for f in longdir:
                name = f.split("/")[-1].split("_")[0]
                self.samples[name] = FilteringSample(name)
                self.samples[name].longread = "{}/{}".format(self.londir, f)
                self.samples[name].minlen = self.minlen

        elif self.shortdir:
            shortdir = [f for f in os.listdir(self.shortdir) if os.path.isfile("{}{}".format(self.shortdir, f))]
            # if '/' not in the end of the shortdir, add it!

            names = set(map(lambda x: x.split("_")[0], shortdir))
            for n in names:
                self.samples[n] = FilteringSample(n)
                for r in shortdir:
                    if "R1" in r and n in r:
                        self.samples[n].sr1 = "{}/{}".format(self.shortdir, r)
                    elif "R2" in r and n in r:

                        self.samples[n].sr2 = "{}/{}".format(self.shortdir, r)

        return self.samples

    def filtering(self):
        """
        Calls what kind of read filtering should be used for each Sample object. Call appropriate Sample method.
        :return: None
        """
        dictionary = self.samples_sorting()

        if self.shortdir and self.londir:
            for v in dictionary.values():
                if self.shortdir:
                    # Filter long reads with trimmed short reads as external reference
                    v.genomesize = self.genomesize
                    v.bbduk()
                    v.filtlong()
                else:
                    # Filter long reads without the external reference
                    v.genomesize = self.genomesize
                    v.filtlong()

        elif self.shortdir:
            for v in dictionary.values():
                v.bbduk()

        elif self.londir:
            for v in dictionary.values():
                v.genomesize = self.genomesize
                v.filtlong()


class FilteringSample:
    """
    Represents each sample that has longread/short read data collected
    """

    def __init__(self, name):
        self.name = name
        self.longread = None
        self.sr1 = None
        self.sr2 = None
        self.genomesize = None
        self.minlen = None
        self.chop = None

    def filtlong(self):
        """
        Results in the filtering of longreads in the path, self.longread, and it can be with external references
        or not. Writes a fastq.gz file.
        :return: None
        """
        output_dir = '/'.join(self.longread.split('/')[:-1])
        self.chop = "{}/pore_chopped_{}".format(output_dir, self.name)
        cmd_chop = ["porechop", "-i", self.longread, "-o", self.chop]
        filter_cmd = ["filtlong",
                      "--min_length", self.minlen,
                      "--keep_percent", "90"]

        if self.sr1 and self.sr2:
            filter_cmd.extend(["-1", self.sr1, "-2", self.sr2])
        if self.genomesize:
            filter_cmd.extend(["--target_bases", str(int(self.genomesize) * 100)])
        filter_cmd.append(self.chop)

        subprocess.run(cmd_chop)
        filt = subprocess.Popen(filter_cmd, stdout=subprocess.PIPE)
        gz = subprocess.Popen(["gzip"], stdin=filt.stdout,
                              stdout=subprocess.PIPE)
        output = gz.communicate()[0]
        with open("{}/filtered_{}".format(output_dir, self.name), 'wb') as f:
            f.write(output)
        # os.remove(self.chop)

    def bbduk(self):
        """
        Results in the filtering of shortreads in the path. Trimmed reads are saved in a directory, trimmed_reads,
        within the directory the illumina short reads are kept.
        :return: None
        """
        out_path = '/'.join(self.sr1.split('/')[:-1])
        sr1_name = self.sr1.split('/')[-1]
        sr2_name = self.sr2.split('/')[-1]

        cmd = ["/home/bioinfo/bbmap/bbduk.sh",
               "in={}".format(self.sr1), "in2={}".format(self.sr2),
               "ref=/home/bioinfo/bbmap/resources/adapters.fa", "ktrim=r", "k=23", "mink=11", "qtrim=lr", "trimq=10",
               "hdist=1", "tbo", "tpe",
               "out={}/trimmed_reads/trimmed_{}".format(out_path, sr1_name),
               "out2={}/trimmed_reads/trimmed_{}".format(out_path, sr2_name)]

        self.sr1 = "{}/trimmed_reads/trimmed_{}".format(out_path, sr1_name)
        self.sr2 = "{}/trimmed_reads/trimmed_{}".format(out_path, sr2_name)

        subprocess.run(cmd)

