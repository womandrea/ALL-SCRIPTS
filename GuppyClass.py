import os, subprocess, shutil, argparse
import pandas as pd


class GuppyBaseCallerRun:

    @staticmethod
    def model_file(kit, flowcell):
        cmd = ["guppy_basecaller", "--print_workflows"]
        p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        stdout, stderr = p1.communicate()

        with open("workflows.txt", "wb") as ftw:
            ftw.write(stdout)

        with open("workflows.txt", "r") as ftr:
            lines = ftr.readlines()
            for l in lines:
                try:
                    k = l.split(" ")[1]
                    f = l.split(" ")[0]
                    if kit == k and flowcell == f:
                        model_file = "/home/bioinfo/prog/ont-guppy/data/template" + (l.split(' ')[-1]).split()[0].split('dna')[-1] + ".jsn"
                except IndexError:
                    pass
        return model_file

    def __init__(self, args):
        """
        Initialize the GuppyBaseCallerRun Object.
        """
        self.args = args
        self.output = args.fastq 
        self.input = args.fast5  
        self.barcodekit = args.barcodekit
        self.kit = args.kit
        self.flowcell = args.flowcell
        self.model_file = GuppyBaseCallerRun.model_file(self.kit, self.flowcell)
        print("Model file {} is being used.\n".format(self.model_file))
        if int(args.threads) > len(os.sched_getaffinity(0)):
            self.threads = str(len(os.sched_getaffinity(0)))
            print("You are trying to run this process on more cores than available. We changed it to {} cores.\n".format(
                len(os.sched_getaffinity(0))))
        else:
            self.threads = args.threads

    def guppy(self):
        """
        Calls guppy_basecaller function on run information.


        :return: None.
        """
        cmd = ["guppy_basecaller",
               "-i", self.input,  # MANDATORY INPUT
               "-s", self.output,  # MANDATORY INPUT
               "--flowcell", self.flowcell,  # MANDATORY INPUT
               "--kit", self.kit,  # MANDATORY INPUT
               "--records_per_fastq", self.args.rpfq,
               "--num_barcode_threads", str(self.threads),
               '--calib_reference', self.args.calib_ref,
               "--model_file", self.model_file,  # ADJUST
               "--hp_correct", self.args.hp,
               "--num_caller", self.args.num_caller,
               "--gpu_runners_per_device", self.args.grpd,
               "--chunk_size", self.args.chunksize,
               "--chunks_per_runner", self.args.chunkprun,
               "--device", self.args.device]

        if self.args.barcodekit != "None":
            cmd.extend(['--barcode_kits', self.barcodekit])
        if self.args.qsf == "Y":
            cmd.extend(["--qscore_filtering"])
        if self.args.compress == "Y":
            cmd.extend(["--compress_fastq"])
        if self.args.recursive == "Y":
            cmd.extend(["--recursive"])
        if self.args.trimbarcodes == "Y":
            cmd.extend(["--trim_barcode"])
        if self.args.calib_detect == "Y":
            cmd.extend(["--calib_detect"])
            
        subprocess.run(cmd)

    def concatenate_output(self):
        """
        Meant to be called after base calling takes place. Calls appropriate methods for renaming of
        barcode directories, concatenation of base calling output, and removal of unnecessary files.

        :return: None
        """
        directory_contents = os.listdir(self.output)
        for i in directory_contents:
            if os.path.isdir("{}/{}".format(self.output, i)) and (i == 'pass' or i == 'fail'):
                action_directory = FileTreatment("{}/{}".format(self.output, i))
                if self.args.ros != "None":
                    action_directory.csv = self.args.ros
                    action_directory.file_renaming()

                action_directory.concatenation()


class FileTreatment:
    def __init__(self, directory):
        self.csv = None
        self.directory = directory

    def file_renaming(self):
        """
        If a csv is provided in the input, resulting in renaming subdirectories from barcode to corresponding sample
        name.

        :return: None.
        """
        working_dir = self.directory
        df = pd.read_csv(self.csv, header=None)
        column_names = df.columns.tolist()
        barcode_names = list(set(df[column_names[0]].tolist()))
        replacement_names = list(set(df[column_names[1]].tolist()))

        df_dict = {}
        if len(replacement_names) != len(barcode_names):
            print("It doesn't look like you have enough unique names or identifiers for your runs. Try again, \
                after making all of your rows in column 2 unique")
        else:
            i = 0
            while i < len(barcode_names):
                df_dict[barcode_names[i]] = replacement_names[i]
                i += 1

        print(working_dir)
        subdirectory_list = [f for f in os.listdir(working_dir) if os.path.isdir('{}/{}'.format(working_dir, f))]
        for sd in subdirectory_list:
            if sd.lower() in df_dict.keys():
                os.rename("{}/{}".format(working_dir, sd), "{}/{}".format(working_dir, df_dict[sd.lower()]))
            else:
                print("{} could not be renamed. No replacement found in the csv file provided.".format(sd))

    def concatenation(self):
        """
        If subdirectories are present, concatenated all "fastq_runid" fastq.gz files. Otherwise, concatenates all
        files starting with "fastq_runid" within the directory. Calls for the deletion of copied files.

        EDIT: Make for both non-gzipped and gzipped

         with gzip.open(fasta, 'rt') if fasta.endswith('.gz') else open(fasta, 'r') as f:

        :return: None.
        """
        starting_directory = self.directory
        contents = os.listdir(starting_directory)
        directory_name = self.directory.split('/')[-1] # provides fail or pass
        subdirectories = [c for c in contents if os.path.isdir('{}/{}'.format(starting_directory, c))]
        if len(subdirectories) != 0:
            for sd in subdirectories:
                ref_dir = '{}/{}'.format(starting_directory, sd)
                files = [f for f in os.listdir(ref_dir) if os.path.isfile('{}/{}'.format(ref_dir, f))]
                text = "{}/{}/{}_{}.fastq.gz".format(starting_directory, sd, directory_name, sd)
                with open(text, "wb") as ftw:
                    for f in files:
                        if f.startswith("fastq_runid"):
                            with open("{}/{}/{}".format(starting_directory, sd, f), "rb") as ftr:
                                shutil.copyfileobj(ftr, ftw)

                # checks to see the size of the final concatenated file is greater than zero.
                if os.stat(text).st_size != 0:
                    FileTreatment.clean_up("{}/{}".format(starting_directory, sd))
                else:
                    print("The concatenated file seems to be empty. We won't delete any of the files in this folder.")

        else:
            text = "{}/{}.fastq.gz".format(starting_directory, directory_name)
            files = [f for f in os.listdir(starting_directory) if os.path.isfile('{}/{}'.format(starting_directory, f))]
            with open(text, "wb") as ftw:
                for f in files:
                    if f.startswith("fastq_runid"):
                        with open("{}/{}".format(starting_directory, f), "rb") as ftr:
                            shutil.copyfileobj(ftr, ftw)

            # checks to see the size of the final concatenated file is greater than zero.
            if os.stat(text).st_size != 0:
                FileTreatment.clean_up(starting_directory)
            else:
                print("The concatenated file seems to be empty. We won't delete any of the files in this folder.")

    @staticmethod
    def clean_up(directory):
        """
        Removes all files that have been copied by concatenation.

        :return: None.
        """
        for f in os.listdir(directory):
            if f.startswith("fastq_runid"):
                os.remove("{}/{}".format(directory, f))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="To be used for guppybase calling. Provide full file paths for all the inputs.")

    parser.add_argument("-fast5", metavar='PATH', type=str,
                        help="Input directory for guppy_basecaller.", dest='fast5', required=True)

    parser.add_argument("-output", metavar='PATH', type=str,
                        help="Directory path for guppy_basecaller output.", dest='fastq', required=True)

    parser.add_argument("-t", metavar="INT", type=str, help="Number of threads to be used for barcoding. Default 1.",
                        dest='threads', default="1")

    parser.add_argument("-flowcell", metavar='F', type=str, required=True
                        help="Flowcell ID of run. Required, write as a string in single quotations", dest='flowcell')

    parser.add_argument("-kit", metavar='K', type=str, required=True
                        help="Kit ID. Required, write as a string in single quotations", dest='kit')

    parser.add_argument("-barcodekit", metavar="KIT", dest="barcodekit", type=str, default="None", required=True
                        help="Barcode appropriate for used kit. Default is None. If there is a barcode input, write it within single quotation marks")

    parser.add_argument("-chunk_size", metavar="INT", type=str, help="Default is 1000", default="1000", dest="chunksize")  # chunk size

    parser.add_argument("-chunk_per_run", metavar="INT", type=str, help="Default is 1000", default="1000", dest="chunkprun")  # chunk per runner

    parser.add_argument("-rosetta", metavar="PATH", dest="ros", default= "None",
                        help="A .csv file you can provide to rename barcodes.\
                        Requires two columns, the first containing barcodes,\
                        and the second adjacent should be the intended name for that assigned barcode. No header should\
                        be provided in this csv.")  #

    parser.add_argument("-gpu_rpd", type=str, help="GPU runner per device. Default 2.", dest="grpd", default="2")

    parser.add_argument("-records_per_fastq", help="Default 0.", type=str, default="0", dest="rpfq", metavar="0")

    parser.add_argument("-recursive", metavar="Y/N", default="Y", dest="recursive", help="Signify yes (Y) or no (N). Default Y.")

    parser.add_argument("-trim_barcodes", metavar="Y/N", default="Y", dest="trimbarcodes", help="Signify yes (Y) or no (N). Default Y.")

    parser.add_argument("-calib_detect", metavar="Y/N", default="Y", dest="calib_detect", help="Signify yes (Y) or no (N). Default Y.")

    parser.add_argument("-compress", metavar="Y/N", default="Y", dest="compress", help="Signify yes (Y) or no (N). Default Y.")

    parser.add_argument("-calib-ref", dest="calib_ref", metavar="arg", default="lambda_3.6kb.fasta")

    parser.add_argument("-hp_correct", dest="hp", metavar="arg", type=str, help="Default 1", default="1")

    parser.add_argument("-qscore_filt", dest="qsf", help="Signify yes (Y) or no (N). Default Y.", default="Y")

    parser.add_argument("-num_caller", metavar="arg", dest="num_caller",help="Default 4.", default="4")

    parser.add_argument("-device", metavar="arg", dest="device", help="Default cudo:0", default="cuda:0")

    args = parser.parse_args()

    CurrentRun = GuppyBaseCallerRun(args)
    CurrentRun.guppy()  # results in base calling
    CurrentRun.concatenate_output()  # results in files being renamed, concatenated, and deleted if redundant
