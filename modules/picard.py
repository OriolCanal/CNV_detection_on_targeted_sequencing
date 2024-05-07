import subprocess
import os
import pandas as pd
from modules.log import logger

class Picard():
    def __init__(self, docker_conf, ref_conf):
        self.docker_path = "/usr/bin/docker"
        self.docker_conf = docker_conf    
        self.bed = ref_conf.bed
        self.bed_dir = os.path.dirname(self.bed)
        self.bed_filename = os.path.basename(self.bed)
        self.reference_dir = ref_conf.hg19.dir_path
        self.ref_dict_file = ref_conf.hg19.dict
        self.ref_fasta_file = ref_conf.hg19.fasta

    def create_fasta_dict(self):
        self.ref_dict_file = self.ref_fasta_file.replace("fasta", "dict")
        picard_cmd = [
            "docker", "run",
            "-v", f"{self.bed_dir}:/bed_dir",
            "-v", f"{self.reference_dir}:/ref_dir",
            self.docker_conf.picard["image"],
            "java", "-Xmx60g", "-jar",
            "/usr/picard/picard.jar",
            "CreateSequenceDictionary",
            "-R", f"/ref_dir/{self.ref_fasta_file}",
            "-O", f"/ref_dir/{self.ref_dict_file}"
        ]

        cmd_str = " ".join(picard_cmd)
        logger.info(
            f"Creating dict file from fasta, if already exists it won't be created:\n{cmd_str}"
        )
        subprocess.run(picard_cmd, encoding="utf-8", capture_output=True)
        return(self.ref_dict_file)

    def run_bed_to_interval_list(self, Bed):
        
        Bed.interval_list_file = f"{Bed.filename}.interval_list"
        Bed.interval_list_path = os.path.join(Bed.dir, Bed.interval_list_file)
        if os.path.isfile(Bed.output_list_path):
            logger.info(
                f"Listfile for bed {Bed.filename} already created in {Bed.interval_list_path}"
            )
            return(Bed.interval_list_path)
        picard_cmd = [
            "docker", "run",
            "-v", f"{Bed.bed_dir}:/bed_dir",
            "-v", f"{self.reference_dir}:/ref_dir",
            self.docker_conf.picard["image"],
            "java", "-Xmx60g", "-jar",
            "/usr/picard/picard.jar",
            "BedToIntervalList",
            "-I", f"/bed_dir/{Bed.filename}",
            "-O", f"/bed_dir/{self.interval_list_file}",
            "-SD", f"/ref_dir/{self.ref_dict_file}"
        ]

        command_str = " ".join(picard_cmd)
        logger.info(
            f"Creating picard interval list:\n{command_str}"
        )
        subprocess.run(picard_cmd, encoding="utf-8", capture_output=True)
        return(self.output_list_path)
    
    def run_collectHsMetrics(self, bam_file):
        if not hasattr(self, "output_list_path"):
            self.run_bed_to_interval_list()

        bam_filename = os.path.basename(bam_file)
        self.hsmetrics_filename = f"{bam_filename}_hs_metrics.txt"
        picard_dirname = "Picard"
        picard_dir = os.path.join(os.path.dirname(bam_file), picard_dirname)
        self.hsmetrics_path = os.path.join(picard_dir, self.hsmetrics_filename)
        if os.path.isfile(self.hsmetrics_path):
            logger.info(
                f"HsMetrics output file already available: {self.hsmetrics_path}"
            )
            return(self.hsmetrics_path)

        if not os.path.isdir(picard_dir):
            os.mkdir(picard_dir)
        picard_cmd = [
            "docker", "run",
            "-v", f"{self.bed_dir}:/bed_dir",
            "-v", f"{self.reference_dir}:/ref_dir",
            "-v", f"{os.path.dirname(bam_file)}:/bam_dir",
            self.docker_conf.picard["image"],
            "java", "-Xmx60g", "-jar",
            "/usr/picard/picard.jar",
            "CollectHsMetrics",
            "--INPUT", f"/bam_dir/{os.path.basename(bam_file)}",
            "--OUTPUT", f"/bam_dir/{picard_dirname}/{self.hsmetrics_filename}",
            "--BAIT_INTERVALS", f"/bed_dir/{self.output_list_file}",
            "--TARGET_INTERVALS", f"/bed_dir/{self.output_list_file}"
        ]

        command_str = " ".join(picard_cmd)
        logger.info(
            f"Running picard CollectHsMetrics:\n{command_str}"
        )
        subprocess.run(picard_cmd, encoding="utf-8", capture_output=True)
        return(self.hsmetrics_path)

    def get_picard_metrics(self):

        with open(self.hsmetrics_path) as f:
            for line in f:
                if not line.startswith("#"):
                    continue
                if line.startswith("## METRICS CLASS"):
                    metrics_info_line = next(f)
                    
                    metrics_value_line = next(f)
                
                    return(metrics_info_line, metrics_value_line)


class Metrics_Df():
    def __init__(self):
        self.metrics_df = pd.DataFrame()
        self.df_has_header = False

    def has_df_header(self):
        return not self.metrics_df.columns.empty

    def add_metrics_header(self, line):
        header = line.strip().split("\t")
        # this items are not specified and should be removed
        items_to_remove = ["LIBRARY", "SAMPLE", "READ_GROUP"]
        for item in items_to_remove:
            if item in header:
                header.remove(item)
        
            
        self.metrics_df = pd.DataFrame(columns=header)  # Initialize DataFrame with header
        self.df_has_header = True

    def add_metrics_line(self, values_line):
        metrics_values = values_line.strip().split("\t")

        if len(metrics_values) != len(self.metrics_df.columns):
            raise ValueError("Length of metrics values does not match the length of columns.")
        
        self.metrics_df = self.metrics_df._append(dict(zip(self.metrics_df.columns, metrics_values)), ignore_index=True)
    