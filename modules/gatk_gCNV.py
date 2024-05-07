import os
import subprocess

from modules.log import logger

class Gatk_gCNV():
    """
    Run all the steps of GATK gCNV in:
    https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-rare-germline-copy-number-variants
    """
    def __init__(self, docker_conf, reference_conf, Bed):

        self.docker_path = "/usr/bin/docker"
        self.Bed = Bed
        self.bed = docker_conf.bed
        self.bed_dir = os.path.dirname(self.bed)
        self.interval_list_path = f"{self.bed}.interval_list"
        self.interval_list_filename = os.path.basename(self.interval_list_path)
        self.gatk_dirname = "GATK_gCNV"
        self.gc_annotated_bed_path = f"{self.bed}.interval_list.annotated.tsv"
        self.gatk_image = docker_conf.gatk["image"]
        self.gatk_version = docker_conf.gatk["version"]
        self.reference_fasta = reference_conf.hg19.fasta_path
        self.preprocessed_interval_list = None
        self.mappability_folder = os.path.join(self.bed_dir, "mappability_track")
        self.mappability_track_path = os.path.join(self.mappability_folder, "k36.umap.bed.gz")

    def get_intervals_from_bed(self):
        pass

    def run_preprocess_intervals(self, Bed):
        
        fasta_dir = os.path.dirname(self.reference_fasta)
        fasta_filename = os.path.basename(self.reference_fasta)


        bed_volume = "/bed_dir"
        fasta_volume = "/fasta_dir"

        output_preprocessed_interval_list = f"preprocessed_{self.interval_list_filename}"
        self.preprocessed_interval_list = os.path.join(Bed.dir, output_preprocessed_interval_list)

        if os.path.exists(self.preprocessed_interval_list):
            logger.info(f"{self.preprocessed_interval_list} already exists, Gatk PreprocessIntervals won't be run")

        cmd = [
            self.docker_path, "run",
            "-v", f"{fasta_dir}:{fasta_volume}",
            "-v", f"{Bed.dir}:{Bed.volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "PreprocessIntervals",
            "-R", f"{fasta_volume}/{fasta_filename}",
            "-L", f"{Bedvolume}/{Bed.filename}",
            "--bin-length", "0",
            "-imr", "OVERLAPPING_ONLY",
            "-O", f"{Bed.volume}/{output_preprocessed_interval_list}"
        ]
        cmd_str = " ".join(cmd)
        logger.info(f"Running GATK Preprocess Intervals:\n{cmd_str}")
        subprocess.run(cmd, encoding="utf-8", capture_output=True)

        return(self.preprocessed_interval_list)
    
    def create_gatk_folder(self, Bam):
        
        self.gatk_dir = os.path.join(Bam.path, self.gatk_dirname)

        if not os.path.isdir(self.gatk_dir):
            logger.info(
                f"Creating new folder to store GATK gCNV files: {self.gatk_dir}"
            )
            os.mkdir(self.gatk_dir)
        return(self.gatk_dir)
    
    
    def run_collect_read_counts(self, Bam):
        """
        Run GATK CollectReadCounts.
        """

        self.create_gatk_folder(Bam)

        if not self.preprocessed_interval_list:
            raise ValueError(
                f"You need to preprocess the interval list file first! Run the funciton Gatk_gCNV.run_preprocess_intervals()"
            )


        fasta_dir = os.path.dirname(self.reference_fasta)
        fasta_filename = os.path.basename(self.reference_fasta)

        bed_dir = os.path.dirname(self.bed)

        hdf5_filename = f"{Bam.filename}.hdf5"
        hdf5_path = os.path.join(self.gatk_dir, hdf5_filename)
        Bam.set_hdf5_path_and_filename(hdf5_path)
        if os.path.exists(Bam.hdf5_path):
            logger.info(
                f"GATK CollectReadCounts won't be executed as {Bam.hdf5_path} file already exists."
            )
            return(Bam.hdf5_path)

        bed_volume = "/bed_dir"
        fasta_volume = "/fasta_dir"

        cmd = [
            self.docker_path, "run",
            "-v", f"{fasta_dir}:{fasta_volume}",
            "-v", f"{bed_dir}:{bed_volume}",
            "-v", f"{Bam.dir}:{Bam.volume}"
            f"{self.gatk_image}:{self.gatk_version}",
            "CollectReadCounts",
            "-L", f"{bed_volume}/{self.interval_list_filename}",
            "-R", f"{fasta_volume}/{fasta_filename}",
            "-imr", "OVERLAPPING_ONLY",
            "-I", f"{Bam.volume}/{Bam.filename}"
            "-O", f"{Bam.volume}/{self.gatk_dirname}/{Bam.hdf5_filename}"
        ]

        cmd_str = " ".join(cmd)
        logger.info(f"Running GATK CollectReadCounts:\n{cmd_str}")
        subprocess.run(cmd, encoding="utf-8", capture_output=True)

        return(Bam.hdf5_path)
    
    def run_annotate_intervals(self, Bed):
        """
        Exclude problematic reagions based on GC content and mappability of intervals
        """

        if not Bam.preprocessed_interval_list:
            raise ValueError(
                f"You need to preprocess the interval list file first! Run the funciton Gatk_gCNV.run_preprocess_intervals()"
            )
        
        fasta_dir = os.path.dirname(self.reference_fasta)
        fasta_filename = os.path.basename(self.reference_fasta)

        bed_dir = os.path.dirname(self.bed)

        mappability_folder = os.path.basename(self.mappability_folder)
        mappability_file = os.path.basename(self.mappability_track_path)


        gc_annotated_bed_filename = os.path.basename(self.gc_annotated_bed_path)
        if os.path.exists(self.gc_annotated_bed_path):
            logger.info(
                f"GATK Annotate intevals won't be executed as {self.gc_annotated_bed_path} file already exists."
            )
            return(self.gc_annotated_bed_path)
        
        bed_volume = "/bed_dir"
        fasta_volume = "/fasta_dir"

        cmd = [
            self.docker_path, "run",
            "-v", f"{fasta_dir}:{fasta_volume}",
            "-v", f"{bed_dir}:{bed_volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "AnnotateIntervals",
            "-L", f"{bed_volume}/{self.preprocessed_interval_list}",
            "-R", f"{fasta_volume}/{fasta_filename}",
            "-imr", "OVERLAPPING_ONLY",
            "--mappability-track", f"{bed_volume}/{mappability_folder}/{mappability_file}",
            "-O", f"{bed_volume}/{gc_annotated_bed_filename}"
        ]
    
        cmd_str = " ".join(cmd)
        logger.info(f"Running GATK AnnotateIntervals:\n{cmd_str}")

        subprocess.run(cmd, encoding="utf-8", capture_output=True)

        return(self.gc_annotated_bed_path)