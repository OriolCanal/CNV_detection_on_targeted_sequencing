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
        self.gatk_dirname = "GATK_gCNV"
        self.gatk_image = docker_conf.gatk["image"]
        self.gatk_version = docker_conf.gatk["version"]
        self.reference_fasta = reference_conf.hg19.fasta_path
        self.gatk_folder = os.path.join(reference_conf.main_dir, "GATK")
        self.mappability_folder = os.path.join(self.gatk_folder, "mappability_track")
        self.mappability_track_path = os.path.join(self.mappability_folder, "k50.umap.bed")
        self.gatk_volume = "/gatk_vol"
        self.contig_ploidy = os.path.join(self.gatk_folder, "contig_ploidy_prior.tsv")

    def run_preprocess_intervals(self, Bed):
        
        fasta_dir = os.path.dirname(self.reference_fasta)
        fasta_filename = os.path.basename(self.reference_fasta)

        fasta_volume = "/fasta_dir"

        preprocessed_interval_list_filename = f"preprocessed_{Bed.interval_list_filename}"
        preprocessed_interval_list = os.path.join(Bed.dir, preprocessed_interval_list_filename)

        if os.path.exists(preprocessed_interval_list):
            logger.info(f"{preprocessed_interval_list} already exists, Gatk PreprocessIntervals won't be run")
            Bed.set_preprocessed_intervals_list_path(preprocessed_interval_list)
        cmd = [
            self.docker_path, "run",
            "-v", f"{fasta_dir}:{fasta_volume}",
            "-v", f"{Bed.dir}:{Bed.volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "PreprocessIntervals",
            "-R", f"{fasta_volume}/{fasta_filename}",
            "-L", f"{Bed.volume}/{Bed.interval_list_filename}",
            "--bin-length", "0",
            "-imr", "OVERLAPPING_ONLY",
            "-O", f"{Bed.volume}/{preprocessed_interval_list_filename}"
        ]
        cmd_str = " ".join(cmd)
        logger.info(f"Running GATK Preprocess Intervals:\n{cmd_str}")
        subprocess.run(cmd, encoding="utf-8", capture_output=True)
        
        Bed.set_preprocessed_intervals_list_path(preprocessed_interval_list)
        return (preprocessed_interval_list)
    
    def create_gatk_folder(self, Bam):
        self.runs_dir = os.path.dirname(Bam.dir)
        self.gatk_dir = os.path.join(self.runs_dir, self.gatk_dirname)

        if not os.path.isdir(self.gatk_dir):
            logger.info(
                f"Creating new folder to store GATK gCNV files: {self.gatk_dir}"
            )
            os.mkdir(self.gatk_dir)
        return(self.gatk_dir)
    
    
    def run_collect_read_counts(self, Bam, Bed):
        """
        Run GATK CollectReadCounts.
        """

        self.create_gatk_folder(Bam)

        if not Bed.preprocessed_intervals_path:
            raise ValueError(
                f"You need to preprocess the interval list file first! Run the funciton Gatk_gCNV.run_preprocess_intervals()"
            )


        fasta_dir = os.path.dirname(self.reference_fasta)
        fasta_filename = os.path.basename(self.reference_fasta)


        hdf5_filename = f"{Bam.filename}.hdf5"
        hdf5_path = os.path.join(self.gatk_dir, hdf5_filename)
        
        if os.path.exists(hdf5_path):
            logger.info(
                f"GATK CollectReadCounts won't be executed as {hdf5_path} file already exists."
            )
            Bam.set_hdf5_path_and_filename(hdf5_path)
            return(Bam.hdf5_path)

        fasta_volume = "/fasta_dir"

        cmd = [
            self.docker_path, "run",
            "-v", f"{self.runs_dir}:/runs_dir",
            "-v", f"{fasta_dir}:{fasta_volume}",
            "-v", f"{Bed.dir}:{Bed.volume}",
            "-v", f"{Bam.dir}:{Bam.volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "CollectReadCounts",
            "-L", f"{Bed.volume}/{Bed.preprocessed_intervals_filename}",
            "-R", f"{fasta_volume}/{fasta_filename}",
            "-imr", "OVERLAPPING_ONLY",
            "-I", f"{Bam.volume}/{Bam.filename}",
            "-O", f"/runs_dir/{self.gatk_dirname}/{hdf5_filename}"
        ]

        cmd_str = " ".join(cmd)
        logger.info(f"Running GATK CollectReadCounts:\n{cmd_str}")
        subprocess.run(cmd, encoding="utf-8", capture_output=True)
        Bam.set_hdf5_path_and_filename(hdf5_path)
        return(Bam.hdf5_path)
    
    def run_index_feature_file(self, Bed):
        """
        Indexing ./bed/mappability_track/k36.umap.bed.gz required to run annotate intervals
        """

        mappability_folder = os.path.basename(self.mappability_folder)
        mappability_filename = os.path.basename(self.mappability_track_path)

        index_map_track_path = f"{self.mappability_track_path}.idx"
        if os.path.isfile(index_map_track_path):
            return(index_map_track_path)
        cmd = [
            self.docker_path, "run",
            "-v", f"{Bed.dir}:{Bed.volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "IndexFeatureFile",
            "-I", f"{Bed.volume}/{mappability_folder}/{mappability_filename}"
        ]
        cmd_str = " ".join(cmd)
        logger.info(
            f"Indexing mappability track:\n{cmd_str}"
        )
        subprocess.run(cmd, encoding="utf-8", capture_output=True)
        return(index_map_track_path)


    def run_annotate_intervals(self, Bed):
        """
        Exclude problematic reagions based on GC content and mappability of intervals
        """

        if not Bed.preprocessed_intervals_path:
            raise ValueError(
                f"You need to preprocess the interval list file first! Run the funciton Gatk_gCNV.run_preprocess_intervals()"
            )
        
        fasta_dir = os.path.dirname(self.reference_fasta)
        fasta_filename = os.path.basename(self.reference_fasta)


        mappability_folder = os.path.basename(self.mappability_folder)
        mappability_filename = os.path.basename(self.mappability_track_path)


        gc_annotated_bed_filename = f"annotated_{Bed.preprocessed_intervals_filename}"
        gc_annotated_bed_path = os.path.join(Bed.dir, gc_annotated_bed_filename)
        if os.path.exists(gc_annotated_bed_path):
            logger.info(
                f"GATK Annotate intevals won't be executed as {gc_annotated_bed_path} file already exists."
            )
            Bed.set_annotated_intervals_path(gc_annotated_bed_path)
            return(gc_annotated_bed_path)
        
        fasta_volume = "/fasta_dir"

        cmd = [
            self.docker_path, "run",
            "-v", f"{fasta_dir}:{fasta_volume}",
            "-v", f"{Bed.dir}:{Bed.volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "AnnotateIntervals",
            "-L", f"{Bed.volume}/{Bed.preprocessed_intervals_filename}",
            "-R", f"{fasta_volume}/{fasta_filename}",
            "-imr", "OVERLAPPING_ONLY",
            "--mappability-track", f"{Bed.volume}/{mappability_folder}/{mappability_filename}",
            "-O", f"{Bed.volume}/{gc_annotated_bed_filename}"
        ]
    
        cmd_str = " ".join(cmd)
        logger.info(f"Running GATK AnnotateIntervals:\n{cmd_str}")

        subprocess.run(cmd, encoding="utf-8", capture_output=True)
        Bed.set_annotated_intervals_path(gc_annotated_bed_path)

        return(gc_annotated_bed_path)
    
    def get_input_read_count_files(self):
        """
        Get a list of input read counts files found in runs/GATK-gCNV7 to be given to gatk docker.
        E.g. ["-I", "gatk_vol/RB35645.hdf5", "-I", "gatk_vol/RB34532.hdf5]"""
        read_counts_input = list()
        sample_hdf5_counts_paths = os.listdir(self.gatk_dir)
        sample_hdf5_counts_filenames = [os.path.basename(counts_filename) for counts_filename in sample_hdf5_counts_paths]
        for sample_hdf5_filename in sample_hdf5_counts_filenames:
            input_files_cmd = ["-I", f"{self.gatk_volume}/{sample_hdf5_filename}"]
            read_counts_input.extend(input_files_cmd)
        
        return(read_counts_input)
    def run_filter_intervals(self, Bed):
        """Given specific intervals (annotated intervals), and counts output by CollectReadCoutns,
        outputs a filtered Picard interval list"""

        sample_hdf5_counts_paths = os.listdir(self.gatk_dir)
        sample_hdf5_counts_filenames = [os.path.basename(counts_filename) for counts_filename in sample_hdf5_counts_paths]

        filtered_intervals_filename = f"filtered_{Bed.preprocessed_intervals_filename}"
        filtered_intervals_path = os.path.join(Bed.dir, filtered_intervals_filename)

        if os.path.isfile(filtered_intervals_path):
            logger.info(
                f"Filtered intervals file already exists: {filtered_intervals_path}, GATK FilteredIntervals won't be run"
            )
            Bed.set_filtered_intervals_path(filtered_intervals_path)
            return(filtered_intervals_path)

        cmd = [
            self.docker_path, "run", "-u", "$(id -u):$(id -g)",
            "-v", f"{Bed.dir}:{Bed.volume}",
            "-v", f"{self.gatk_dir}:{self.gatk_volume}",
            f"{self.gatk_image}:{self.gatk_version}",
            "gatk", "FilterIntervals",
            "-L", f"{Bed.volume}/{Bed.preprocessed_intervals_filename}",
            "--annotated-intervals", f"{Bed.volume}/{Bed.annotated_intevals_filename}",
        ]

        read_counts_cmd = self.get_input_read_count_files()
        
        cmd.extend(read_counts_cmd)

        cmd2 = [
            "-imr", "OVERLAPPING_ONLY",
            "-O", f"{Bed.volume}/{filtered_intervals_filename}",
        ]

        cmd.extend(cmd2)

        cmd_str = " ".join(cmd)
        logger.info(
            f"Running GATK FilterIntervals:\n {cmd_str}"
        )
        subprocess.run(cmd, encoding="utf-8", capture_output=True)
        Bed.set_filtered_intervals_path(filtered_intervals_path)

        return(filtered_intervals_path)
    
    def run_determine_germline_contig_ploidy(self, Bed, Bam):
        cmd = [
            self.docker_path, "run", "-u", "$(id -u):$(id -g)",
            "-v", f"{Bed.dir}:{Bed.volume}",
            "-v", f"{self.gatk_dir}:{self.gatk_volume}",
            f"{self.gatk_image}:{self.gatk_version}",
        ]