import os

class Bam():
    def __init__(self, bam_path):
        self.path = bam_path
        self.filename = os.path.basename(self.path)
        self.dir = os.path.dirname(self.path)
        self.bai_path = self.path.replace("bam", "bai")
        self.bai_filename = self.filename.replace("bam", "bai")
        self.validate_bam_bai()
        self.bam_volume = "/bam_dir"

        self.hdf5_path = None
        self.hdf5_filename = None

    def validate_bam_bai(self):
        bam_size = os.path.getsize(self.path)
        bai_size = os.path.getsize(self.bai_path)
        if os.path.exists(self.path) and os.path.exists(self.bai_path):
            if bam_size > 0 and bai_size > 0:
                return True
        
        raise ValueError(
            f"Bam and bai files are not present or its size is 0: \n\t{self.path}\n\t{self.bai_path}"
        )
    

    def set_hdf5_path_and_filename(self, hdf5_path):
        self.hdf5_path = hdf5_path
        self.hdf5_filename = os.path.basename(hdf5_path)