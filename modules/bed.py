import os

class Bed():
    def __init__(
            self,
            bed_path,
            interval_list_path=None,
            annotated_intervals=None,
            preprocessed_intervals_path=None
        ):
        self.path = bed_path
        self.filename = os.path.basename(self.bed_path)
        self.dir = os.path.dirname(self.path)
        self.interval_list_path = interval_list_path
        self.preprocessed_intervals_path = preprocessed_intervals_path
        self.annotated_intervals = annotated_intervals
        self.bed_volume = "/bed_dir

    def get_interval_list_path(self):
        if self.interval_list_path:
            return(self.interval_list_path)
        else:
            raise(ValueError(
                f"Specify a interval_list_path before trying to obtain it, you can obtain by:\n\
                \tSpecifying in Bed class as input parameter\n\
                \tRunning Picard.get_intervals_from_bed()."))