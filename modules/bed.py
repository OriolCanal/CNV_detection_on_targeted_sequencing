import os

class Bed():
    def __init__(
            self,
            bed_path,
            interval_list_path="",
            annotated_intervals_path="",
            preprocessed_intervals_path=""
        ):
        self.path = bed_path
        self.filename = os.path.basename(self.path)
        self.dir = os.path.dirname(self.path)

        self.interval_list_path = interval_list_path
        self.interval_list_filename = os.path.basename(self.interval_list_path)

        self.preprocessed_intervals_path = preprocessed_intervals_path
        self.preprocessed_intervals_filename = os.path.basename(self.preprocessed_intervals_path)

        self.annotated_intervals_path = annotated_intervals_path
        self.annotated_intevals_filename = os.path.basename(annotated_intervals_path)

        self.volume = "/bed_dir"

    def get_interval_list_path(self):
        if self.interval_list_path:
            return(self.interval_list_path)
        else:
            raise(ValueError(
                f"Specify a interval_list_path before trying to obtain it, you can obtain by:\n\
                \tSpecifying in Bed class as input parameter\n\
                \tRunning Picard.get_intervals_from_bed()."))
        
    def set_interval_list_path(self, interval_list_path):
        if os.path.isfile(interval_list_path) and os.path.getsize(interval_list_path) > 0:
            self.interval_list_path = interval_list_path
        else:
            raise ValueError(
                f"interval list path does not exist: {interval_list_path}"
            )
        
        self.interval_list_filename = os.path.basename(interval_list_path)
    
    def set_preprocessed_intervals_list_path(self, preprocessed_intervals_path):
        if os.path.isfile(preprocessed_intervals_path) and os.path.getsize(preprocessed_intervals_path) > 0:
            self.preprocessed_intervals_path = preprocessed_intervals_path
        else:
            raise ValueError(
                f"Preprocessed intervals path does not exists: {preprocessed_intervals_path}"
            )
        self.preprocessed_intervals_filename = os.path.basename(preprocessed_intervals_path)
    
    def set_annotated_intervals_path(self, annotated_intervals_path):
        if os.path.isfile(annotated_intervals_path) and os.path.getsize(annotated_intervals_path) > 0:
            self.annotated_intervals_path = annotated_intervals_path
        else:
            raise ValueError(
                f"Annotated intervals path does not exists: {annotated_intervals_path}"
            )
        self.annotated_intevals_filename = os.path.basename(annotated_intervals_path)

    def set_filtered_intervals_path(self, filtered_intervals_path):
        if os.path.isfile(filtered_intervals_path) and os.path.getsize(filtered_intervals_path):
            self.filtered_intervals_path = filtered_intervals_path
        else:
            raise ValueError(
                f"Filtered intervals path does not exists: {filtered_intervals_path}"
            )
        self.filtered_intervals_filename = os.path.basename(filtered_intervals_path)