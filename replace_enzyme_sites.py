import csv

class Enzyme:
    def __init__(self, name, fwd_recognition_site, rev_recognition_site, spacer_length, OH_length):
        self.name = name
        self.fwd_recognition_site = fwd_recognition_site
        self.rev_recognition_site = rev_recognition_site
        self.spacer_length = int(spacer_length)
        self.OH_length = int(OH_length)

    def __repr__(self):
        return (f"Enzyme(name={self.name}, fwd_recognition_site={self.fwd_recognition_site}, "
                f"rev_recognition_site={self.rev_recognition_site}, spacer_length={self.spacer_length}, "
                f"OH_length={self.OH_length})")

def load_enzymes_from_csv(csv_file_path):
    """Reads a CSV file and creates Enzyme objects for each row."""
    enzymes = []
    with open(csv_file_path, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            enzyme = Enzyme(
                name=row['Enzyme'],
                fwd_recognition_site=row['Fwd_recognition_site'],
                rev_recognition_site=row['Rev_recognition_site'],
                spacer_length=row['spacer_length'],
                OH_length=row['OH_length']
            )
            enzymes.append(enzyme)
    return enzymes


#def find_IIs_sites(DNA, enzyme):
