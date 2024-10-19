# take a sequence and an enzyme and find cut sites

import replace_enzyme_sites

enzymes_fp = './test_data/enzyme_sites.csv'

enzymes = replace_enzyme_sites.load_enzymes_from_csv(enzymes_fp)

print(enzymes)