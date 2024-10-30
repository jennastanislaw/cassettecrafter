# take a sequence and an enzyme and find cut sites

import replace_enzyme_sites


enzymes_fp = './test_data/enzyme_sites.csv'

enzyme_class_objs = replace_enzyme_sites.load_enzymes_from_csv(enzymes_fp)

enzyme_dict = replace_enzyme_sites.create_enzyme_dict(enzyme_class_objs)

enzyme = enzyme_dict.get('BbsI')

DNA = "ATCGAAGACGTCGTCTTCAGCTGAAGACGTCTTCCAGT"

print(len(DNA))

DNA_rec_sites_rem = list(DNA)

fwd_matches, rev_matches = replace_enzyme_sites.find_matching_sites(enzyme, DNA)

print(fwd_matches, rev_matches)

fwd_codons = replace_enzyme_sites.get_affected_codons_by_recognition_sites(DNA, fwd_matches, enzyme.fwd_recognition_site)

rev_codons = replace_enzyme_sites.get_affected_codons_by_recognition_sites(DNA, rev_matches, enzyme.rev_recognition_site)

all_codons = {**fwd_codons, **rev_codons}

print(all_codons)

for key, value in all_codons.items():
    print(key, value)
    index_vals = [codon_tuple[1] for codon_tuple in value]
    flattened_index_values = [index for sublist in index_vals for index in sublist]
    max_idx = max(flattened_index_values)
    min_idx = min(flattened_index_values)

    replaced = False

    for codon_tuple in value:
        codon = codon_tuple[0]
        codon_idx = codon_tuple[1][0]
        alt_codons = replace_enzyme_sites.generate_synonymous_codons_dna(codon)
        for alt_codon in alt_codons:
            DNA_rec_sites_rem[codon_idx:codon_idx+3] = alt_codon
            # check if this removed the site
            DNA_check = ''.join(DNA_rec_sites_rem)
            fwd_matches, rev_matches = replace_enzyme_sites.find_matching_sites(enzyme, DNA_rec_sites_rem)
            match_idx=fwd_matches+rev_matches
            if key not in match_idx:
                print(f"Found replacement for codon {codon} for recognition site at index {key}, replaced with {alt_codon}")
                replaced = True
                break

        if not replaced:
            continue

        if replaced:
            break

    if not replaced:
        print(f"Alternative codon for recognition site at index {key} could not be found")

print(''.join(DNA_rec_sites_rem))

new_DNA=''.join(DNA_rec_sites_rem)
