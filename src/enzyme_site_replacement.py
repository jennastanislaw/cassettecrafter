# take a sequence and an enzyme and find cut sites
import enzyme_site_replacement_utils

def replace_enzyme_site(enzyme_fp, enzyme_name, DNA):

    enzyme_class_objs = enzyme_site_replacement_utils.load_enzymes_from_csv(enzyme_fp)
    enzyme_dict = enzyme_site_replacement_utils.create_enzyme_dict(enzyme_class_objs)

    enzyme = enzyme_dict.get(enzyme_name)

    # Add this check to raise KeyError for invalid enzyme names
    if enzyme is None:
        raise KeyError(f"Enzyme '{enzyme_name}' not found.")

    DNA_rec_sites_rem = list(DNA)

    fwd_matches, rev_matches = enzyme_site_replacement_utils.find_matching_sites(enzyme, DNA)
    fwd_codons = enzyme_site_replacement_utils.get_affected_codons_by_recognition_sites(DNA, fwd_matches, enzyme.fwd_recognition_site)
    rev_codons = enzyme_site_replacement_utils.get_affected_codons_by_recognition_sites(DNA, rev_matches, enzyme.rev_recognition_site)

    all_codons = {**fwd_codons, **rev_codons}

    # print(all_codons)

    for key, value in all_codons.items():
        # print(key, value)
        index_vals = [codon_tuple[1] for codon_tuple in value]
        flattened_index_values = [index for sublist in index_vals for index in sublist]
        max_idx = max(flattened_index_values)
        min_idx = min(flattened_index_values)

        replaced = False

        for codon_tuple in value:
            codon = codon_tuple[0]
            codon_idx = codon_tuple[1][0]
            alt_codons = enzyme_site_replacement_utils.generate_synonymous_codons_dna(codon)
            for alt_codon in alt_codons:
                DNA_rec_sites_rem[codon_idx:codon_idx+3] = alt_codon
                # check if this removed the site
                DNA_check = ''.join(DNA_rec_sites_rem)
                fwd_matches, rev_matches = enzyme_site_replacement_utils.find_matching_sites(enzyme, DNA_rec_sites_rem)
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
            raise ValueError(f"No suitable replacement found for recognition site at index {key}.")

    new_DNA=''.join(DNA_rec_sites_rem)

    return new_DNA
