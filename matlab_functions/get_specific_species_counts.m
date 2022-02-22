function [new_normalized_count_matrix_condensed] = get_specific_species_counts(taxons_to_get,taxon_list, counts_matrix_normalized)

species_not_to_count = zeros(length(taxon_list),1);
new_normalized_count_matrix_condensed = zeros(length(taxons_to_get)+1, size(counts_matrix_normalized,2));

for ii=1:length(taxons_to_get)
    contains_spec_name = contains(taxon_list, taxons_to_get(ii)) & ~species_not_to_count;
    species_not_to_count = species_not_to_count | contains_spec_name;
    new_normalized_count_matrix_condensed(ii,:) = sum(counts_matrix_normalized(contains_spec_name,:),1);
end

new_normalized_count_matrix_condensed(end,:) = sum(counts_matrix_normalized(~species_not_to_count,:),1);

end

