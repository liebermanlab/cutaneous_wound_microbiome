function [new_taxons, new_counts, new_counts_normalized,indexes_of_bad_taxons]= remove_taxons_index(bad_taxa,taxon_labels, counts_matrix)

    indexes_of_bad_taxons = zeros(size(taxon_labels));
    for ii=1:length(bad_taxa)
        indexes_of_bad_taxons(taxon_labels==bad_taxa(ii)) = 1;
    end
    
    new_taxons = taxon_labels(~indexes_of_bad_taxons);
    new_counts = counts_matrix(~indexes_of_bad_taxons,:);
    new_counts_normalized = new_counts./sum(new_counts);

end