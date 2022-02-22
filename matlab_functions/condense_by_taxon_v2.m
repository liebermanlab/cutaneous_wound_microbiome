function condensed_matrix = condense_by_taxon_v2(input_matrix,taxon_labels, taxons_to_show,toggle_add_other)
    % Function that combines ASVs that match a label in taxon_labels
    % Input: input_matrix, an ASV matrix
    %        taxon_labels, taxon labels that match the ASV matrix
    %        taxons_to_show, list of taxons to condense 
    %        toggle_add_other, a logical of whether or not to include a 
    %        single line containing all other taxon rows not contained in taxons_to_show
    condensed_matrix = zeros(length(taxons_to_show)+toggle_add_other,size(input_matrix,2));
    taxon_labels = string(taxon_labels);
    in_taxon_labels = zeros(length(taxon_labels),1);
    for ij =1:length(taxons_to_show)
        locations_of_taxons = taxon_labels==taxons_to_show{ij};
        condensed_matrix(ij,:) = sum(input_matrix(locations_of_taxons,:),1);
        in_taxon_labels = in_taxon_labels | locations_of_taxons;    
    end
    % Add a line containing all other taxons not included in the
    % taxon_labels list
    if toggle_add_other
        condensed_matrix(length(taxons_to_show)+1,:) = sum(input_matrix(~in_taxon_labels,:),1);
    end
end