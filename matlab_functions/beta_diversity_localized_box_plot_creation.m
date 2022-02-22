function holding_mat = beta_diversity_localized_box_plot_creation(input_mat, matched_surgical_control_locations,n_term)
holding_mat = nan(length(matched_surgical_control_locations)*2, 3);

% Intragroup comparison - control
control_distance = create_distance_matrix(input_mat(matched_surgical_control_locations(:,1),1:2),input_mat(matched_surgical_control_locations(:,1),1:2));
for ii=1:size(control_distance,2)
    sorted_distances = sort(control_distance(:,ii),'ascend');
    holding_mat(ii,1) = mean(sorted_distances(1:n_term));
end
% Intragroup comparison - surgical
surgical_distance = create_distance_matrix(input_mat(matched_surgical_control_locations(:,2),1:2),input_mat(matched_surgical_control_locations(:,2),1:2));
for ii=1:size(surgical_distance,2)
    sorted_distances = sort(surgical_distance(:,ii),'ascend');
    holding_mat(ii,2) = mean(sorted_distances(1:n_term));
end

% Intergroup comparison
surg_cont_distance = create_distance_matrix(input_mat(matched_surgical_control_locations(:,1),1:2),input_mat(matched_surgical_control_locations(:,2),1:2));
for ii=1:size(surg_cont_distance,2)
    sorted_distances = sort(surg_cont_distance(:,ii),'ascend');
    holding_mat(ii,3) = mean(sorted_distances(1:n_term));
end
for ii=1:size(surg_cont_distance,1)
    sorted_distances = sort(control_distance(ii,:),'ascend');
    holding_mat(ii+size(surg_cont_distance,2),3) = mean(sorted_distances(1:5));
end

end
