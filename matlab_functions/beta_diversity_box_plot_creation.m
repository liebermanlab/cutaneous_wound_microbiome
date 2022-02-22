function holding_mat = beta_diversity_box_plot_creation(input_mat, matched_surgical_control_locations)

control_distance = create_distance_matrix(input_mat(matched_surgical_control_locations(:,1),1:2),input_mat(matched_surgical_control_locations(:,1),1:2));
control_distance = control_distance(triu(ones(size(control_distance)),1)==1);
surgical_distance = create_distance_matrix(input_mat(matched_surgical_control_locations(:,2),1:2),input_mat(matched_surgical_control_locations(:,2),1:2));
surgical_distance = surgical_distance(triu(ones(size(surgical_distance)),1)==1);
surg_cont_distance = create_distance_matrix(input_mat(matched_surgical_control_locations(:,1),1:2),input_mat(matched_surgical_control_locations(:,2),1:2));
surg_cont_distance = surg_cont_distance(:);

holding_mat = nan(length(surg_cont_distance), 3);
holding_mat(1:length(control_distance),1) = control_distance;
holding_mat(1:length(surgical_distance),2) = surgical_distance;
holding_mat(1:length(surg_cont_distance),3) = surg_cont_distance;
end
