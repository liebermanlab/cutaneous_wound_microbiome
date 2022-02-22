function dist_mat = create_distance_matrix(x,y)
dist_mat = nan(size(x,1), size(y,1));
for ii =1:size(x,1)
    x_temp = x(ii,:);
    for ij= 1:size(y,1)
        y_temp = y(ij,:);
        %dist_mat(ii,ij) = sqrt(sum((x_temp-y_temp).^2));
        dist_mat(ii,ij) = (sum((x_temp-y_temp).^2));
    end
end
end