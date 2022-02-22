function diversity_metrics = Shannon_diversity(shannon_matrix)
diversity_metrics = nansum(-log(shannon_matrix).*shannon_matrix);
end