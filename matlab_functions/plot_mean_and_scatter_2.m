function plot_mean_and_scatter_2(datasets_to_plot,x_values,color_vals,figure_toggle, std_toggle,sem_toggle,ci_toggle, batch_num)

if figure_toggle
    figure;
end
for ii =1:length(datasets_to_plot)
    metric_to_plot = datasets_to_plot{ii};
    mean_to_plot = mean(metric_to_plot);
    scatter(x_values(ii).*ones(1,length(metric_to_plot(batch_num==1))),metric_to_plot(batch_num==1),[],color_vals(ii,:),'o','filled', 'jitter', 'on', 'jitterAmount', 0.3);
    hold on
    scatter(x_values(ii).*ones(1,length(metric_to_plot(batch_num==2))),metric_to_plot(batch_num==2),[],color_vals(ii,:),'^','filled', 'jitter', 'on', 'jitterAmount', 0.3);
    hold on
    plot([x_values(ii)-.3,x_values(ii)+.3],[mean_to_plot,mean_to_plot],'k','Linewidth',1.5);
    if std_toggle
        standard_dev = std(metric_to_plot);
        errlow = standard_dev;
        errhigh = standard_dev;
    end
    if sem_toggle
        sem_val = get_sem(metric_to_plot);
        errlow = sem_val;
        errhigh = sem_val;        
    end
    if ci_toggle
        confidence_int_vals = confidence_int(metric_to_plot);
        errlow = mean_to_plot - confidence_int_vals(1);
        errhigh = confidence_int_vals(2) - mean_to_plot;          
    end
    if std_toggle || sem_toggle || ci_toggle
        er = errorbar(x_values(ii),mean_to_plot,errlow,errhigh,'linewidth',1.5);    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  
    end
    hold on
end
set(gca,'xtick',[]); set(gca,'fontsize',16); set(gca,'TickDir','out');set(gca,'linewidth',1)

end