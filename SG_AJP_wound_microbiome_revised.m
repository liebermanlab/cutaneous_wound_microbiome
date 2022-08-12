%% Preprocess files
clear all; close all;

% Load in all raw data files
metadata_file  = "patient_sequencing_metadata.xlsx";
batch_2_file = "batch_3_QIIME2_output.xlsx";
batch_1_file = "batch_1_2_QIIME2_output.xlsx";
date_for_file_names = date;

% Import data into useable data structures 
rawdata_1 = readtable(batch_1_file); rawdata_1 = rawdata_1(:,[1:4,22:end]); 
rawdata_2 = readtable(batch_2_file); rawdata_2 = rawdata_2(:,[1:4,26:end]);

% Import metadata 
metadata = readtable(metadata_file); 
metadata_string_holder = string(metadata{:,[2,3,4]});
metadata_string_holder(metadata_string_holder(:,2)==metadata_string_holder(:,3),3) = "";
time_of_sampling = str2double(string(metadata{:,5}));batches = str2double(string(metadata{:,6}));
 
% Create a new naming scheme for every sample 
% Format (connected by dashes): 
%    Patient Number and Sample Distinguisher (A/B/C...etc) 
%    Sample type ( X = surgical, C = contralateral control, B = other skin control, oral, air, and nares samples lack a label)  
%    Site description (ex. right nose, temple, scalp, air, oral, etc)

new_names = join(metadata_string_holder,"-"); new_names= strrep(new_names,"-c-","-B-");
new_names= strrep(new_names,"control","C"); new_names= strrep(new_names,"surgery","X"); 
remove_last_char = cellfun(@(x) x(end),new_names)=='-';
new_names(remove_last_char) = string(cellfun(@(x) x(1:end-1),new_names(remove_last_char),'UniformOutput',false));

column_names_1 = rawdata_1.Properties.VariableNames; 
column_names_2 = rawdata_2.Properties.VariableNames;
samplenames_1 = string(column_names_1(5:end)); 
samplenames_2 = strrep(string(column_names_2(5:end)),'_','-');
samplenames_1 = string(cellfun(@(x) x(1:11),samplenames_1,'UniformOutput',false));

new_samplenames_1 = cellfun(@(x) new_names(contains(table2array(metadata(:,1)),x)),samplenames_1);
new_samplenames_2 = cellfun(@(x) new_names(contains(table2array(metadata(:,1)),x)),samplenames_2);
sample_time_2 = cellfun(@(x) time_of_sampling(contains(table2array(metadata(:,1)),x)),samplenames_2);
sample_time_1 = cellfun(@(x) time_of_sampling(contains(table2array(metadata(:,1)),x)),samplenames_1);

collection_time = [sample_time_1,sample_time_2];
new_samplenames = [new_samplenames_1, new_samplenames_2];
batch_num= zeros(size(new_samplenames));
batch_num(contains(new_samplenames, string([14:43]))) = 1;
batch_num(contains(new_samplenames, string([44:69]))) = 2;
batch_num(contains(new_samplenames, string([70:88]))) = 3;

% Create a new table structure that combines both raw ASV count datasets 
raw_counts_1 = table2array(rawdata_1(:,5:end)); raw_counts_2 = table2array(rawdata_2(:,5:end));
all_counts_1 = [raw_counts_1, zeros(size(raw_counts_1,1), size(raw_counts_2,2))]; all_counts_2 = [zeros(size(raw_counts_2,1), size(raw_counts_1,2)), raw_counts_2];
combined_all_counts = [all_counts_1;all_counts_2];

raw_combined_table = table([rawdata_1.x_OTUID;rawdata_2.x_OTUID],[rawdata_1.FASTA;rawdata_2.Fasta],[rawdata_1.Confidence;rawdata_2.Confidence],[rawdata_1.Taxon;rawdata_2.Taxon],'VariableNames', rawdata_1.Properties.VariableNames(1:4)); % Combine FATSA and OTU labels
raw_combined_table = [raw_combined_table,array2table(combined_all_counts,'VariableNames', new_samplenames)]; % Reform this into a table 

% Condense table by OTU ID (technically redundant)
[unique_OTU_ID, unique_OTU_ID_indexes] = unique(raw_combined_table.x_OTUID,'stable');
condensed_ASV_matrix = condense_by_taxon_v2(combined_all_counts,raw_combined_table.x_OTUID,unique_OTU_ID ,0);
combined_table = [raw_combined_table(unique_OTU_ID_indexes,1:4),array2table(condensed_ASV_matrix,'VariableNames', new_samplenames)];

%% Remove ASVs that are batch contaminants
% Sequencing note:
% We had 2 patients that did not return for sampling (69 and 73), as well
% as some patients that returned too late (after 8 days). These samples are
% removed from this analysis. 

control_sample = new_samplenames(contains(new_samplenames, "-C-") & ~contains(new_samplenames,["69","73"]) & collection_time<=8); 
% Patients "69","73" didn't complete sampling protocols and are removed
% from all analysis; similarly, all patients sampled after a week are
% removed as they didn't complete the nessesary study protocol. 

[~,first_sample_index] = unique(string(cellfun(@(x) x(1:2), control_sample, 'UniformOutput', false)));
first_sample = zeros(size(new_samplenames)); first_sample(contains(new_samplenames,control_sample(first_sample_index))) = 1;
% We grab one control sample from healthy skin from each patient to
% identify batch contaminants 

% Determine contaminant samples by batch
condensed_ASV_matrix_norm = condensed_ASV_matrix./sum(condensed_ASV_matrix);
p_vals = nan(size(condensed_ASV_matrix_norm,1), 3);
for ii=1:size(condensed_ASV_matrix_norm,1)
    p_vals(ii,1) = ranksum(condensed_ASV_matrix_norm(ii,batch_num==1 & first_sample),condensed_ASV_matrix_norm(ii,batch_num~=1 & first_sample),'tail','right');
    p_vals(ii,2) = ranksum(condensed_ASV_matrix_norm(ii,batch_num==2 & first_sample),condensed_ASV_matrix_norm(ii,batch_num~=2 & first_sample),'tail','right');
    p_vals(ii,3) = ranksum(condensed_ASV_matrix_norm(ii,batch_num==3 & first_sample),condensed_ASV_matrix_norm(ii,batch_num~=3 & first_sample),'tail','right');
end
min_p_vals = min(p_vals,[],2);

p_value_corrected_cutoff = 4e-5;

sig_batch_taxons = min_p_vals<p_value_corrected_cutoff;

figure; histogram(log10(min_p_vals)); hold on; plot([log10(p_value_corrected_cutoff),log10(p_value_corrected_cutoff)],[.5,10^4.3],'r');
set(gca,'fontsize',16); set(gca,'TickDir','out'); set(gca,'linewidth',1); box off;
xlabel("Log10 p-value"); ylabel("Count"); set(gca, 'yscale','log'); ylim([.5,10^4.1])

av_abundance_across_batch = [mean(condensed_ASV_matrix_norm(sig_batch_taxons>0,batch_num==1),2),mean(condensed_ASV_matrix_norm(sig_batch_taxons>0,batch_num==2),2),mean(condensed_ASV_matrix_norm(sig_batch_taxons>0,batch_num==3),2)];
view_contams = [combined_table.Taxon(sig_batch_taxons>0),num2cell(av_abundance_across_batch),num2cell(min_p_vals(sig_batch_taxons>0))];

indexes_of_removed_taxons = find(sig_batch_taxons>0);
mean_removed_taxons_not_bacteria = mean(sum(condensed_ASV_matrix_norm(indexes_of_removed_taxons(combined_table.Taxon(indexes_of_removed_taxons)~="Bacteria"),:)));
mean_salinimic_taxons = mean(sum(condensed_ASV_matrix_norm(indexes_of_removed_taxons(contains(combined_table.Taxon(indexes_of_removed_taxons),"Salinimicrobium")),batch_num==3)));
mean_not_salinimic_taxons = mean(sum(condensed_ASV_matrix_norm(indexes_of_removed_taxons(~contains(combined_table.Taxon(indexes_of_removed_taxons),"Salinimicrobium")),batch_num==3)));

mean_removed_taxons_all_bacteria = mean(sum(condensed_ASV_matrix_norm(indexes_of_removed_taxons,:)));


ASVs_to_remove = sig_batch_taxons;

%% Remove human mapping phyla and spike ins
phyla_taxon = cellfun(@(x) count(x, ";"), combined_table.Taxon)==1;
all_phyla_taxons = string(unique(combined_table.Taxon(phyla_taxon)));

bad_taxon_labels = ["Bacteria","Unassigned","Bacteria;Proteobacteria;Gammaproteobacteria;Thiomicrospirales;Thioglobaceae;SUP05 cluster"];
[~, ~, ~,bacteria_index]= remove_taxons_index(bad_taxon_labels(1),combined_table.Taxon, condensed_ASV_matrix); % Get indexes of bad taxons
[~, ~, ~,unassigned_index]= remove_taxons_index(bad_taxon_labels(2),combined_table.Taxon, condensed_ASV_matrix); % Get indexes of bad taxons
[~, ~, ~,SUP05_index]= remove_taxons_index(bad_taxon_labels(3),combined_table.Taxon, condensed_ASV_matrix); % Get indexes of bad taxons

[~, ~, ~,proteobacteria_index]= remove_taxons_index(all_phyla_taxons(3),combined_table.Taxon, condensed_ASV_matrix); % Get indexes of bad taxons
[~, ~, ~,firmicutes_index]= remove_taxons_index(all_phyla_taxons(2),combined_table.Taxon, condensed_ASV_matrix); % Get indexes of bad taxons
[~, ~, ~,actinobacteria_index]= remove_taxons_index(all_phyla_taxons(1),combined_table.Taxon, condensed_ASV_matrix); % Get indexes of bad taxons

blastable_datastruct = [">" +  combined_table.x_OTUID, combined_table.FASTA];
rng(13579)
random_seq_num = 25;
bacteria_blast = blastable_datastruct(logical(bacteria_index),:)'; bacteria_blast = bacteria_blast(:,randperm(sum(bacteria_index),random_seq_num)); bacteria_blast = bacteria_blast(:);
unassigned_blast = blastable_datastruct(logical(unassigned_index),:)'; unassigned_blast = unassigned_blast(:,randperm(sum(unassigned_index),random_seq_num)); unassigned_blast = unassigned_blast(:);
SUP05_blast = blastable_datastruct(logical(SUP05_index),:)'; SUP05_blast = SUP05_blast(:,:); SUP05_blast = SUP05_blast(:);

proteobacteria_blast = blastable_datastruct(logical(proteobacteria_index),:)'; proteobacteria_blast = proteobacteria_blast(:,randperm(sum(proteobacteria_index),random_seq_num)); proteobacteria_blast = proteobacteria_blast(:);
firmicutes_blast = blastable_datastruct(logical(firmicutes_index),:)'; firmicutes_blast = firmicutes_blast(:,:); firmicutes_blast = firmicutes_blast(:);
actinobacteria_blast = blastable_datastruct(logical(actinobacteria_index),:)'; actinobacteria_blast = actinobacteria_blast(:,:); actinobacteria_blast = actinobacteria_blast(:);

figure; 
scatter(batch_num, sum(condensed_ASV_matrix(SUP05_index==1,:))./sum(condensed_ASV_matrix),'k','filled',jitter=.2)
set(gca,'fontsize',16); set(gca,'TickDir','out'); set(gca,'linewidth',1); box off;xticks(1:3);
xlabel("Batch number"); xlim([.5,3.5]); ylabel("Fraction Thioglobaceae;SUP05 in Sample")

figure('Renderer', 'painters', 'Position', [10 10 700 400])
scatter(batch_num, sum(condensed_ASV_matrix(bacteria_index==1,:))./sum(condensed_ASV_matrix),...
    [],'k','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,jitter=.2)
hold on
scatter(batch_num, sum(condensed_ASV_matrix(proteobacteria_index==1,:))./sum(condensed_ASV_matrix),...
    [],'b','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,jitter=.2)
hold on
scatter(batch_num, sum(condensed_ASV_matrix(unassigned_index==1,:))./sum(condensed_ASV_matrix),...
    [],'r','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,jitter=.2)
set(gca,'fontsize',16); set(gca,'TickDir','out'); set(gca,'linewidth',1); box off;xticks(1:3);
xlabel("Batch number"); xlim([.5,3.5]); ylabel("Fraction ASV in Sample")
legend(["Bacteria","Proteobacteria","Unassigned"],'Location','eastoutside')

%% Remove contaminant and human taxons and create final arrays for analysis
ASVs_to_remove = ASVs_to_remove | bacteria_index | unassigned_index | proteobacteria_index | SUP05_index;

toggle_remove_ASVs = 1;
if toggle_remove_ASVs
    combined_table(ASVs_to_remove>0,:) = [];
end

% Condense data by taxon 
unique_taxons =  unique(combined_table.Taxon,'stable');
condensed_counts = condense_by_taxon_v2(table2array(combined_table(:,5:end)),combined_table.Taxon, unique_taxons,0);
condensed_counts_normalized = condensed_counts./sum(condensed_counts);

% Create an excite file with the raw ASV counts for sharing
writeable_count_table = array2table(table2array(combined_table(:,5:end)), 'VariableNames', new_samplenames,'RowNames',string(combined_table.x_OTUID));
writetable(writeable_count_table,"saved_excel_files/" + date +  "_all_unique_ASVs.xlsx", 'WriteVariableNames', true,'WriteRowNames', true);

fasta_list = [">" + string(combined_table.x_OTUID),string(combined_table.FASTA)]';
fasta_list = fasta_list(:);

% Create an excite file with the raw taxon counts for sharing
writeable_count_table = array2table(condensed_counts, 'VariableNames', new_samplenames,'RowNames',unique_taxons);
writetable(writeable_count_table,"saved_excel_files/" + date +  "_all_unique_taxons.xlsx", 'WriteVariableNames', true,'WriteRowNames', true);

%% Get genus and species level classifications 
% Split taxons by level
split_by_classification = cellfun(@(x) split(x, ';'), unique_taxons, 'UniformOutput', false);
taxonomic_classes_split = repmat({strings(size(split_by_classification))}, 7,1);
for ii=1:7
    is_classified_in_class = cellfun(@length, split_by_classification) >=ii;
    taxonomic_classes_split{ii}(is_classified_in_class) = cellfun(@(x) x(ii), split_by_classification(is_classified_in_class));
end

unique_genuses = unique(taxonomic_classes_split{6}, "stable"); unique_genuses(unique_genuses == "") = [];
[genus_counts_unnormalized] = get_specific_species_counts(unique_genuses',taxonomic_classes_split{6}, condensed_counts);
unique_genuses = [unique_genuses; 'Other'];

writeable_genus = array2table(genus_counts_unnormalized, 'VariableNames', new_samplenames,'RowNames',[unique_genuses]);
writetable(writeable_genus,"saved_excel_files/" + date +  "_unique_genus_counts.xlsx", 'WriteVariableNames', true,'WriteRowNames', true);

unique_species = unique(taxonomic_classes_split{7}, "stable"); unique_species(unique_species == "") = [];
[species_counts_unnormalized] = get_specific_species_counts(unique_species',unique_taxons, condensed_counts);
unique_species = [unique_species; 'Other'];

writeable_species = array2table(species_counts_unnormalized, 'VariableNames', new_samplenames,'RowNames',[unique_species]);
writetable(writeable_species,"saved_excel_files/" + date +  "_species_counts.xlsx", 'WriteVariableNames', true,'WriteRowNames', true);

%% Further quality control metrics 
clinical_infections = ["30","41"]; prophylactic_ab = ["36","75","86","78","80","82","84"]; patients_to_remove = [clinical_infections,prophylactic_ab,"69","73"];
% clinical_infections removes patients with clinical infections
% prophylactic_ab patients removes those who took prophylactic antibiotics
% patients_to_remove removes patients "69" and "73", both of whom did not
% complete sampling protocols (missing day ~7 samples and a control site,
% respectively) 


number_of_pairs_remaining = (~contains(new_samplenames,patients_to_remove) & collection_time<=8 & contains(new_samplenames,"-X-") &  collection_time~=0);

read_cutoff = 500;
good_samples = sum(condensed_counts)>read_cutoff & ~contains(new_samplenames,'-air') &  ~contains(new_samplenames,patients_to_remove) & collection_time<=8;

control_samps = contains(new_samplenames,"-C-") & good_samples & collection_time~=0;
surgical_samps = contains(new_samplenames,"-X-") & good_samples & collection_time~=0;
day_0_control_samps = contains(new_samplenames,"-C-") & good_samples & collection_time==0;
day_0_surg_samps = contains(new_samplenames,"-X-") & good_samples & collection_time==0;

% reason_for_removal can be manually parsed to see why samples were removed
% reason_for_removal(1,:) = sample failed to pass read cutoff
% reason_for_removal(2,:) = sample comes from removed patient
% reason_for_removal(3,:) = collection time outside cutoff bound (above 8 days)
%   note: some patients (which were not included in the study description)
%   failed to return in time sampling <8 days later. These are not included
%   in any study plots. 
% reason_for_removal(4,:) = pair sample failed to pass cutoff
% reason_for_removal(5,:) = sample doubled from the same patient

reason_for_removal = zeros(5,length(new_samplenames));
reason_for_removal(1,:) = sum(condensed_counts)<=read_cutoff;
reason_for_removal(2,:) = contains(new_samplenames,patients_to_remove);
reason_for_removal(3,:) = collection_time>8 ;

%% Match up samples by body site 
surgical_samps_split = split(new_samplenames(surgical_samps)','-');
patients_to_parse = string(cellfun(@(x) x(1:2), surgical_samps_split(:,1), 'UniformOutput', false));
locations_to_parse = erase(surgical_samps_split(:,3),["L ","R ","inferior ","superior "," wound"," normal"]);

matched_surgical_control_locations = [];
matched_surgical_control_patient_list = [];
for ij = 1:length(patients_to_parse)
    possible_control_samps = contains(new_samplenames, patients_to_parse(ij)) & contains(new_samplenames, locations_to_parse(ij)) & control_samps;
    possible_surgical_samps = contains(new_samplenames, patients_to_parse(ij)) & contains(new_samplenames, locations_to_parse(ij)) & surgical_samps;
    if sum(possible_surgical_samps)>1
        disp(ij)
    end
    if sum(possible_control_samps)==1 && sum(possible_surgical_samps)==1
        matched_surgical_control_locations = [matched_surgical_control_locations; find(possible_control_samps), find(possible_surgical_samps)];
        matched_surgical_control_patient_list = [matched_surgical_control_patient_list; repelem(patients_to_parse(ij),1,sum(possible_control_samps))];
    elseif sum(possible_control_samps)>1 || sum(possible_surgical_samps)>1
        disp(new_samplenames(possible_control_samps))
        disp(new_samplenames(possible_surgical_samps))
    else
        reason_for_removal(4,:) = reason_for_removal(4,:)| possible_control_samps | possible_surgical_samps;
    end
end

% Manually add in a right left nose sample pair and glabella-nose pair that would be confusing to parse otherwise
matched_surgical_control_locations = [matched_surgical_control_locations; find(contains(new_samplenames, '19B')),find(contains(new_samplenames, '19A'))];
matched_surgical_control_patient_list = [matched_surgical_control_patient_list;[19]];

% Removes the second instance of a pair match from the same patient 
[~,unique_matches] = unique(matched_surgical_control_patient_list,'stable');
unique_matches_2 = ones(length(matched_surgical_control_patient_list),1); unique_matches_2(unique_matches) = 0;
reason_for_removal(5,matched_surgical_control_locations(unique_matches_2>0,:)) = 1;
matched_surgical_control_locations = matched_surgical_control_locations(unique_matches,:);

%% Find genuses and species present in matched samples

% Genuses
genus_counts_normalized = genus_counts_unnormalized./sum(genus_counts_unnormalized);
[abundant_wound_genuses,abundant_wound_genus_index] = sort(sum(genus_counts_normalized(:,matched_surgical_control_locations(:,1:2)),2),'descend');
sorted_genuses = unique_genuses(abundant_wound_genus_index);

% Species
species_counts_normalized = species_counts_unnormalized./sum(species_counts_unnormalized);
[abundant_wound_species,abundant_wound_species_index] = sort(sum(species_counts_unnormalized(:,matched_surgical_control_locations(:)),2),'descend');
sorted_species = unique_species(abundant_wound_species_index);

%% Grab the 8 most abundant genuses

genuses_to_show_figure = ["Cutibacterium";"Corynebacterium";"Staphylococcus";"Streptococcus";"Pseudomonas";"Neisseria";"Anaerococcus";"Finegoldia"];
figure_colors= ["a6cee3","1f78b4","b2df8a","52974d","cca175","ffba60","fdabaa","ffe2e3","757575"]; %figure_colors= ["a6cee3", "1f78b4", "b2df8a", "33a02c", "fb9a99", "e31a1c", "fdbf6f", "ff7f00","757575"];

all_colors_in_order = zeros(length(figure_colors), 3);
for ii=1:length(figure_colors)
    all_colors_in_order(ii,:) = reshape(sscanf(figure_colors(ii).','%2x'),3,[]).'/255;
end

sorted_genus_matrix_figure = get_specific_species_counts(genuses_to_show_figure,unique_genuses, genus_counts_normalized);
count_matrix_control = sorted_genus_matrix_figure(:,matched_surgical_control_locations(:,1));
count_matrix_surgical = sorted_genus_matrix_figure(:,matched_surgical_control_locations(:,2));


%% Grab genuses and additionally S. aureus
% This code is a little redundant, but it is useful to seperate out in
% order to identify how seperating out S. aureus impacts the Staph genus as a whole  

genuses_to_show_figure_2 = ["Cutibacterium";"Corynebacterium";"Staphylococcus aureus";"Staphylococcus";"Streptococcus";"Pseudomonas";"Neisseria";"Anaerococcus";"Finegoldia"];
figure_colors_2= ["a6cee3","1f78b4","52974d","b2df8a","d9d9d9","cca175","ffba60","fdabaa","ffe2e3","757575"]; %figure_colors= ["a6cee3", "1f78b4", "b2df8a", "33a02c", "fb9a99", "e31a1c", "fdbf6f", "ff7f00","757575"];

all_colors_in_order_2 = zeros(length(figure_colors_2), 3);
for ii=1:length(figure_colors_2)
    all_colors_in_order_2(ii,:) = reshape(sscanf(figure_colors_2(ii).','%2x'),3,[]).'/255;
end

sorted_genus_matrix_figure_2 = get_specific_species_counts(genuses_to_show_figure_2,unique_taxons, condensed_counts_normalized);
count_matrix_control_2 = sorted_genus_matrix_figure_2(:,matched_surgical_control_locations(:,1));
count_matrix_surgical_2 = sorted_genus_matrix_figure_2(:,matched_surgical_control_locations(:,2));
[~,display_index_control_2] = sort(count_matrix_control_2(1,:),'descend');
[~,display_index_surgical_2] = sort(count_matrix_surgical_2(1,:),'descend');



%% Grab species to explore breakdowns of common skin species
species_to_show_figure = ["Cutibacterium acnes", "Staphylococcus aureus", "Staphylococcus epidermidis","Staphylococcus capitis"...
"Staphylococcus", "Corynebacterium kroppenstedtii", "__TAXCLUSTER3__Corynebacterium tuberculostearicum", "__TAXCLUSTER2__Corynebacterium accolens", "Corynebacterium jeikeium", "Corynebacterium amycolatum", "Corynebacterium"];

sorted_species_matrix_figure = get_specific_species_counts(species_to_show_figure,unique_species, species_counts_normalized);

%% Make bar plots breaking down genus distribution )(Figures 1C and 3B) 

% Figure 1C
fontsize_holder = 24; linewidth_holder = 2;

figure('Renderer', 'painters', 'Position', [10 10 1050 300])
colororder(all_colors_in_order_2())
bar(count_matrix_control_2(:,display_index_control_2)','stacked','EdgeColor','none')
set(gca,'xtick',[]); set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder)
ylim([0 1])
box off

figure('Renderer', 'painters', 'Position', [10 10 1050 300])
colororder(all_colors_in_order_2)
bar(count_matrix_surgical_2(:,display_index_surgical_2)','stacked','EdgeColor','none')
set(gca,'xtick',[]); set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder)
ylim([0 1])
box off

% Note how this uses count_matrix_control and not count_matrix_control_2,
% the second of which seperates out S. aureus 
wilranksum_cuti = signrank(count_matrix_control(1,:),count_matrix_surgical(1,:));
wilranksum_coryne = signrank(count_matrix_control(2,:),count_matrix_surgical(2,:));
wilranksum_staph = signrank(count_matrix_control(3,:),count_matrix_surgical(3,:));

%% Create one large scatter plot
fontsize_holder = 24; linewidth_holder = 2; 
box_plot_batch_nums = batch_num(matched_surgical_control_locations(:,1));
coryne_colors = [[1, 104, 111]./255; [246, 45, 0]./255;[1, 104, 111]./255;[246, 45, 0]./255];

sorted_genus_matrix_figure_no_cuti = sorted_genus_matrix_figure(2:end,:)./sum(sorted_genus_matrix_figure(2:end,:));
count_matrix_control_no_cuti = sorted_genus_matrix_figure_no_cuti(:,matched_surgical_control_locations(:,1));
count_matrix_surgical_no_cuti = sorted_genus_matrix_figure_no_cuti(:,matched_surgical_control_locations(:,2));

% Grab p-values for Corynebacterium and Staph when cutibacterium is removed 
% Note how this uses the non-S. aureus containing matrix
wilranksum_coryne_no_cuti = signrank(count_matrix_control_no_cuti(1,:),count_matrix_surgical_no_cuti(1,:));
wilranksum_staph_no_cuti = signrank(count_matrix_control_no_cuti(2,:),count_matrix_surgical_no_cuti(2,:));

sorted_species_matrix_figure_no_cuti = sorted_species_matrix_figure(2:end,:)./sum(sorted_species_matrix_figure(2:end,:));
count_matrix_control_species = sorted_species_matrix_figure_no_cuti(:,matched_surgical_control_locations(:,1));
count_matrix_surgical_species = sorted_species_matrix_figure_no_cuti(:,matched_surgical_control_locations(:,2));

figure('Renderer', 'painters', 'Position', [10 10 700 300])
scatter_batch(1, count_matrix_control_no_cuti(1,:), box_plot_batch_nums, coryne_colors(1,:), 40)
scatter_batch(2, count_matrix_surgical_no_cuti(1,:), box_plot_batch_nums, coryne_colors(2,:), 40)
scatter_batch(3.5, count_matrix_control_no_cuti(2,:), box_plot_batch_nums, coryne_colors(1,:), 40)
scatter_batch(4.5, count_matrix_surgical_no_cuti(2,:), box_plot_batch_nums, coryne_colors(2,:), 40)

scatter_batch(6+1, count_matrix_control_species(1,:), box_plot_batch_nums, coryne_colors(1,:), 40)
scatter_batch(7+1, count_matrix_surgical_species(1,:), box_plot_batch_nums, coryne_colors(2,:), 40)
scatter_batch(8.5+1, count_matrix_control_species(2,:), box_plot_batch_nums, coryne_colors(1,:), 40)
scatter_batch(9.5+1, count_matrix_surgical_species(2,:), box_plot_batch_nums, coryne_colors(2,:), 40)
scatter_batch(11+1, count_matrix_control_species(3,:), box_plot_batch_nums, coryne_colors(1,:), 40)
scatter_batch(12+1, count_matrix_surgical_species(3,:), box_plot_batch_nums, coryne_colors(2,:), 40)
box on;xlim([0,14]); ylim([0,1.25])
set(gca,'xtick',[]); set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder);set(gca,'ytick',0:.25:1); box off


%% Create box plots to demonstrate changes in species and genus 
fontsize_holder = 24; linewidth_holder = 2; 

sorted_species_matrix_figure_no_cuti = sorted_species_matrix_figure(2:end,:)./sum(sorted_species_matrix_figure(2:end,:));
wilranksum_aureus = signrank(sorted_species_matrix_figure_no_cuti(1,matched_surgical_control_locations(:,1)),sorted_species_matrix_figure_no_cuti(1,matched_surgical_control_locations(:,2)));
wilranksum_epidermidis = signrank(sorted_species_matrix_figure_no_cuti(2,matched_surgical_control_locations(:,1)),sorted_species_matrix_figure_no_cuti(2,matched_surgical_control_locations(:,2)));
wilranksum_capitis = signrank(sorted_species_matrix_figure_no_cuti(3,matched_surgical_control_locations(:,1)),sorted_species_matrix_figure_no_cuti(3,matched_surgical_control_locations(:,2)));


%% Now we explore Unifrac and Bray Curtis diversity metrics 

% Bray curtis
bray_curtis_matrix = f_braycurtis(condensed_counts_normalized(:,:));
[Y_braycurtis,BC_percent_variability] = cmdscale(bray_curtis_matrix);
BC_percent_variability = BC_percent_variability([1,2])./sum(BC_percent_variability);

control_bray_curt = Y_braycurtis(matched_surgical_control_locations(:,1),:);
surgical_bray_curt = Y_braycurtis(matched_surgical_control_locations(:,2),:);

fontsize_holder = 24; linewidth_holder = 2; 

% Figure 1A
figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(control_bray_curt(box_plot_batch_nums==1,1),control_bray_curt(box_plot_batch_nums==1,2),40,[1, 104, 111]./255,'o','filled');hold on
scatter(surgical_bray_curt(box_plot_batch_nums==1,1),surgical_bray_curt(box_plot_batch_nums==1,2),40,[246, 45, 0]./255,'o','filled');hold on
scatter(control_bray_curt(box_plot_batch_nums==2,1),control_bray_curt(box_plot_batch_nums==2,2),40,[1, 104, 111]./255,'v','filled');hold on
scatter(surgical_bray_curt(box_plot_batch_nums==2,1),surgical_bray_curt(box_plot_batch_nums==2,2),40,[246, 45, 0]./255,'v','filled')
scatter(control_bray_curt(box_plot_batch_nums==3,1),control_bray_curt(box_plot_batch_nums==3,2),40,[1, 104, 111]./255,'d','filled');hold on
scatter(surgical_bray_curt(box_plot_batch_nums==3,1),surgical_bray_curt(box_plot_batch_nums==3,2),40,[246, 45, 0]./255,'d','filled')
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); 
ylim([-.3,.5]);xlim([-.4,.4]);set(gca,'ytick',[-.3:.4:.5]); set(gca,'xtick',[-.4:.4:.4]);
box on; 

% Unifrac 
rawdata_unifrac = readtable('distance-matrix.xlsx');
unifrac = table2array(rawdata_unifrac(:,2:end));
unifrac_sample_name_index = cellfun(@(x) find(contains(new_samplenames, x(1:3))),rawdata_unifrac.Var1);
confirm_sample_order = sum(unifrac_sample_name_index' == 1:344)==344;
[Y_unifrac,Unifrac_percent_variability] = cmdscale(unifrac);
Uni_percent_variability = Unifrac_percent_variability([1,2])./sum(Unifrac_percent_variability);

control_unifrac = [Y_unifrac(matched_surgical_control_locations(:,1),1),Y_unifrac(matched_surgical_control_locations(:,1),2)];
surgical_unifrac = [Y_unifrac(matched_surgical_control_locations(:,2),1),Y_unifrac(matched_surgical_control_locations(:,2),2)];

% Figure 1A
figure('Renderer', 'painters', 'Position', [10 10 500 300]);
scatter(control_unifrac(box_plot_batch_nums==1,1),control_unifrac(box_plot_batch_nums==1,2),40,[1, 104, 111]./255,'o','filled');hold on
scatter(surgical_unifrac(box_plot_batch_nums==1,1),surgical_unifrac(box_plot_batch_nums==1,2),40,[246, 45, 0]./255,'o','filled');hold on
scatter(control_unifrac(box_plot_batch_nums==2,1),control_unifrac(box_plot_batch_nums==2,2),40,[1, 104, 111]./255,'v','filled');hold on
scatter(surgical_unifrac(box_plot_batch_nums==2,1),surgical_unifrac(box_plot_batch_nums==2,2),40,[246, 45, 0]./255,'v','filled')
scatter(control_unifrac(box_plot_batch_nums==3,1),control_unifrac(box_plot_batch_nums==3,2),40,[1, 104, 111]./255,'d','filled');hold on
scatter(surgical_unifrac(box_plot_batch_nums==3,1),surgical_unifrac(box_plot_batch_nums==3,2),40,[246, 45, 0]./255,'d','filled')
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); ylim([-.8,.6]);xlim([-.8,.6])
box on; set(gca,'xtick',[-.8:.7:.6]);set(gca,'ytick',[-.8:.7:.6]); 


%% Generating plots comparing surgical to control sample diversity metrics
bc_cont_surg = bray_curtis_matrix(matched_surgical_control_locations(:,1), matched_surgical_control_locations(:,2));
bc_cont_surg = bc_cont_surg(:);
bc_cont_cont = bray_curtis_matrix(matched_surgical_control_locations(:,1), matched_surgical_control_locations(:,1));
bc_cont_cont = bc_cont_cont(triu(ones(size(bc_cont_cont)),1)==1);
bc_surg_surg = bray_curtis_matrix(matched_surgical_control_locations(:,2), matched_surgical_control_locations(:,2));
bc_surg_surg = bc_surg_surg(triu(ones(size(bc_surg_surg)),1)==1);

bray_curt_beta_div = nan(length(bc_cont_surg),3);
bray_curt_beta_div(:,3) = bc_cont_surg;
bray_curt_beta_div(1:length(bc_surg_surg),2) = bc_surg_surg;
bray_curt_beta_div(1:length(bc_cont_cont),1) = bc_cont_cont;

unifrac_cont_surg = unifrac(matched_surgical_control_locations(:,1), matched_surgical_control_locations(:,2));
unifrac_cont_surg = unifrac_cont_surg(:);
unifrac_cont_cont = unifrac(matched_surgical_control_locations(:,1), matched_surgical_control_locations(:,1));
unifrac_cont_cont = unifrac_cont_cont(triu(ones(size(unifrac_cont_cont)),1)==1);
unifrac_surg_surg = unifrac(matched_surgical_control_locations(:,2), matched_surgical_control_locations(:,2));
unifrac_surg_surg = unifrac_surg_surg(triu(ones(size(unifrac_surg_surg)),1)==1);

unifrac_beta_div = nan(length(unifrac_cont_surg),3);
unifrac_beta_div(:,3) = unifrac_cont_surg;
unifrac_beta_div(1:length(unifrac_surg_surg),2) = unifrac_surg_surg;
unifrac_beta_div(1:length(unifrac_cont_cont),1) = unifrac_cont_cont;

% Figure 1B - Bray Curtis
bray_curt_beta_div_x = ones(length(bc_cont_surg),3); bray_curt_beta_div_x(:,2) =  bray_curt_beta_div_x(:,2) .*2; bray_curt_beta_div_x(:,3) =  bray_curt_beta_div_x(:,3) .*3;

figure('Renderer', 'painters', 'Position', [10 10 300 300])
b = boxchart(bray_curt_beta_div_x(:),bray_curt_beta_div(:),'GroupByColor',bray_curt_beta_div_x(:),'markerstyle', 'none','BoxWidth',2);
colororder([[1, 104, 111]./255;[246, 45, 0]./255;[.25 .25 .25]]); 
set(gca,'xtick',[]); set(gca,'fontsize',16); set(gca,'TickDir','out');set(gca,'linewidth',2);ylim([0 1.2]);set(gca,'ytick',0:.25:1);


% Figure 1C - Unifrac
unifrac_beta_div_x = ones(length(bc_cont_surg),3); unifrac_beta_div_x(:,2) =  unifrac_beta_div_x(:,2) .*2; unifrac_beta_div_x(:,3) =  unifrac_beta_div_x(:,3) .*3;

figure('Renderer', 'painters', 'Position', [10 10 300 300])
u = boxchart(unifrac_beta_div_x(:),unifrac_beta_div(:),'GroupByColor',unifrac_beta_div_x(:),'markerstyle', 'none','BoxWidth',2);
colororder([[1, 104, 111]./255;[246, 45, 0]./255;[.25 .25 .25]]); 
set(gca,'xtick',[]); set(gca,'fontsize',16); set(gca,'TickDir','out');set(gca,'linewidth',2);ylim([0 1.5]);set(gca,'ytick',0:.25:1.25);

bray_curt_C_S = ranksum(bray_curt_beta_div(:,1),bray_curt_beta_div(:,2)).*3;
bray_curt_C_CS = ranksum(bray_curt_beta_div(:,1),bray_curt_beta_div(:,3)).*3;
bray_curt_S_CS = ranksum(bray_curt_beta_div(:,2),bray_curt_beta_div(:,3)).*3;

unifrac_C_S = ranksum(unifrac_beta_div(:,1),unifrac_beta_div(:,2));
unifrac_C_CS = ranksum(unifrac_beta_div(:,1),unifrac_beta_div(:,3));
unifrac_S_CS = ranksum(unifrac_beta_div(:,2),unifrac_beta_div(:,3));

%% Create a sample summarizing all samples 

first_batch_ordering = cellfun(@(x) find(contains(table2array(metadata(:,1)),x)),samplenames_1);
second_batch_ordering = cellfun(@(x) find(contains(table2array(metadata(:,1)),x)),samplenames_2);
table_for_supplementary = metadata([first_batch_ordering,second_batch_ordering],:);
table_for_supplementary = [array2table(new_samplenames','VariableName',"Alternate Sample Name"),table_for_supplementary,array2table(~good_samples','VariableName',"Discarded from analysis"),array2table(sum(condensed_counts)','VariableName',"Cleaned Total Reads")];

cancer_metadata = readtable('cancer_type_metadata.xlsx'); 
new_samplenames_to_patient = cellfun(@(x) str2num(x(1:2)),new_samplenames);
index_cancer_meta_match_new_samp = arrayfun(@(x) find(cancer_metadata.Var1==x),new_samplenames_to_patient);

index_cancer_meta_match_new_samp(~contains(new_samplenames, '-X-')) = [];
cancer_type_holder = cell(length(new_samplenames),4);
cancer_type_holder(contains(new_samplenames, '-X-'),:) = table2cell(cancer_metadata(index_cancer_meta_match_new_samp,2:5));

table_for_supplementary = [table_for_supplementary, array2table(cancer_type_holder,'VariableName',["Tumor Type","Type of Closure","Unna Boot (Y/N)","SI Defect size (cm or cm2)"])];

writetable(table_for_supplementary,'sample_list_supplementary_table.xlsx');

%% Create table 2 which contains all p-values 
species_to_show_figure = ["Cutibacterium acnes", "Staphylococcus aureus", "Staphylococcus epidermidis","Staphylococcus capitis"...
"Staphylococcus", "Corynebacterium kroppenstedtii", "__TAXCLUSTER3__Corynebacterium tuberculostearicum", "__TAXCLUSTER2__Corynebacterium accolens", "Corynebacterium jeikeium", "Corynebacterium amycolatum", "Corynebacterium"];
species_to_show_legend = ["C. acnes", "S. aureus", "S. epidermidis","S. capitis","Staphylococcus"...
"C. kroppenstedtii", "TC3 C. tuberculostearicum","TC2 C. accolens", "C. jeikeium", "C. amycolatum","Corynebacterium"]; 

% Species table
species_p_val_holder = nan(length(species_to_show_figure),5);
for ii=1:length(species_to_show_figure)
    wil_sign_rank_spec = signrank(sorted_species_matrix_figure(ii,matched_surgical_control_locations(:,1)),sorted_species_matrix_figure(ii,matched_surgical_control_locations(:,2)));
    average_abun_control = mean(sorted_species_matrix_figure(ii,matched_surgical_control_locations(:,1)));
    average_abun_surgical = mean(sorted_species_matrix_figure(ii,matched_surgical_control_locations(:,2)));
    median_abun_control = median(sorted_species_matrix_figure(ii,matched_surgical_control_locations(:,1)));
    median_abun_surgical = median(sorted_species_matrix_figure(ii,matched_surgical_control_locations(:,2)));

    species_p_val_holder(ii,:) = [wil_sign_rank_spec,average_abun_control,average_abun_surgical,median_abun_control,median_abun_surgical];
end

% Genus table
genus_p_val_holder = nan(length(genuses_to_show_figure),5);
for ii=1:length(genuses_to_show_figure)
    wil_sign_rank_spec = signrank(sorted_genus_matrix_figure(ii,matched_surgical_control_locations(:,1)),sorted_genus_matrix_figure(ii,matched_surgical_control_locations(:,2)));
    average_abun_control = mean(sorted_genus_matrix_figure(ii,matched_surgical_control_locations(:,1)));
    average_abun_surgical = mean(sorted_genus_matrix_figure(ii,matched_surgical_control_locations(:,2)));
    median_abun_control = median(sorted_genus_matrix_figure(ii,matched_surgical_control_locations(:,1)));
    median_abun_surgical = median(sorted_genus_matrix_figure(ii,matched_surgical_control_locations(:,2)));

    genus_p_val_holder(ii,:) = [wil_sign_rank_spec,average_abun_control,average_abun_surgical,median_abun_control,median_abun_surgical];
end

% Cutibacterium removed genus
no_cuti_genus_p_val_holder = nan(length(genuses_to_show_figure)-1,5);
for ii=1:length(genuses_to_show_figure)-1
    wil_sign_rank_spec = signrank(sorted_genus_matrix_figure_no_cuti(ii,matched_surgical_control_locations(:,1)),sorted_genus_matrix_figure_no_cuti(ii,matched_surgical_control_locations(:,2)));
    average_abun_control = mean(sorted_genus_matrix_figure_no_cuti(ii,matched_surgical_control_locations(:,1)));
    average_abun_surgical = mean(sorted_genus_matrix_figure_no_cuti(ii,matched_surgical_control_locations(:,2)));
    median_abun_control = median(sorted_genus_matrix_figure_no_cuti(ii,matched_surgical_control_locations(:,1)));
    median_abun_surgical = median(sorted_genus_matrix_figure_no_cuti(ii,matched_surgical_control_locations(:,2)));

    no_cuti_genus_p_val_holder(ii,:) = [wil_sign_rank_spec,average_abun_control,average_abun_surgical,median_abun_control,median_abun_surgical];
end

% Cutibacterium removed species
sorted_species_matrix_no_cuti = sorted_species_matrix_figure(2:end,:)./sum(sorted_species_matrix_figure(2:end,:));
species_p_val_holder_no_cuti = nan(length(species_to_show_figure)-1,5);
for ii=1:length(species_to_show_figure)-1
    wil_sign_rank_spec = signrank(sorted_species_matrix_no_cuti(ii,matched_surgical_control_locations(:,1)),sorted_species_matrix_no_cuti(ii,matched_surgical_control_locations(:,2)));
    average_abun_control = mean(sorted_species_matrix_no_cuti(ii,matched_surgical_control_locations(:,1)));
    average_abun_surgical = mean(sorted_species_matrix_no_cuti(ii,matched_surgical_control_locations(:,2)));
    median_abun_control = median(sorted_species_matrix_no_cuti(ii,matched_surgical_control_locations(:,1)));
    median_abun_surgical = median(sorted_species_matrix_no_cuti(ii,matched_surgical_control_locations(:,2)));

    species_p_val_holder_no_cuti(ii,:) = [wil_sign_rank_spec,average_abun_control,average_abun_surgical,median_abun_control,median_abun_surgical];
end
%% Create corynebacterium abundance metrics (used for figure 3C)
coyrnebacterium_loc = find(contains(unique_species, 'Corynebacterium'));
coyrnebacterium_labels = unique_species(contains(unique_species, 'Corynebacterium'));
coyrnebacterium_control = mean(species_counts_normalized(coyrnebacterium_loc,matched_surgical_control_locations(:,1)),2);
corynebacterium_surgery = mean(species_counts_normalized(coyrnebacterium_loc,matched_surgical_control_locations(:,2)),2);

corynebacterium_abun_sum = coyrnebacterium_control + corynebacterium_surgery;
coyrnebacterium_species = coyrnebacterium_labels(corynebacterium_abun_sum>0);
control_sum_coryne = coyrnebacterium_control(corynebacterium_abun_sum>0);
surgical_sum_coryne = corynebacterium_surgery(corynebacterium_abun_sum>0);

%% Make table that begisn to summarize all metadata

first_batch_ordering = cellfun(@(x) find(contains(table2array(metadata(:,1)),x)),samplenames_1);
second_batch_ordering = cellfun(@(x) find(contains(table2array(metadata(:,1)),x)),samplenames_2);
table_for_supplementary = metadata([first_batch_ordering,second_batch_ordering],:);
table_for_supplementary = [array2table(new_samplenames','VariableName',"Alternate Sample Name"),table_for_supplementary,array2table(~good_samples','VariableName',"Discarded from analysis"),array2table(sum(condensed_counts)','VariableName',"Cleaned Total Reads")];

cancer_metadata = readtable('clinical_metadata.xlsx'); 
new_samplenames_to_patient = cellfun(@(x) str2num(x(1:2)),new_samplenames);
index_cancer_meta_match_new_samp = arrayfun(@(x) find(contains(string(cancer_metadata.SwabID),num2str(x))),new_samplenames_to_patient,'UniformOutput',false)';
index_cancer_meta_match_new_samp(cellfun(@isempty,index_cancer_meta_match_new_samp)) = {NaN};
index_cancer_meta_match_new_samp = cell2mat(index_cancer_meta_match_new_samp);

cancer_type_holder = cell(length(new_samplenames),9);
cancer_type_holder(:,:) = table2cell(cancer_metadata(index_cancer_meta_match_new_samp,[4:7,12:16]));

table_for_supplementary = [table_for_supplementary, array2table(cancer_type_holder,'VariableName',cancer_metadata.Properties.VariableNames([4:7,12:16]))];

%% Compare correlation between S. aureus in wounds and in day zero control skin

% Find any other control samples associated with the above patient
accompanying_control_samp_day_0 = cell(length(matched_surgical_control_locations),1);
accompanying_control_samp_day_7 = cell(length(matched_surgical_control_locations),1);

for ii=1:length(matched_surgical_control_locations)
    samp_to_parse = new_samplenames(matched_surgical_control_locations(ii,1)); samp_to_parse= samp_to_parse{:};
    holder_parsed = ones(size(new_samplenames)); holder_parsed(matched_surgical_control_locations(ii,1)) = 0;
    accompanying_control_samp_day_0(ii) = {find(contains(new_samplenames,samp_to_parse(1:2)) & good_samples & holder_parsed & (collection_time==0) & ~contains(new_samplenames, '-X-'))};
    accompanying_control_samp_day_7(ii) = {matched_surgical_control_locations(ii,2)};  
end

has_extra_d0_samples = cellfun(@length, accompanying_control_samp_day_0)>0;

staph_control = (sorted_species_matrix_figure(2,matched_surgical_control_locations(:,1))>0)';
staph_surgical = (sorted_species_matrix_figure(2,matched_surgical_control_locations(:,2))>0)';
has_d7_staph_cont = cellfun(@(x) sum(sorted_species_matrix_figure(2,x)>0)>0,accompanying_control_samp_day_7);
has_d0_staph_cont = cellfun(@(x) sum(sorted_species_matrix_figure(2,x)>0)>0,accompanying_control_samp_day_0);

d0_no_staph_d7_no_staph = sum(has_d0_staph_cont ==0 & staph_surgical==0 & has_extra_d0_samples);
d0_staph_d7_no_staph = sum(has_d0_staph_cont~=0 & staph_surgical==0 & has_extra_d0_samples);
d0_no_staph_d7_staph = sum(has_d0_staph_cont==0 & staph_surgical~=0 & has_extra_d0_samples);
d0_staph_d7_staph = sum(has_d0_staph_cont~=0 & staph_surgical~=0 & has_extra_d0_samples);

staph_fisher_setup = [d0_no_staph_d7_no_staph,d0_staph_d7_no_staph;d0_no_staph_d7_staph,d0_staph_d7_staph];
[~,p_staph_day0_pres_absense,~] = fishertest(staph_fisher_setup);

%% Looking at correlations between genuses - Supplementary figure 7
corynebacterium_control_values = genus_counts_normalized(contains(unique_genuses, 'Corynebacterium'),matched_surgical_control_locations(:,1));
corynebacterium_surgical_values = genus_counts_normalized(contains(unique_genuses, 'Corynebacterium'),matched_surgical_control_locations(:,2));

S_aureus_control_values = species_counts_normalized(contains(unique_species, "Staphylococcus aureus"),matched_surgical_control_locations(:,1));
S_aureus_surgical_values = species_counts_normalized(contains(unique_species, "Staphylococcus aureus"),matched_surgical_control_locations(:,2));

coryne_species_corr_wound = nan(size(unique_species));
aureus_genus_corr_wound = nan(size(unique_genuses));
coryne_species_corr_health = nan(size(unique_species));
aureus_genus_corr_health = nan(size(unique_genuses));


for ii=1:length(coryne_species_corr_wound)
    if contains(unique_species(ii), "Corynebacterium")
        continue
    else
        corr_cov_mat= corrcoef(corynebacterium_surgical_values, species_counts_normalized(ii,matched_surgical_control_locations(:,2)));
        coryne_species_corr_wound(ii) = corr_cov_mat(2);
    end
end

for ii=1:length(coryne_species_corr_health)
    if contains(unique_species(ii), "Corynebacterium")
        continue
    else
        corr_cov_mat= corrcoef(corynebacterium_control_values, species_counts_normalized(ii,matched_surgical_control_locations(:,1)));
        coryne_species_corr_health(ii) = corr_cov_mat(2);
    end
end

for ii=1:length(aureus_genus_corr_wound)
    if contains(unique_genuses(ii), "Staphylococcus")
        continue
    else
        corr_cov_mat= corrcoef(S_aureus_surgical_values, genus_counts_normalized(ii,matched_surgical_control_locations(:,2)));
        aureus_genus_corr_wound(ii) = corr_cov_mat(2);
    end
end

for ii=1:length(aureus_genus_corr_health)
    if contains(unique_genuses(ii), "Staphylococcus")
        continue
    else
        corr_cov_mat= corrcoef(S_aureus_control_values, genus_counts_normalized(ii,matched_surgical_control_locations(:,1)));
        aureus_genus_corr_health(ii) = corr_cov_mat(2);
    end
end

figure; histogram(coryne_species_corr_wound,40,'Facecolor',[.5,.5,.5]);
set(gca,'fontsize',16); set(gca,'TickDir','out'); set(gca,'linewidth',1); box off; xlim([-.7,.7]);

figure; histogram(coryne_species_corr_health,40,'Facecolor',[.5,.5,.5]);
set(gca,'fontsize',16); set(gca,'TickDir','out'); set(gca,'linewidth',1); box off; xlim([-.7,.7]);

figure; histogram(aureus_genus_corr_wound,40,'Facecolor',[.5,.5,.5]);
set(gca,'fontsize',16); set(gca,'TickDir','out'); set(gca,'linewidth',1); box off; xlim([-.7,.7]);

figure; histogram(aureus_genus_corr_health,40,'Facecolor',[.5,.5,.5]);
set(gca,'fontsize',16); set(gca,'TickDir','out'); set(gca,'linewidth',1); box off; xlim([-.7,.7]);

coryne_aureus_wound_corr = coryne_species_corr_wound(contains(unique_species, "Staphylococcus aureus"));
coryne_aureus_health_corr = coryne_species_corr_health(contains(unique_species, "Staphylococcus aureus"));
aureus_coryne_wound_corr = aureus_genus_corr_wound(contains(unique_genuses, "Corynebacterium"));
aureus_coryne_health_corr = aureus_genus_corr_health(contains(unique_genuses, "Corynebacterium"));

[~,min_aureus_coryne_health_index] = min(aureus_genus_corr_health);
[~,min_coryne_aureus_health_index] = min(coryne_species_corr_health);

aureus_coryne_max_neg_health_corr = unique_genuses(min_aureus_coryne_health_index);
coryne_aureus_max_neg_health_corr = unique_species(min_coryne_aureus_health_index);

%% Generate scatter plot of S. aureus vs Corynebacterium

figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(corynebacterium_surgical_values,S_aureus_surgical_values,[],'k','filled');
hold on
plot([0,1],[1,0],'k');
set(gca,'fontsize',18); set(gca,'TickDir','out'); set(gca,'linewidth',1); box off; yticks([0:.25:1]);xticks([0:.25:1]);


corrcoef_coryne_aureus = corrcoef(corynebacterium_surgical_values,S_aureus_surgical_values);

%% Gender
gender = string(table_for_supplementary.Gender);
gender_control = (gender(matched_surgical_control_locations(:,1))=="F")';
gender_surgical = (gender(matched_surgical_control_locations(:,2))=="M")';
fontsize_holder = 18;
linewidth_holder = 2; 

% Bray curtis

figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(control_bray_curt(box_plot_batch_nums==1 & gender_control==1,1),control_bray_curt(box_plot_batch_nums==1 & gender_control==1,2),35,[1, 104, 111]./255,'o','filled');hold on
scatter(control_bray_curt(box_plot_batch_nums==2 & gender_control==1,1),control_bray_curt(box_plot_batch_nums==2 & gender_control==1,2),35,[1, 104, 111]./255,'v','filled');hold on
scatter(control_bray_curt(box_plot_batch_nums==1 & gender_control==0,1),control_bray_curt(box_plot_batch_nums==1 & gender_control==0,2),35,[1, 104, 111]./255,'o','LineWidth',1.5);hold on
scatter(control_bray_curt(box_plot_batch_nums==2 & gender_control==0,1),control_bray_curt(box_plot_batch_nums==2 & gender_control==0,2),35,[1, 104, 111]./255,'v','LineWidth',1.5);hold on
scatter(control_bray_curt(box_plot_batch_nums==3 & gender_control==0,1),control_bray_curt(box_plot_batch_nums==3 & gender_control==0,2),35,[1, 104, 111]./255,'d','LineWidth',1.5);hold on
scatter(control_bray_curt(box_plot_batch_nums==3 & gender_control==1,1),control_bray_curt(box_plot_batch_nums==3 & gender_control==1,2),35,[1, 104, 111]./255,'d','filled');hold on
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); 
ylim([-.3,.5]);xlim([-.4,.4]);set(gca,'ytick',[-.3:.4:.5]); set(gca,'xtick',[-.4:.4:.4]);box on; 

figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(surgical_bray_curt(box_plot_batch_nums==2 & gender_control==1,1),surgical_bray_curt(box_plot_batch_nums==2 & gender_control==1,2),35,[246, 45, 0]./255,'v','filled'); hold on
scatter(surgical_bray_curt(box_plot_batch_nums==1 & gender_control==1,1),surgical_bray_curt(box_plot_batch_nums==1 & gender_control==1,2),35,[246, 45, 0]./255,'o','filled');hold on
scatter(surgical_bray_curt(box_plot_batch_nums==2 & gender_control==0,1),surgical_bray_curt(box_plot_batch_nums==2 & gender_control==0,2),35,[246, 45, 0]./255,'v','LineWidth',1.5); hold on
scatter(surgical_bray_curt(box_plot_batch_nums==1 & gender_control==0,1),surgical_bray_curt(box_plot_batch_nums==1 & gender_control==0,2),35,[246, 45, 0]./255,'o','LineWidth',1.5);hold on
scatter(surgical_bray_curt(box_plot_batch_nums==3 & gender_control==1,1),surgical_bray_curt(box_plot_batch_nums==3 & gender_control==1,2),35,[246, 45, 0]./255,'d','filled');hold on
scatter(surgical_bray_curt(box_plot_batch_nums==3 & gender_control==0,1),surgical_bray_curt(box_plot_batch_nums==3 & gender_control==0,2),35,[246, 45, 0]./255,'d','LineWidth',1.5);hold on
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); 
ylim([-.3,.5]);xlim([-.4,.4]);set(gca,'ytick',[-.3:.4:.5]); set(gca,'xtick',[-.4:.4:.4]);box on; 

% Unifrac

figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(control_unifrac(box_plot_batch_nums==1 & gender_control==1,1),control_unifrac(box_plot_batch_nums==1 & gender_control==1,2),35,[1, 104, 111]./255,'o','filled');hold on
scatter(control_unifrac(box_plot_batch_nums==2 & gender_control==1,1),control_unifrac(box_plot_batch_nums==2 & gender_control==1,2),35,[1, 104, 111]./255,'v','filled');hold on
scatter(control_unifrac(box_plot_batch_nums==1 & gender_control==0,1),control_unifrac(box_plot_batch_nums==1 & gender_control==0,2),35,[1, 104, 111]./255,'o','LineWidth',1.5);hold on
scatter(control_unifrac(box_plot_batch_nums==2 & gender_control==0,1),control_unifrac(box_plot_batch_nums==2 & gender_control==0,2),35,[1, 104, 111]./255,'v','LineWidth',1.5);hold on
scatter(control_unifrac(box_plot_batch_nums==3 & gender_control==0,1),control_unifrac(box_plot_batch_nums==3 & gender_control==0,2),35,[1, 104, 111]./255,'d','LineWidth',1.5);hold on
scatter(control_unifrac(box_plot_batch_nums==3 & gender_control==1,1),control_unifrac(box_plot_batch_nums==3 & gender_control==1,2),35,[1, 104, 111]./255,'d','filled');hold on
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); ylim([-.8,.6]);xlim([-.8,.6])
box on; set(gca,'xtick',[-.8:.7:.6]);set(gca,'ytick',[-.8:.7:.6]); 

figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(surgical_unifrac(box_plot_batch_nums==2 & gender_control==1,1),surgical_unifrac(box_plot_batch_nums==2 & gender_control==1,2),35,[246, 45, 0]./255,'v','filled'); hold on
scatter(surgical_unifrac(box_plot_batch_nums==1 & gender_control==1,1),surgical_unifrac(box_plot_batch_nums==1 & gender_control==1,2),35,[246, 45, 0]./255,'o','filled');hold on
scatter(surgical_unifrac(box_plot_batch_nums==2 & gender_control==0,1),surgical_unifrac(box_plot_batch_nums==2 & gender_control==0,2),35,[246, 45, 0]./255,'v','LineWidth',1.5); hold on
scatter(surgical_unifrac(box_plot_batch_nums==1 & gender_control==0,1),surgical_unifrac(box_plot_batch_nums==1 & gender_control==0,2),35,[246, 45, 0]./255,'o','LineWidth',1.5);hold on
scatter(surgical_unifrac(box_plot_batch_nums==3 & gender_control==1,1),surgical_unifrac(box_plot_batch_nums==3 & gender_control==1,2),35,[246, 45, 0]./255,'d','filled');hold on
scatter(surgical_unifrac(box_plot_batch_nums==3 & gender_control==0,1),surgical_unifrac(box_plot_batch_nums==3 & gender_control==0,2),35,[246, 45, 0]./255,'d','LineWidth',1.5);hold on
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); ylim([-.8,.6]);xlim([-.8,.6])
box on; set(gca,'xtick',[-.8:.7:.6]);set(gca,'ytick',[-.8:.7:.6]); 

%% Cancer type 
color_1 = [166,206,227]./255;
color_2 = [31,120,180]./255;
color_3 = [33, 176, 57]./255;

cancertype = string(table_for_supplementary.TumorType_Simplified);
cancer_log = double(cancertype(matched_surgical_control_locations(:,1))=="BCC")';
cancer_log(cancertype(matched_surgical_control_locations(:,1))=="SCC") = 0;
cancer_log(cancertype(matched_surgical_control_locations(:,1))=="Other") = 2;

figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(control_bray_curt(box_plot_batch_nums==1 & cancer_log==1,1),control_bray_curt(box_plot_batch_nums==1 & cancer_log==1,2),35,color_1,'o','filled');hold on
scatter(control_bray_curt(box_plot_batch_nums==2 & cancer_log==1,1),control_bray_curt(box_plot_batch_nums==2 & cancer_log==1,2),35,color_1,'v','filled');hold on
scatter(control_bray_curt(box_plot_batch_nums==1 & cancer_log==0,1),control_bray_curt(box_plot_batch_nums==1 & cancer_log==0,2),35,color_2,'o','filled');hold on
scatter(control_bray_curt(box_plot_batch_nums==2 & cancer_log==0,1),control_bray_curt(box_plot_batch_nums==2 & cancer_log==0,2),35,color_2,'v','filled');hold on
scatter(control_bray_curt(box_plot_batch_nums==1 & cancer_log==2,1),control_bray_curt(box_plot_batch_nums==1 & cancer_log==2,2),35,color_3,'o','filled');hold on
scatter(control_bray_curt(box_plot_batch_nums==2 & cancer_log==2,1),control_bray_curt(box_plot_batch_nums==2 & cancer_log==2,2),35,color_3,'v','filled');hold on
scatter(control_bray_curt(box_plot_batch_nums==3 & cancer_log==1,1),control_bray_curt(box_plot_batch_nums==3 & cancer_log==1,2),35,color_1,'d','filled');hold on
scatter(control_bray_curt(box_plot_batch_nums==3 & cancer_log==0,1),control_bray_curt(box_plot_batch_nums==3 & cancer_log==0,2),35,color_2,'d','filled');hold on
scatter(control_bray_curt(box_plot_batch_nums==3 & cancer_log==2,1),control_bray_curt(box_plot_batch_nums==3 & cancer_log==2,2),35,color_3,'d','filled');hold on
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); 
ylim([-.3,.5]);xlim([-.4,.4]);set(gca,'ytick',[-.3:.4:.5]); set(gca,'xtick',[-.4:.4:.4]);box on;  

figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(surgical_bray_curt(box_plot_batch_nums==2 & cancer_log==1,1),surgical_bray_curt(box_plot_batch_nums==2 & cancer_log==1,2),35,color_1,'v','filled'); hold on
scatter(surgical_bray_curt(box_plot_batch_nums==1 & cancer_log==1,1),surgical_bray_curt(box_plot_batch_nums==1 & cancer_log==1,2),35,color_1,'o','filled');hold on
scatter(surgical_bray_curt(box_plot_batch_nums==2 & cancer_log==0,1),surgical_bray_curt(box_plot_batch_nums==2 & cancer_log==0,2),35,color_2,'v','filled'); hold on
scatter(surgical_bray_curt(box_plot_batch_nums==1 & cancer_log==0,1),surgical_bray_curt(box_plot_batch_nums==1 & cancer_log==0,2),35,color_2,'o','filled');hold on
scatter(surgical_bray_curt(box_plot_batch_nums==2 & cancer_log==2,1),surgical_bray_curt(box_plot_batch_nums==2 & cancer_log==2,2),35,color_3,'v','filled'); hold on
scatter(surgical_bray_curt(box_plot_batch_nums==1 & cancer_log==2,1),surgical_bray_curt(box_plot_batch_nums==1 & cancer_log==2,2),35,color_3,'o','filled');hold on
scatter(surgical_bray_curt(box_plot_batch_nums==3 & cancer_log==1,1),surgical_bray_curt(box_plot_batch_nums==3 & cancer_log==1,2),35,color_1,'d','filled'); hold on
scatter(surgical_bray_curt(box_plot_batch_nums==3 & cancer_log==0,1),surgical_bray_curt(box_plot_batch_nums==3 & cancer_log==0,2),35,color_2,'d','filled'); hold on
scatter(surgical_bray_curt(box_plot_batch_nums==3 & cancer_log==2,1),surgical_bray_curt(box_plot_batch_nums==3 & cancer_log==2,2),35,color_3,'d','filled'); hold on
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); 
ylim([-.3,.5]);xlim([-.4,.4]);set(gca,'ytick',[-.3:.4:.5]); set(gca,'xtick',[-.4:.4:.4]);box on; 

figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(control_unifrac(box_plot_batch_nums==1 & cancer_log==1,1),control_unifrac(box_plot_batch_nums==1 & cancer_log==1,2),35,color_1,'o','filled');hold on
scatter(control_unifrac(box_plot_batch_nums==2 & cancer_log==1,1),control_unifrac(box_plot_batch_nums==2 & cancer_log==1,2),35,color_1,'v','filled');hold on
scatter(control_unifrac(box_plot_batch_nums==1 & cancer_log==0,1),control_unifrac(box_plot_batch_nums==1 & cancer_log==0,2),35,color_2,'o','filled');hold on
scatter(control_unifrac(box_plot_batch_nums==2 & cancer_log==0,1),control_unifrac(box_plot_batch_nums==2 & cancer_log==0,2),35,color_2,'v','filled');hold on
scatter(control_unifrac(box_plot_batch_nums==1 & cancer_log==2,1),control_unifrac(box_plot_batch_nums==1 & cancer_log==2,2),35,color_3,'o','filled');hold on
scatter(control_unifrac(box_plot_batch_nums==2 & cancer_log==2,1),control_unifrac(box_plot_batch_nums==2 & cancer_log==2,2),35,color_3,'v','filled');hold on
scatter(control_unifrac(box_plot_batch_nums==3 & cancer_log==1,1),control_unifrac(box_plot_batch_nums==3 & cancer_log==1,2),35,color_1,'d','filled');hold on
scatter(control_unifrac(box_plot_batch_nums==3 & cancer_log==0,1),control_unifrac(box_plot_batch_nums==3 & cancer_log==0,2),35,color_2,'d','filled');hold on
scatter(control_unifrac(box_plot_batch_nums==3 & cancer_log==2,1),control_unifrac(box_plot_batch_nums==3 & cancer_log==2,2),35,color_3,'d','filled');hold on
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); ylim([-.8,.6]);xlim([-.8,.6])
box on; set(gca,'xtick',[-.8:.7:.6]);set(gca,'ytick',[-.8:.7:.6]); 

figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(surgical_unifrac(box_plot_batch_nums==2 & cancer_log==1,1),surgical_unifrac(box_plot_batch_nums==2 & cancer_log==1,2),35,color_1,'v','filled'); hold on
scatter(surgical_unifrac(box_plot_batch_nums==1 & cancer_log==1,1),surgical_unifrac(box_plot_batch_nums==1 & cancer_log==1,2),35,color_1,'o','filled');hold on
scatter(surgical_unifrac(box_plot_batch_nums==2 & cancer_log==0,1),surgical_unifrac(box_plot_batch_nums==2 & cancer_log==0,2),35,color_2,'v','filled'); hold on
scatter(surgical_unifrac(box_plot_batch_nums==1 & cancer_log==0,1),surgical_unifrac(box_plot_batch_nums==1 & cancer_log==0,2),35,color_2,'o','filled');hold on
scatter(surgical_unifrac(box_plot_batch_nums==2 & cancer_log==2,1),surgical_unifrac(box_plot_batch_nums==2 & cancer_log==2,2),35,color_3,'v','filled'); hold on
scatter(surgical_unifrac(box_plot_batch_nums==1 & cancer_log==2,1),surgical_unifrac(box_plot_batch_nums==1 & cancer_log==2,2),35,color_3,'o','filled');hold on
scatter(surgical_unifrac(box_plot_batch_nums==3 & cancer_log==1,1),surgical_unifrac(box_plot_batch_nums==3 & cancer_log==1,2),35,color_1,'d','filled');hold on
scatter(surgical_unifrac(box_plot_batch_nums==3 & cancer_log==0,1),surgical_unifrac(box_plot_batch_nums==3 & cancer_log==0,2),35,color_2,'d','filled');hold on
scatter(surgical_unifrac(box_plot_batch_nums==3 & cancer_log==2,1),surgical_unifrac(box_plot_batch_nums==3 & cancer_log==2,2),35,color_3,'d','filled');hold on
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); ylim([-.8,.6]);xlim([-.8,.6])
box on; set(gca,'xtick',[-.8:.7:.6]);set(gca,'ytick',[-.8:.7:.6]); 

%% Wound closure type

closure_type = string(table_for_supplementary.CompleteVsPartialSIH);
closure_type_log = (closure_type(matched_surgical_control_locations(:,1))=="Partial" | contains(closure_type(matched_surgical_control_locations(:,1)),"Complex"))';

figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(control_bray_curt(box_plot_batch_nums==1 & closure_type_log==1,1),control_bray_curt(box_plot_batch_nums==1 & closure_type_log==1,2),35,[1, 104, 111]./255,'o','filled');hold on
scatter(control_bray_curt(box_plot_batch_nums==2 & closure_type_log==1,1),control_bray_curt(box_plot_batch_nums==2 & closure_type_log==1,2),35,[1, 104, 111]./255,'v','filled');hold on
scatter(control_bray_curt(box_plot_batch_nums==1 & closure_type_log==0,1),control_bray_curt(box_plot_batch_nums==1 & closure_type_log==0,2),35,[1, 104, 111]./255,'o','LineWidth',1.5);hold on
scatter(control_bray_curt(box_plot_batch_nums==2 & closure_type_log==0,1),control_bray_curt(box_plot_batch_nums==2 & closure_type_log==0,2),35,[1, 104, 111]./255,'v','LineWidth',1.5);hold on
scatter(control_bray_curt(box_plot_batch_nums==3 & closure_type_log==1,1),control_bray_curt(box_plot_batch_nums==3 & closure_type_log==1,2),35,[1, 104, 111]./255,'d','filled');hold on
scatter(control_bray_curt(box_plot_batch_nums==3 & closure_type_log==0,1),control_bray_curt(box_plot_batch_nums==3 & closure_type_log==0,2),35,[1, 104, 111]./255,'d','LineWidth',1.5);hold on
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); 
ylim([-.3,.5]);xlim([-.4,.4]);set(gca,'ytick',[-.3:.4:.5]); set(gca,'xtick',[-.4:.4:.4]);box on; 

figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(surgical_bray_curt(box_plot_batch_nums==2 & closure_type_log==1,1),surgical_bray_curt(box_plot_batch_nums==2 & closure_type_log==1,2),35,[246, 45, 0]./255,'v','filled'); hold on
scatter(surgical_bray_curt(box_plot_batch_nums==1 & closure_type_log==1,1),surgical_bray_curt(box_plot_batch_nums==1 & closure_type_log==1,2),35,[246, 45, 0]./255,'o','filled');hold on
scatter(surgical_bray_curt(box_plot_batch_nums==2 & closure_type_log==0,1),surgical_bray_curt(box_plot_batch_nums==2 & closure_type_log==0,2),35,[246, 45, 0]./255,'v','LineWidth',1.5); hold on
scatter(surgical_bray_curt(box_plot_batch_nums==1 & closure_type_log==0,1),surgical_bray_curt(box_plot_batch_nums==1 & closure_type_log==0,2),35,[246, 45, 0]./255,'o','LineWidth',1.5);hold on
scatter(surgical_bray_curt(box_plot_batch_nums==3 & closure_type_log==1,1),surgical_bray_curt(box_plot_batch_nums==3 & closure_type_log==1,2),35,[246, 45, 0]./255,'d','filled'); hold on
scatter(surgical_bray_curt(box_plot_batch_nums==3 & closure_type_log==0,1),surgical_bray_curt(box_plot_batch_nums==3 & closure_type_log==0,2),35,[246, 45, 0]./255,'d','LineWidth',1.5); hold on
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); 
ylim([-.3,.5]);xlim([-.4,.4]);set(gca,'ytick',[-.3:.4:.5]); set(gca,'xtick',[-.4:.4:.4]);box on;  

figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(control_unifrac(box_plot_batch_nums==1 & closure_type_log==1,1),control_unifrac(box_plot_batch_nums==1 & closure_type_log==1,2),35,[1, 104, 111]./255,'o','filled');hold on
scatter(control_unifrac(box_plot_batch_nums==2 & closure_type_log==1,1),control_unifrac(box_plot_batch_nums==2 & closure_type_log==1,2),35,[1, 104, 111]./255,'v','filled');hold on
scatter(control_unifrac(box_plot_batch_nums==1 & closure_type_log==0,1),control_unifrac(box_plot_batch_nums==1 & closure_type_log==0,2),35,[1, 104, 111]./255,'o','LineWidth',1.5);hold on
scatter(control_unifrac(box_plot_batch_nums==2 & closure_type_log==0,1),control_unifrac(box_plot_batch_nums==2 & closure_type_log==0,2),35,[1, 104, 111]./255,'v','LineWidth',1.5);hold on
scatter(control_unifrac(box_plot_batch_nums==3 & closure_type_log==1,1),control_unifrac(box_plot_batch_nums==3 & closure_type_log==1,2),35,[1, 104, 111]./255,'d','filled');hold on
scatter(control_unifrac(box_plot_batch_nums==3 & closure_type_log==0,1),control_unifrac(box_plot_batch_nums==3 & closure_type_log==0,2),35,[1, 104, 111]./255,'d','LineWidth',1.5);hold on
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); ylim([-.8,.6]);xlim([-.8,.6])
box on; set(gca,'xtick',[-.8:.7:.6]);set(gca,'ytick',[-.8:.7:.6]); 

figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(surgical_unifrac(box_plot_batch_nums==2 & closure_type_log==1,1),surgical_unifrac(box_plot_batch_nums==2 & closure_type_log==1,2),35,[246, 45, 0]./255,'v','filled'); hold on
scatter(surgical_unifrac(box_plot_batch_nums==1 & closure_type_log==1,1),surgical_unifrac(box_plot_batch_nums==1 & closure_type_log==1,2),35,[246, 45, 0]./255,'o','filled');hold on
scatter(surgical_unifrac(box_plot_batch_nums==2 & closure_type_log==0,1),surgical_unifrac(box_plot_batch_nums==2 & closure_type_log==0,2),35,[246, 45, 0]./255,'v','LineWidth',1.5); hold on
scatter(surgical_unifrac(box_plot_batch_nums==1 & closure_type_log==0,1),surgical_unifrac(box_plot_batch_nums==1 & closure_type_log==0,2),35,[246, 45, 0]./255,'o','LineWidth',1.5);hold on
scatter(surgical_unifrac(box_plot_batch_nums==3 & closure_type_log==0,1),surgical_unifrac(box_plot_batch_nums==3 & closure_type_log==0,2),35,[246, 45, 0]./255,'d','LineWidth',1.5); hold on
scatter(surgical_unifrac(box_plot_batch_nums==3 & closure_type_log==1,1),surgical_unifrac(box_plot_batch_nums==3 & closure_type_log==1,2),35,[246, 45, 0]./255,'d','filled'); hold on
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); ylim([-.8,.6]);xlim([-.8,.6])
box on; set(gca,'xtick',[-.8:.7:.6]);set(gca,'ytick',[-.8:.7:.6]); 

%% Wound closure size
closure_size = string(table_for_supplementary.SIHArea_mm2_);
closure_size_log = double(closure_size(matched_surgical_control_locations(:,1)))';
closure_size_log_round = round(log(closure_size_log)*100);closure_size_log_round(isnan(closure_size_log_round))=max(closure_size_log_round)+1;
color_1 = 255.*[1,0,1]./255;
color_2 = 255.*[0,1,1]./255;
grad_num = max(closure_size_log_round)-1;
colors_grad = [linspace(color_1(1),color_2(1),grad_num)', linspace(color_1(2),color_2(2),grad_num)', linspace(color_1(3),color_2(3),grad_num)'];
colors_grad = [colors_grad;[.5,.5,.5]];
colors_to_plot = colors_grad(closure_size_log_round,:);

figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(control_bray_curt(box_plot_batch_nums==1 ,1),control_bray_curt(box_plot_batch_nums==1,2),35,colors_to_plot(box_plot_batch_nums==1,:),'o','filled');hold on
scatter(control_bray_curt(box_plot_batch_nums==2,1),control_bray_curt(box_plot_batch_nums==2,2),35,colors_to_plot(box_plot_batch_nums==2,:),'v','filled');hold on
scatter(control_bray_curt(box_plot_batch_nums==3,1),control_bray_curt(box_plot_batch_nums==3,2),35,colors_to_plot(box_plot_batch_nums==3,:),'d','filled');hold on
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); 
ylim([-.3,.5]);xlim([-.4,.4]);set(gca,'ytick',[-.3:.4:.5]); set(gca,'xtick',[-.4:.4:.4]);box on; 

figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(surgical_bray_curt(box_plot_batch_nums==1 ,1),surgical_bray_curt(box_plot_batch_nums==1,2),35,colors_to_plot(box_plot_batch_nums==1,:),'o','filled');hold on
scatter(surgical_bray_curt(box_plot_batch_nums==2,1),surgical_bray_curt(box_plot_batch_nums==2,2),35,colors_to_plot(box_plot_batch_nums==2,:),'v','filled');hold on
scatter(surgical_bray_curt(box_plot_batch_nums==3,1),surgical_bray_curt(box_plot_batch_nums==3,2),35,colors_to_plot(box_plot_batch_nums==3,:),'d','filled');hold on
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); 
ylim([-.3,.5]);xlim([-.4,.4]);set(gca,'ytick',[-.3:.4:.5]); set(gca,'xtick',[-.4:.4:.4]);box on; 

figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(control_unifrac(box_plot_batch_nums==1 ,1),control_unifrac(box_plot_batch_nums==1,2),35,colors_to_plot(box_plot_batch_nums==1,:),'o','filled');hold on
scatter(control_unifrac(box_plot_batch_nums==2,1),control_unifrac(box_plot_batch_nums==2,2),35,colors_to_plot(box_plot_batch_nums==2,:),'v','filled');hold on
scatter(control_unifrac(box_plot_batch_nums==3,1),control_unifrac(box_plot_batch_nums==3,2),35,colors_to_plot(box_plot_batch_nums==3,:),'d','filled');hold on
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); ylim([-.8,.6]);xlim([-.8,.6])
box on; set(gca,'xtick',[-.8:.7:.6]);set(gca,'ytick',[-.8:.7:.6]); 

figure('Renderer', 'painters', 'Position', [10 10 500 300])
scatter(surgical_unifrac(box_plot_batch_nums==1 ,1),surgical_unifrac(box_plot_batch_nums==1,2),35,colors_to_plot(box_plot_batch_nums==1,:),'o','filled');hold on
scatter(surgical_unifrac(box_plot_batch_nums==2,1),surgical_unifrac(box_plot_batch_nums==2,2),35,colors_to_plot(box_plot_batch_nums==2,:),'v','filled');hold on
scatter(surgical_unifrac(box_plot_batch_nums==3,1),surgical_unifrac(box_plot_batch_nums==3,2),35,colors_to_plot(box_plot_batch_nums==3,:),'d','filled');hold on
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); ylim([-.8,.6]);xlim([-.8,.6])
box on; set(gca,'xtick',[-.8:.7:.6]);set(gca,'ytick',[-.8:.7:.6]); 

%% Seperating samples by location
loc_string = string(table_for_supplementary.Location);
locations_to_parse_for = {["eyelid","canthus"],["postauricular","scalp"],["ear","helix"],["nasal","nose","alar crease"],"lip",["forehead","temple","glabella"],["cheek","chin"],"neck","shin"};
samplenames_control = new_samplenames(matched_surgical_control_locations(:,1));
samplenames_surgical = new_samplenames(matched_surgical_control_locations(:,2));

control_loc_holder = zeros(size(matched_surgical_control_locations(:,1)));
surgical_loc_holder = zeros(size(matched_surgical_control_locations(:,1)));

number_of_loc_occurances = nan(9,1);
for ii = 1:length(locations_to_parse_for)
    control_loc = cellfun(@(x) contains(x,locations_to_parse_for{ii}), samplenames_control);
    surg_loc = cellfun(@(x) contains(x,locations_to_parse_for{ii}), samplenames_surgical);
    if sum(control_loc'>0 & control_loc_holder)>0
        disp("WARNING")
        break
    end
    control_loc_holder(control_loc) = ii;
    surgical_loc_holder(surg_loc) = ii;
    number_of_loc_occurances(ii) = sum(control_loc);
end

figure('Renderer', 'painters', 'Position', [10 10 500 300])
for ii=1:9
    scatter(control_bray_curt(box_plot_batch_nums==1 & control_loc_holder'==ii,1),control_bray_curt(box_plot_batch_nums==1 & control_loc_holder'==ii,2),35,all_colors_in_order(ii,:),'o','filled');hold on
    scatter(control_bray_curt(box_plot_batch_nums==2 & control_loc_holder'==ii,1),control_bray_curt(box_plot_batch_nums==2 & control_loc_holder'==ii,2),35,all_colors_in_order(ii,:),'v','filled');hold on
    scatter(control_bray_curt(box_plot_batch_nums==3 & control_loc_holder'==ii,1),control_bray_curt(box_plot_batch_nums==3 & control_loc_holder'==ii,2),35,all_colors_in_order(ii,:),'d','filled');hold on

end
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); 
ylim([-.3,.5]);xlim([-.4,.4]);set(gca,'ytick',[-.3:.4:.5]); set(gca,'xtick',[-.4:.4:.4]);box on; 

figure('Renderer', 'painters', 'Position', [10 10 500 300])
for ii=1:9
    scatter(surgical_bray_curt(box_plot_batch_nums==1 & surgical_loc_holder'==ii,1),surgical_bray_curt(box_plot_batch_nums==1 & surgical_loc_holder'==ii,2),35,all_colors_in_order(ii,:),'o','filled');hold on
    scatter(surgical_bray_curt(box_plot_batch_nums==2 & surgical_loc_holder'==ii,1),surgical_bray_curt(box_plot_batch_nums==2 & surgical_loc_holder'==ii,2),35,all_colors_in_order(ii,:),'v','filled');hold on
    scatter(surgical_bray_curt(box_plot_batch_nums==3 & surgical_loc_holder'==ii,1),surgical_bray_curt(box_plot_batch_nums==3 & surgical_loc_holder'==ii,2),35,all_colors_in_order(ii,:),'d','filled');hold on

end
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); 
ylim([-.3,.5]);xlim([-.4,.4]);set(gca,'ytick',[-.3:.4:.5]); set(gca,'xtick',[-.4:.4:.4]);box on;  

figure('Renderer', 'painters', 'Position', [10 10 500 300])
for ii=1:9
    scatter(control_unifrac(box_plot_batch_nums==1 & control_loc_holder'==ii,1),control_unifrac(box_plot_batch_nums==1 & control_loc_holder'==ii,2),35,all_colors_in_order(ii,:),'o','filled');hold on
    scatter(control_unifrac(box_plot_batch_nums==2 & control_loc_holder'==ii,1),control_unifrac(box_plot_batch_nums==2 & control_loc_holder'==ii,2),35,all_colors_in_order(ii,:),'v','filled');hold on
    scatter(control_unifrac(box_plot_batch_nums==3 & control_loc_holder'==ii,1),control_unifrac(box_plot_batch_nums==3 & control_loc_holder'==ii,2),35,all_colors_in_order(ii,:),'d','filled');hold on

end
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); ylim([-.8,.6]);xlim([-.8,.6])
box on; set(gca,'xtick',[-.8:.7:.6]);set(gca,'ytick',[-.8:.7:.6]); 

figure('Renderer', 'painters', 'Position', [10 10 500 300])
for ii=1:9
    scatter(surgical_unifrac(box_plot_batch_nums==1 & surgical_loc_holder'==ii,1),surgical_unifrac(box_plot_batch_nums==1 & surgical_loc_holder'==ii,2),35,all_colors_in_order(ii,:),'o','filled');hold on
    scatter(surgical_unifrac(box_plot_batch_nums==2 & surgical_loc_holder'==ii,1),surgical_unifrac(box_plot_batch_nums==2 & surgical_loc_holder'==ii,2),35,all_colors_in_order(ii,:),'v','filled');hold on
    scatter(surgical_unifrac(box_plot_batch_nums==3 & surgical_loc_holder'==ii,1),surgical_unifrac(box_plot_batch_nums==3 & surgical_loc_holder'==ii,2),35,all_colors_in_order(ii,:),'d','filled');hold on
end
set(gca,'fontsize',fontsize_holder); set(gca,'TickDir','out');set(gca,'linewidth',linewidth_holder); ylim([-.8,.6]);xlim([-.8,.6])
box on; set(gca,'xtick',[-.8:.7:.6]);set(gca,'ytick',[-.8:.7:.6]); 
