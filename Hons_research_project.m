%% CRISPR Cancer Dependency Project

%{
Mapping Breast Cancer Cell Signalling Pathways with CRISPR-Cas9 Gene 
Knockout Data

Breast cancer is a heterogeneous disease, and the mechanisms that drive
its progression are still not fully understood. Cell signalling pathways
play a crucial role in the development and progression of breast cancer,
and targeting them has emerged as a promising therapeutic approach. 
However, the specific signalling pathways that are critical for breast 
cancer development and progression are not fully elucidated. CRISPR-Cas9 
gene knockout technology has emerged as a powerful tool for studying gene 
function and identifying key regulators of cellular processes. CRISPR is a 
genome editing technology that utilises RNA-guided enzymes to create 
double-strand breaks in DNA, enabling the precise and targeted knockout of 
specific genes. In this project, we will use CRISPR-Cas9 gene knockout 
data to identify the dependencies of breast cancer on various cell 
signalling pathways. Furthermore, by analysing the patterns of gene 
essentiality across breast cancer cell lines, we aim to infer the relative 
importance of different signalling pathways in breast cancer development 
and progression. The findings of this study could have significant 
implications for the development of targeted therapies for breast cancer 
by identifying critical signalling pathways that can be targeted for 
therapeutic intervention.
%}

% have a new slate 
close all force; clc; clear 

% change the directory -- choose the directory on your computer
cd /Users/jclan/OneDrive/Desktop/Code_and_files

%% Loading the data and cleaning the data

% load the Reactome pathway database
reactome = readtable('Reactome_2016.xlsx','ReadVariableNames',false) ;

% clean up the data, remove the second column because it is of no use and
% remove information of human pathway that are not necessarily pathways
reactome.Properties.VariableNames(1) = "Pathway" ;
reactome.Var2 = [] ;

reactome( ~contains( reactome.Pathway,'Homo sapiens', ...
    'IgnoreCase', true), :) = [];
reactome.Pathway = extractBefore( reactome.Pathway ,' Homo sapiens') ;

% load the CRISPR data 
% the data is found at the following directory 
% /Users/sinkala/Documents/MATLAB/Depmap_data
crispr = readtable('CRISPR_(DepMap_22Q2_Public+Score,_Chronos).csv') ;
crispr.Properties.VariableNames(1) = "DepMap_ID" ;

% get the sample information from the achille project so i can return only
% the cell line of breast cancer 
% Info found here '/Users/sinkala/Documents/MATLAB/DarkGenomeProject/
sample_info = readtable('Achilles_sample_info_20Q4.csv') ;
sample_info = sample_info( ismember( sample_info.lineage, 'breast'), :) ;

% get only the breast cancer cell lines in the crispr data and include the
% cell line names. then remove the depmap ID
crispr = innerjoin( sample_info(: , ...
    {'DepMap_ID','stripped_cell_line_name'}) , crispr ) ;
crispr.DepMap_ID = [] ;

% change the first variable names 
crispr.Properties.VariableNames(1) = "cell_line" ;


%% Perform Welch's t-test for Pathway Enrichment Analysis

% turn of off the warnings 
warning('off')

% here is the folder for saving the figures 
if ~exist('effect_plots','dir')
    mkdir effect_plots
end

% here is the data for the welch test 
crispr_welch = crispr(:,2:end) ;

% Initialise the results table and the mean values 
dep_pathways_welch = reactome(:,1);
dep_pathways_welch.mean_inPathway = zeros(size(reactome, 1), 1);
dep_pathways_welch.mean_notInPathway = zeros(size(reactome, 1), 1);

% Add new columns to dep_pathways_welch for Glass Delta and its confidence
% intervals
dep_pathways_welch.Glass_delta = zeros(size(reactome, 1), 1);
dep_pathways_welch.CI_lower = zeros(size(reactome, 1), 1);
dep_pathways_welch.CI_upper = zeros(size(reactome, 1), 1);

% Initialize a table to store t-values, p-values,
dep_pathways_welch.tValue = zeros(size(reactome, 1), 1);
dep_pathways_welch.pValue = zeros(size(reactome, 1), 1);

% here are the colours for the boxplots
theColours = [ 0.9019   0.2941  0.20784 ; 0.3019   0.7333   0.835 ] ;

% Iterate over each pathway in Reactome
for ii = 1:size(reactome, 1)
    
    % print something to screen
    if rem(ii,20) == 0
        fprintf('\nRunning Welch t-test for pathway %s number %d of %d\n', ...
            reactome.Pathway{ii}, ii, height(reactome))
    end
    
    % Extract the entire row for the current pathway as a table
    pathway_row = reactome(ii, 2:end);
    
    % Convert the table to a cell array
    pathway_cells = table2cell(pathway_row);
    
    % Remove any empty cells or non-string values
    genes_in_pathway = pathway_cells(cellfun(@(x) ischar(x) || ...
        isstring(x), pathway_cells));
    
    % Remove empty cells (if any)
    genes_in_pathway(cellfun('isempty',genes_in_pathway)) = [] ;
    
    % Return only genes that are present in the CRISPR in the
    % genes_in_pathways
    genes_in_pathway = genes_in_pathway( ismember( genes_in_pathway , ...
        crispr_welch.Properties.VariableNames) ) ;
    
    % Extract CRISPR scores for genes in the pathway and not in the pathway
    scores_in_pathway = crispr_welch{:, genes_in_pathway};
    scores_not_in_pathway = crispr_welch{:, ~ismember( ...
        crispr_welch.Properties.VariableNames, genes_in_pathway) };
    
    % convert to arrays
    scores_in_pathway = scores_in_pathway(:) ;
    scores_not_in_pathway = scores_not_in_pathway(:) ;
    
    % Perform Welch's t-test
    [h,p,~,stats] = ttest2( scores_in_pathway, ...
        scores_not_in_pathway, 'Vartype','unequal');
    
    % Store the computed t-value and p-value
    dep_pathways_welch.tValue(ii) = stats.tstat;
    dep_pathways_welch.pValue(ii) = p;
        
    % add the mean values to the table
    dep_pathways_welch.mean_inPathway(ii) = mean(scores_in_pathway);
    dep_pathways_welch.mean_notInPathway(ii) = mean(scores_not_in_pathway);
   
    % *****************  Glass's delta compututations *******************
    
    % Glass's delta
    % It measures the size of the difference between two groups.
    % A rule of thumb for interpreting the magnitude:
    % Small effect: delta = 0.2
    % Medium effect: delta = 0.5
    % Large effect: delta = 0.8

    % Calculate the mean of scores for genes in the pathway (treatment
    % group)
    mean1 = mean(scores_in_pathway);
    
    % Calculate the mean of scores for genes not in the pathway (control
    % group)
    mean2 = mean(scores_not_in_pathway);
    
    % Calculate the standard deviation of scores for genes not in the
    % pathway (control group)
    std_control = std(scores_not_in_pathway);
    
    % Calculate Glass's delta using the formula
    delta = (mean1 - mean2) / std_control;
    
    % Calculate the number of scores for genes in the pathway
    n1 = length(scores_in_pathway);
    
    % Calculate the number of scores for genes not in the pathway
    n2 = length(scores_not_in_pathway);
    
    % Calculate the variance of scores for genes in the pathway
    variance1 = var(scores_in_pathway) ;
    
    % Calculate the variance of scores for genes not in the pathway
    variance2 = var(scores_not_in_pathway);

    % Confidence interval for Glass's delta

    % Calculate the degrees of freedom for Welch's t-test
    df = ((variance1 / n1 + variance2 / n2)^2) / ((variance1 / n1)^2 / ...
        (n1 - 1) + (variance2 / n2)^2 / (n2 - 1));
    
    % Calculate the critical t-value for a 95% CI
    t_crit = tinv(0.975, df); 
    
    % Calculate the margin of the confidence interval using the control
    % group standard deviation
    CI_margin = t_crit * sqrt(1/n1 + 1/n2) * std_control;
    
    % Determine the lower and upper bounds of the 95% confidence interval
    % for Glass's delta
    CI_lower = delta - CI_margin;
    CI_upper = delta + CI_margin;

    % Store Glass delta and its CI in the results table
    dep_pathways_welch.Glass_delta(ii) = delta;
    dep_pathways_welch.CI_lower(ii) = CI_lower;
    dep_pathways_welch.CI_upper(ii) = CI_upper;
    
    % *******************************************************************
    
    
    % *********** plot the boxplots and save them to a folder ***********

    % Match any character that is not a letter or a number and the
    % replacement using regular expression to see if the file already
    % exists on the computer 
    pattern = '[^a-zA-Z0-9]';
    figure_path = fullfile('effect_plots', ...
        [ regexprep(reactome.Pathway{ii},  pattern, '_'), '.png']) ;

    % only for signicant results and only when the plot does not exist
    if p < 1e-10 && ~exist(figure_path,'file')
        
        % get the scores and the groups
        all_scores = [scores_in_pathway ; scores_not_in_pathway ] ;
        theGroups =  ...
            [ repmat({'In Pathway'},length(scores_in_pathway), 1) ; ...
            repmat({'Not In Pathway'},length(scores_not_in_pathway), 1) ] ;
        theGroups = categorical( theGroups) ;
        
        % here is the figure
        t = tiledlayout(1,3,'TileSpacing','compact');
        
        % plot the box plots for the for the crispr scores
        nexttile([1, 2]);
        colourBoxPlot(all_scores,  theGroups , theColours, true)
        hold on
        
        % set the line width of the box plots
        set(findobj(gca,'type','line'),'linew',1)
        
        % anotate the graphs
        if p == 0
            text(1.2, max(all_scores), ...
                ['p < ',convertPValue2SuperScript(p)], 'FontSize',13)
        else
            text(1.2, max(all_scores), ...
                ['p = ',convertPValue2SuperScript(p)], 'FontSize',13)
        end
        
        % change some plot properties
        set(gca,'FontSize',12,'FontWeight','normal','LineWidth',1,...
            'Box','off')
        
        % add the y-label and the figure title ot the figure
        ylabel('Gene Effect','FontWeight','normal')
        
        
        % get the limits of the current boxplots 
        curYlim = get(gca,'YLim') ;
        
        % here is the histgram
        ax2 = nexttile(t);
        hold on
        
        % Estimate the kernel density
        [x, f] = ksdensity( scores_not_in_pathway);
        
        % Plot the density with filled color
        fill( ax2, x, f, theColours(1,:) , 'FaceAlpha', 0.7); 
        
        % Estimate the kernel density
        [x, f] = ksdensity(scores_in_pathway);
        
        % Plot the density with filled color
        fill(ax2, x, f, theColours(2,:),'FaceAlpha', 0.7 ); 
        
        % change the plot properties.
        set(gca,'FontSize',12,'FontWeight','normal','LineWidth',1,...
            'Box','off','YLim',curYlim)

        % here is the title for both plots
        title(t,reactome.Pathway{ii} ,'FontSize',14 ,'FontWeight' ,'bold')
   
        hold off
        
        % let get the figure names we will use regular expression to remove
        % characters that window does not want in the names of files 

        % here is the filenema
        file_name =  regexprep(reactome.Pathway{ii} , pattern ,'_') ;

        % save the figure
         try   
             % for Linux systems
             saveas(gcf,['effect_plots/',file_name  ,'.png'],'png');
         catch
             try 
                 % for windows system
                  saveas(gcf,['effect_plots\',file_name  ,'.png'],'png');
             catch ME
                 fprintf('\nCould not save the file to the drive\n') 
             end
         end
        
        % close the figure
        close(gcf);
        
    end
    % ************************************************************
    
end

% Adjust p-values for multiple testing using Benjamini-Hochberg method and
% then sort the table according to the p-value
dep_pathways_welch.adjPvalue = mafdr(dep_pathways_welch.pValue, ...
    'BHFDR', true);
dep_pathways_welch = sortrows(dep_pathways_welch, 'adjPvalue', 'ascend');

% Display the top pathways from the Welch's t-test analysis
head(dep_pathways_welch)

% save the results to excel 
writetable( dep_pathways_welch, 'dep_pathways_welch_results.xlsx') 

% Assuming your dataset is stored in a variable named 'crispr'
% and the column containing cell lines is named 'cellLines'

% Extract the column containing cell lines
cellLines = crispr(:, 'cell_line');

% Use the 'unique' function to get unique cell lines
uniqueCellLines = unique(cellLines);

% Get the total number of unique cell lines
numUniqueCellLines = numel(uniqueCellLines);
% Display the result
disp(['Number of unique cell lines: ', num2str(numUniqueCellLines)]);


%% Plot the graph 

% check that the variable exists
if ~exist('dep_pathways_welch','var')
    dep_pathways_welch = readtable('dep_pathways_welch_results.xlsx')  ;
end

% let see the top pathways on which the cancer cell lines are most depended
% on for their fitness

% we will use the t-values ot plot the figures 
dep_pathways_welch = sortrows(dep_pathways_welch,'tValue','ascend');

% get the required column 
dep_pathways_welch1 = dep_pathways_welch(:,[7,1]) ;

% get the top 20 pathways 
dep_pathways_welch11 = dep_pathways_welch1(1:20,1:2);

% convert to categorical arrays for easy plotting of bar graphs 
dep_pathways_welch11.Pathway = categorical( ...
    dep_pathways_welch11.Pathway, ...
    flipud(dep_pathways_welch11.Pathway ) );

% here is the horizontal bargraph 
barh(dep_pathways_welch11.Pathway,dep_pathways_welch11.tValue, ...
   'DisplayName','dep_pathways_welch11.tValue')
set(gca, 'FontSize', 12); % Change the font size to your desired value


% Save as a PNG file for write up
saveas(gcf, 'dep_pathways_welch_plot.png');
%% Get genes involved in transcription 

transcription = dep_pathways_welch( contains( ...
    dep_pathways_welch.Pathway, 'transcription', 'IgnoreCase', true), :) ;


%% Load GDSC data 

% Load the FDR data with explicit variable names
fdr = readtable('FDR2.csv','ReadVariableNames',true);

% Rename the 'PUTATIVE_TARGET' column to 'TARGET'
fdr.Properties.VariableNames{'PUTATIVE_TARGET'} = 'TARGETS';

% Filter for rows where TCGA_DESC is 'BRCA'
fdr_brca = fdr(strcmp(fdr.TCGA_DESC, 'BRCA'), :);

% Rename pathway_name to target_pathway
fdr_brca.Properties.VariableNames{'PATHWAY_NAME'} = 'TARGET_PATHWAY';

% Define the columns to remove
columns_to_remove = {'DATASET', 'NLME_RESULT_ID', 'NLME_CURVE_ID', 'COSMIC_ID', ...
    'SANGER_MODEL_ID', 'TCGA_DESC', 'DRUG_ID', 'COMPANY_ID', 'WEBRELEASE', ...
    'MIN_CONC', 'MAX_CONC', 'AUC', 'RMSE'};

% Remove the specified columns
fdr_brca(:, columns_to_remove) = [];

% Load the drugs2 dataset with explicit variable names
drugs = readtable('Drugs2.csv','Delimiter',',','ReadVariableNames',true);

% Define the columns to remove
columns_to_remove = {'drug_id', 'synonyms', 'pubchem'};

% Remove the specified columns
drugs(:, columns_to_remove) = [];

% Convert all column names to uppercase
drugs.Properties.VariableNames = upper(drugs.Properties.VariableNames);

% Rename pathway_name to target_pathway
drugs.Properties.VariableNames{'PATHWAY_NAME'} = 'TARGET_PATHWAY';

% Merge drugs and FDR_BRCA
% Convert 'DRUG_NAME' column in both tables to uppercase for consistency
drugs.DRUG_NAME = upper(drugs.DRUG_NAME);
fdr_brca.DRUG_NAME = upper(fdr_brca.DRUG_NAME);

% Merge the tables based on DRUG_NAME and TARGET
gdsc = innerjoin(fdr_brca, drugs, 'Keys', {'DRUG_NAME', 'TARGETS', 'TARGET_PATHWAY'});

%% Data Preprocessing

% Assuming gdsc.Z_SCORE is a cell array with potentially mixed data types
% Convert commas to periods in Z_SCORE
gdsc.Z_SCORE = cellfun(@(x) strrep(x, ',', '.'), gdsc.Z_SCORE, 'UniformOutput', false);

% Convert Z_SCORE to a numeric array, replacing non-numeric values with NaN
gdsc.Z_SCORE = cellfun(@str2double, gdsc.Z_SCORE, 'UniformOutput', false);
gdsc.Z_SCORE = cell2mat(gdsc.Z_SCORE);

% Filter out NaN values from Z_SCORE
gdsc = gdsc(~isnan(gdsc.Z_SCORE), :);

% Turn off warnings
warning('off');

% Create a directory for saving the plots
if ~exist('ic50_effect_plots/', 'dir')
    mkdir ic50_effect_plots;
end

%% Perform Analysis and Generate Plots

ic50_results = table();
unique_pathways = unique(gdsc.TARGET_PATHWAY);
theColours = [0.3019 0.7333 0.835; 0.9019 0.2941 0.20784];


for ii = 1:length(unique_pathways)
    current_pathway = unique_pathways{ii};
    
    ic50_targeting_pathway = gdsc.Z_SCORE(strcmp(gdsc.TARGET_PATHWAY, ...
        current_pathway));
    ic50_not_targeting_pathway = gdsc.Z_SCORE(~strcmp(gdsc.TARGET_PATHWAY, ...
        current_pathway));
    
    [~, p, ~, stats] = ttest2(ic50_targeting_pathway, ...
        ic50_not_targeting_pathway, 'Vartype', 'unequal');
    
    ic50_results = [ic50_results; table({current_pathway},...
        stats.tstat, p, 'VariableNames',...
        {'Pathway', 'T_Value', 'P_Value'})];
    
    figure();
    t = tiledlayout(1, 2, 'TileSpacing', 'compact');
    
    nexttile([1, 1]);
    all_scores = [ic50_targeting_pathway; ic50_not_targeting_pathway];
    theGroups = [repmat({'Not Targeting'}, ...
        length(ic50_not_targeting_pathway), 1); ...
        repmat({'Targeting'}, length(ic50_targeting_pathway), 1)];
    theGroups = categorical(theGroups);
    colourBoxPlot(all_scores, theGroups, theColours, true);
    hold on;
    set(gca, 'FontSize', 12, 'LineWidth', 1, 'Box', 'off');
    
    p_str = ['p = ', sprintf('%.4g', p)];
    p_str = strrep(p_str, 'e', ' x 10^');
    text(mean(xlim), max(all_scores) * 0.9 + (max(all_scores) - min(all_scores)) * 0.05, p_str, 'FontSize', 10, ...
        'HorizontalAlignment', 'center');

    
    title([current_pathway], 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    nexttile([1, 1]);
    hold on;
    [f, y] = ksdensity(ic50_targeting_pathway);
    fill(f, y, theColours(2,:), 'FaceAlpha', 0.7);
    [f, y] = ksdensity(ic50_not_targeting_pathway);
    fill(f, y, theColours(1,:), 'FaceAlpha', 0.7);
    hold off;
    
    ax = gca;
    ax.XTickLabelRotation = 0;
    
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) pos(3)*1.5 pos(4)]);
    
    saveas(gcf, ['ic50_effect_plots/',...
        regexprep(current_pathway, '[^a-zA-Z0-9]', '_'), '.png']);
    close(gcf);
end

ic50_results = sortrows(ic50_results,'P_Value','ascend');
disp(ic50_results);
writetable(ic50_results, 'ic50_pathway_comparison_results.xlsx');

% change notation of pvalues
ic50_results.P_Value = arrayfun(@(p) strrep(sprintf('%.4g', p), 'e', 'x10^'), ic50_results.P_Value, 'UniformOutput', false);

% Assuming ic50_results is a table and gdsc is the dataset with drugs and targets
for i = 1:height(ic50_results)
    pathway = ic50_results.Pathway{i};
    
    % Find drugs targeting the pathway
    drugs_targeting_pathway = unique(gdsc.DRUG_NAME(strcmp(gdsc.TARGET_PATHWAY, pathway)));
    ic50_results.Drugs_Targeting{i} = strjoin(drugs_targeting_pathway, ', ');
    
    % Find targets associated with the pathway
    targets_in_pathway = unique(gdsc.TARGETS(strcmp(gdsc.TARGET_PATHWAY, pathway)));
    ic50_results.Targets_In_Pathway{i} = strjoin(targets_in_pathway, ', ');
end


% Assuming ic50_results is a table and you have added Drugs_Targeting and Targets_In_Pathway columns

% Create a new table including the original ic50_results and the new columns
new_ic50_results = ic50_results;

% Rename the new columns to remove the prefix
new_ic50_results.Properties.VariableNames{'Drugs_Targeting'} = 'New_Drugs_Targeting';
new_ic50_results.Properties.VariableNames{'Targets_In_Pathway'} = 'New_Targets_In_Pathway';

% Save the new table to an Excel file
writetable(new_ic50_results, 'ic50_results_with_info.xlsx');

%% similar pathways 

gdsc_pathways = {
    'ABL signaling'
    'Apoptosis regulation'
    'Cell cycle'
    'Chromatin histone acetylation'
    'Chromatin histone methylation'
    'Chromatin other'
    'Cytoskeleton'
    'DNA replication'
    'EGFR signaling'
    'ERK MAPK signaling'
    'Genome integrity'
    'Hormone-related'
    'IGF1R signaling'
    'JNK and p38 signaling'
    'Metabolism'
    'Mitosis'
    'Other'
    'Other, kinases'
    'PI3K/MTOR signaling'
    'Protein stability and degradation'
    'RTK signaling'
    'WNT signaling'
    'p53 pathway'
};
 
 
 
% Define a function to find pathways with similar words
find_similar_pathways = @(query, pathways) pathways(contains(lower(pathways), lower(query)));
 
% Iterate over each pathway in gdsc_pathways and find similar pathways in Reactome
common_pathways = cell(0, 2);
for i = 1:numel(gdsc_pathways)
    query_pathway = gdsc_pathways{i};
    similar_pathways = find_similar_pathways(query_pathway, dep_pathways_welch1.Pathway);
    common_pathways = [common_pathways; repmat({query_pathway}, numel(similar_pathways), 1), similar_pathways];
end



% Assuming ic50_results is a table with Pathway and P_Value columns


% Initialize a cell array to store the p-values
p_values = cell(size(common_pathways));

% Loop through each pathway in common_pathways
for i = 1:numel(common_pathways)
    % Find the corresponding pathway in ic50_results
    idx = strcmp(ic50_results.Pathway, common_pathways{i});
    
    % Check if the pathway is found
    if any(idx)
        % Extract the corresponding p-value
        p_values{i} = ic50_results.P_Value(idx);
    else
        % If pathway is not found, assign NaN or any other value as needed
        p_values{i} = NaN;
    end
end

% Combine common_pathways and p_values into a table
common_pathways_table = table(common_pathways', p_values', 'VariableNames', {'Pathway', 'P_Value'});

% Display the updated table
disp(common_pathways_table);


% Assuming dep_pathways_welch1 is a table with columns 'Pathway' and 'T_Value'
% Initialize a cell array to store the t-values
t_values = cell(size(common_pathways, 1), 1);

% Loop through each pathway in common_pathways
for i = 1:size(common_pathways, 1)
    % Find the corresponding pathway in dep_pathways_welch1
    idx = strcmp(dep_pathways_welch1.Pathway, common_pathways{i, 2});
    
    % Check if the pathway is found
    if any(idx)
        % Extract the corresponding t-value
        t_values{i} = dep_pathways_welch1.tValue(idx);
    else
        % If pathway is not found, assign NaN or any other value as needed
        t_values{i} = NaN;
    end
end

% Add the t-values to column 3 of common_pathways
common_pathways(:, 3) = t_values;

% Display the updated common_pathways
disp(common_pathways);

% Assuming ic50_results is a table with columns 'Pathway' and 'P_Value'
% Initialize a cell array to store the P_Values
p_values = cell(size(common_pathways, 1), 1);

% Loop through each pathway in common_pathways
for i = 1:size(common_pathways, 1)
    % Find the corresponding pathway in ic50_results
    idx = strcmp(ic50_results.Pathway, common_pathways{i, 1});
    
    % Check if the pathway is found
    if any(idx)
        % Extract the corresponding P_Value
        p_values{i} = ic50_results.P_Value(idx);
    else
        % If pathway is not found, assign NaN or any other value as needed
        p_values{i} = NaN;
    end
end

% Insert the P_Values as a new column next to the first column in common_pathways
common_pathways = [common_pathways(:,1) p_values common_pathways(:,2:end)];

% Define the column headings
headings = {'GDSC_pathways', 'P_Value', 'Reactome_pathways', 'T_Value'};

% Create a cell array to hold the combined data
combined_data = cell(size(common_pathways, 1) + 1, size(common_pathways, 2));

% Insert the original data into the combined array
combined_data(2:end, :) = common_pathways;

% Insert the headings in the first row
combined_data(1, :) = headings;

% Display the updated data
disp(combined_data);

writecell(combined_data, 'combined.csv')


%% *********************** Internal Functions ************************

% ======================= another function =========================
function pSuperScript = convertPValue2SuperScript(p)

    % converts the number to scientific superscript for printing on a
    % figure
    pS = num2str(p) ;

    % get the first number
    firstNumbers = extractBefore(pS,'e') ;

    % check if there is a decimal place. then only get the first 4 numbers
    if contains( firstNumbers  ,'.')
        try
            firstNumbers = firstNumbers(1:4) ;
        catch
            firstNumbers = firstNumbers(1:3) ;
        end
    end
    
    % get the correctly formated p value
    pSuperScript = sprintf('%s x 10^{%d}', firstNumbers, ...
        str2double(extractAfter(pS, 'e') )) ;
    
    % if the p value is large
    if p > 0.0001
       pSuperScript = sprintf('%0.4f', p) ;
    elseif p == 0
         pSuperScript = sprintf('1 x 10^{%d}', -300) ;
    end

end

% *********************** end of function ************************
% ****************************************************************

% ======================= another function =========================

function colourBoxPlot(plotData,groups, color, includeScatter)

% set the color to the box plots
if nargin == 2 || isempty(color)
    rng(6);
    color = rand(length(unique(groups )),3) ;
    if ~exist('includeScatter','var')
        includeScatter = false;
    end
end

% plot the data
boxplot(plotData,groups,'Color', flipud(color),'Symbol','k+', ...
    'OutlierSize',2) ;

% set some figure properties and add title ot the figure
set(gca,'FontSize',14,'LineWidth',1.5,'Box','off')

% set the line width of the box plots
set(findobj(gca,'type','line'),'linew',2)
set(findobj(gca,'Tag','Lower Whisker'),'LineStyle','-')
set(findobj(gca,'Tag','Upper Whisker'),'LineStyle','-')

% set the color of the box plots
h4 = findobj(gca,'Tag','Box') ;
for kk=1:length(h4)
    patch(get(h4(kk),'XData'),get(h4(kk),'YData'),...
        color(kk,:),'FaceAlpha', 0.3,'LineStyle','-');
end

% add a scatter plot if we that is true
if includeScatter
   
    % get the unique groups 
    uniqueGroups = unique(groups) ;

    % add the scatter plots
    hold on
    
    % add scatter plot to the box plots
    groupSc1 = plotData(groups == uniqueGroups(1,1));
    groupSc2 = plotData(groups == uniqueGroups(2,1));
    
    % set the size the size of the marker
    makerSz = 15;
     
    % get only a smaller subset if the data points are too many
    if length(groupSc1) > 600
        groupSc1 = randsample(groupSc1,600) ;
        
        % change the markerSize 
        makerSz = 15 ;
    end
    
    if length(groupSc2) > 600
        groupSc2 = randsample(groupSc2,600) ;
    end
    
    x = ones(length(groupSc1)).*(1+(rand(length(groupSc1))-0.5)/5) ;
    x1 = ones(length(groupSc2)).*(1+(rand(length(groupSc2))-0.5)/10);
    
    % here is the first scatter plot
    scatter(x(:,1), groupSc1, makerSz, color(2,:),'filled', ...
        'MarkerFaceColor',color(2,:),'Marker','o','MarkerFaceAlpha',0.8)
    hold on
    scatter(x1(:,2).*2, groupSc2, makerSz, color(1,:),'filled', ...
        'MarkerFaceColor',color(1,:),'Marker','o','MarkerFaceAlpha',0.8)
    
    hold off
end

end


