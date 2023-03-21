%% Load abundance data
clear all
clc

chem = load('compound_abundance.mat');

chem_C = table2cell(chem.MetaboanalystParticipantIDPlottedOnly);

compounds = chem.MetaboanalystParticipantIDPlottedOnly.Properties.VariableNames;

%% sort data

chem_ALL = {};

count = 1;
humans_all = cell2mat(chem_C(:,2));
humans = unique(humans_all);


for i=3:(length(compounds));
    
    curr_conc = cell2mat(chem_C(:,i));
    
    
    for z=1:length(humans);
        
        curr_hum_ID = humans_all == z;
        
        curr_hum_conc = curr_conc(curr_hum_ID);
        
        chem_ALL{1,count} = curr_hum_conc;
        chem_ALL{2,count} = humans(z,1);
        chem_ALL{3,count} = mean(curr_hum_conc);
        chem_ALL{4,count} = ste(curr_hum_conc,0,1);
        chem_ALL{5,count} = compounds(1,i);
        count=count+1;
    end

end

%% Plot and save figures


for i = 3:(length(compounds));
    
    figure()
   
    curr_comp = cell2mat(cellfun(@(x) strcmp(compounds(1,i),x), chem_ALL(5,:), 'UniformOutput',false));
    
    curr_abundance = cell2mat(chem_ALL(1,curr_comp));
    
    boxplot(curr_abundance,'Color','k','Symbol','+','OutlierSize',5)
   
    h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
    get(gca,'fontname');  % shows you what you are using.
    set(gca,'fontname','Arial');  % Set it to arial
    
      
    ylabel('Abundance')
    xlabel('human')

    xlim([0.5, 6.5]);

    box off
    
    set(gca, 'XTick',1:length(humans));
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'fontsize',8.2526)
    set(gca,'TickDir','out');
    set(gca,'linewidth',0.5)
    
    x0=0;
    y0=0;
    width=4.3;
    height=4.5;
    set(gcf,'units','centimeters','position',[x0,y0,width,height])
    
    
    title(compounds(1,i));
    
    comp = [compounds{1,i} '.eps'];
    
    saveas(gcf, comp);
    close 
end   
%% Run a Fisher's exact permutation test

chem_stats={};

groups_per=string(unique([chem_ALL{2,:}]));

count = 1;

for z = 3:(length(compounds));
    
    curr_comp = cell2mat(cellfun(@(x) strcmp(compounds(1,z),x), chem_ALL(5,:), 'UniformOutput',false));
    
    curr_abundance =chem_ALL(1,curr_comp);
    
    curr_ab = cell2mat(chem_ALL(1,curr_comp));
    
    
    [pVals, pVals_C, H] = EST_multipleTests(curr_abundance,'title',compounds{1,z},...
    'repNb',20000,'correctionType','benjami_P','fileOutput','./mean_percentage_stats.txt',...
    'alpha',0.05,'groupNames',groups_per,...
    'testFunc','EST_permTest_meanDiff');


    chem_stats{1,count} = c(:,1);
    chem_stats{1,count} (:,2) = c(:,2);
    chem_stats{1,count} (:,3)= pVals; % friedman multcomp values
    chem_stats{1,count}(:,4) = pVals_C;
    chem_stats{1,count}(:,5) = H;
    chem_stats{2,count} = compounds{1,z};
    
    count = count+1;
    close all
end