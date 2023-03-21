%% load concatenated data from ivTrace analysis
clc
clear all

data=load('data_field_2020.mat'); % load the concatenated cell containing the semifield data analyzed with ivTrace

data=data.data_AllC;

%% Export all the data for permutation analysis on Python

% this section exports the semi-field data collected in 2020 in a csv file
% that can be imported into python for statistical analyses using a
% permutation apporach

data_new = data(1,:); % get landings
data_new(2,:) = data(5,:); % Get experiment type
data_new(3,:) = data(6,:); % Get tratment type
data_new(4,:) = data(7,:); % Get position in cage
data_new(5,:) = data(8,:); % get data


data_table = cell2table(data_new);
writetable(data_table,'data_permutation.csv') % export a table for input in Python


%% select data type

exp_type=data(5,:);

experiment_list = {"CO2",  "1-human", "2-human"};


selection = listdlg("ListString",experiment_list,"SelectionMode", "single"); % select which experiment type to analyze 
experiment_type = experiment_list{1,selection};


type_ID = cell2mat(cellfun(@(x) strcmp(x,experiment_type), exp_type ,'UniformOutput',false));

data_type = data(:,type_ID); % make a cell containing the subset of data to be analyzed (i.e. the selected expriment type)

type = unique([data_type{6,:}]);




%% calculate mean landings per group

groupNumbers = nan(1,size(data_type,2));

exp_pos = data_type(6,:); % get the treatment information

landings = cell2mat(data_type(1,:)); % get the total landing information

real_pos = data_type(7,:); % get the position relative to the cage information

mean_landings=nan(1,length(type)); % make an empty vector to store the mean landings



for i=1:length(type);
    
    exp_ID=cell2mat(cellfun(@(x) strcmp(x,type(1,i)), exp_pos,'UniformOutput',false)); % find the data from the current treatment
    groupNumbers(exp_ID)=i;
    
    mean_landings(1,i)=mean(landings(exp_ID)); % calculate the mean landings for the treatment
    
end

landings_cell={}; % make a cell to store the total landings sorted by treatment


for i=1:length(unique(groupNumbers));
    
    curr_type_ID=groupNumbers==i;
    landings_cell{1,i}=landings(curr_type_ID)'; % store the total landings for each experiment
end


ste_groups=cellfun(@(x) ste(x,1,1), landings_cell,'UniformOutput',false); % calculate the standard error of the mean


exp_date = data_type(8,:);
dates = unique(exp_date); % get the night dates


night={};

for z=1:length(dates)
    night{1,z}=['night ' num2str(z)];
end



%% Plot total landings


for i=1:length(dates)
    
    date_ID = cell2mat(cellfun(@(x) strcmp(x,dates(1,i)), exp_date,'UniformOutput',false));
    curr_land = landings(date_ID);
    curr_real_pos = string(real_pos(date_ID));
        
    [max_data_pos{1,i},max_ID] = max(curr_land);    
    max_data_pos{2,i} = curr_real_pos(max_ID);
    
    total_land{1,i} = sum(curr_land);
    curr_gn = groupNumbers(date_ID);
    [gn_sorted, gn_order] = sort(curr_gn);
    curr_land = curr_land(gn_order);
    plot(gn_sorted,curr_land,"o-");
    
    hold on
end

hold on 


errorbar((1:length(type)),mean_landings,cell2mat(ste_groups),'s')

xticklabels(type)
xlim([0.5 length(type)+0.5])
ylabel('Number of landings')
xlabel('Position')

legend(night)
box off

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16)


hold off



%% Calculate landing percentages

exp_date=data_type(8,:);
dates=unique(exp_date);
landing_percentages=nan(1,length(exp_date));
count=1;

for i=1:length(dates)
    
    date_ID=cell2mat(cellfun(@(x) strcmp(x,dates(1,i)), exp_date,'UniformOutput',false));
   
    landings_date=landings(date_ID);
    total_date_landings=sum(landings_date);
    landing_percentages(1,count:count+7)=(landings_date/total_date_landings)*100;
    count=count+8;
end

per_ID=~isnan(landing_percentages);
landing_percentages=landing_percentages(per_ID);
exp_pos=data_type(6,:);
exp_pos=exp_pos(per_ID);
    
mean_percentage=nan(1,length(type));
median_percentage=nan(1,length(type));
ste_percentage = nan(1,length(type));

for i=1:length(type);
    
    exp_ID=cell2mat(cellfun(@(x) strcmp(x,type(1,i)), exp_pos,'UniformOutput',false));
    groupNumbers(exp_ID)=i;
    
    mean_percentage(1,i)=mean(landing_percentages(exp_ID));
    median_percentage(1,i)=median(landing_percentages(exp_ID));    
    ste_percentage(1,i) = ste(landing_percentages(exp_ID),1,2);
end

groupNumbers=groupNumbers(per_ID);


%% plot landing percentage
close 

figure()

for i=1:length(dates)
    
    date_ID=cell2mat(cellfun(@(x) strcmp(x,dates(1,i)), exp_date,'UniformOutput',false));
    curr_per=landing_percentages(date_ID);
    curr_gn=groupNumbers(date_ID);
    [gn_sorted, gn_order]=sort(curr_gn);
    curr_per = curr_per(gn_order);
    plot(gn_sorted,curr_per,"o-");
    
    hold on
end

%plot((1:length(type)),mean_percentage,'k*', 'MarkerSize', 8, 'LineWidth', 2)

errorbar((1:length(type)),mean_percentage,ste_percentage,'s')


xticklabels(type)
xlim([0.5 length(type)+0.5])
ylabel('% of landings')
ylim([0 105])
xlabel('Position')
title("landing percentage with mean")
legend(night)

box off

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16)

%% calculate preference index for 2-human data

human_landings= landings_cell(1,8:9);

PI = (human_landings{1,1}-human_landings{1,2})./(human_landings{1,1}+human_landings{1,2});
PI_mean = mean(PI);
PI_sem = ste(PI,0,1);

for i=1:length(PI);
    
    
    scatter(0,PI(i,1),100)
    hold on
end

limits_x=[-0.1 0.1];
set(0,'DefaultLegendAutoUpdate','off')
legend(night)

errorbar(0, PI_mean, PI_sem,'CapSize',18)
ylim([-1 1])
xlim(limits_x)
ylabel("Preference index")
plot(limits_x,[0 0],"k")

box off
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16)


%% distance between humans and PI

close all

land_diff=landings_cell{1,9}-landings_cell{1,8};

[pos pos_id]=sort([7 1 3 6 2 4 5]);

PI_pos = (PI(pos_id))';

dist = [5.32 9.83 12.84 13.9 12.84 9.83 5.32];

land_diff_s=land_diff(pos_id);


p=polyfit(dist,PI_pos,1);

yfit = polyval(p,dist);
yresid = PI_pos - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(PI_pos)-1) * var(PI_pos);
rsq = 1 - SSresid/SStotal;

figure()


for i=1:length(PI);
    
    
    scatter(dist, PI_pos')
    hold on
end
scatter(dist, PI_pos,'o')

legend(night)


ylabel("Preference index")
xlabel("Distance (m)")
% a = get(gca,'XTickLabel');  
 set(gca,'fontsize',16)
 xlim([0 max(dist)+1])
 ylim([-1 1])
 
h=lsline;
set(h,"LineWidth",3)
txt=['R2 = ' num2str(rsq)];

t=text(min(dist)+0.2,max(PI_pos),txt);
t.FontSize=16;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16)
 
 box off

%% 2020 weather data analysis


weather_data=load("weather_2020.mat");
weather_data=weather_data.hobo_per_night; % load all hobo data

exp_date = data_type(8,:);

dates=unique(exp_date);

weather_dates=weather_data(6,:);

weather={};

for i=1:length(dates);
    
    weather_dates_ID=cell2mat(cellfun(@(x) strcmp(x,dates(1,i)),weather_dates,'UniformOutput',false));
    weather_curr=weather_data(:,weather_dates_ID);
    weather(:,i)=weather_curr;
end


%% Temperature

temperature = weather(4,:);

exp_temp = cellfun(@(x) x(13:end), temperature, 'UniformOutput',false);

mean_temp = cellfun(@(x) mean(x), exp_temp, 'UniformOutput',false);

landing_total = nan(1,length(dates));
count = 1;

for i=1:length(dates)
    
    date_ID=cell2mat(cellfun(@(x) strcmp(x,dates(1,i)), exp_date,'UniformOutput',false));
   
    landings_date = landings(date_ID);
    landing_total(1,i) = sum(landings_date);
    count = count+8;
end

mean_t = cell2mat(mean_temp);

%% calculate R2 temp

p = polyfit(mean_t,landing_total,1);

yfit = polyval(p,mean_t);
yresid = landing_total - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(landing_total)-1) * var(landing_total);
rsq = 1 - SSresid/SStotal;

%% plot mean temeprature per night vs landings 
close 

figure()
scatter(mean_t, landing_total)
ylabel("Total landings")
xlabel("Mean temperature (°C)")
title('Temeprature VS Landings')
h=lsline;
set(h,"LineWidth",3)
txt=['R2 = ' num2str(rsq)];

xlim([16 22])


t=text(min(mean_t)+0.2,max(landing_total)-100,txt);
t.FontSize=16;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16)
set(gca,'TickDir','out'); % The only other option is 'in'
%% get RH data

RH = weather(5,:);

exp_RH = cellfun(@(x) x(13:end), RH, 'UniformOutput',false);

mean_RH = cellfun(@(x) mean(x), exp_RH, 'UniformOutput',false);
mean_RH = cell2mat(mean_RH);

%% calculate R2 RH

p=polyfit(mean_RH,landing_total,1);

yfit = polyval(p,mean_RH);
yresid = landing_total - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(landing_total)-1) * var(landing_total);
rsq = 1 - SSresid/SStotal;


%% plot mean RH per night vs landings 
close 

figure()
scatter(mean_RH, landing_total)
ylabel("Total landings")
xlabel("Mean relative humidity (%)")
title("RH VS Landings")
h=lsline;
set(h,"LineWidth",3)

xlim([89 98])

txt=['R2 = ' num2str(rsq)];
t=text(min(mean_RH)+0.5,max(landing_total)-100,txt);
t.FontSize=16;

a = get(gca,'XTickLabel') 

set(gca,'XTickLabel',a,'fontsize',16)
set(gca,'TickDir','out'); % The only other option is 'in'


%% 6 person analysis from python tracking

data = load('data_6-human.mat');

data_ALL = data.data_ALL;


%% scatter percentages

human = cell2mat(data_ALL(5,:));

human_ID = unique(human);

percentage = cell2mat(data_ALL(6,:));

means_per = nan(1, length(human_ID));
medians_per = nan(1, length(human_ID));
dev_per = nan(2, length(human_ID));

all_dates = string(cellfun(@(x) unique(x), data_ALL(4,:),'UniformOutput',false));

dates_exp = unique(all_dates); % get dates

percentage_human = {};


for i=1:length(human_ID);
    
    curr_hum_ID = human == human_ID(:,i);
    
    means_per(1,i) = mean(percentage(:,curr_hum_ID));
    
    dev_per(1,i) = ste(percentage(:,curr_hum_ID),0,2); % get the standard error of the mean
    
    dev_per(2,i) = std(percentage(:,curr_hum_ID),0,2); % get the standard error of the mean
    
    
    medians_per(1,i) = median(percentage(:,curr_hum_ID));
    
    
    percentage_human{1,i} = percentage(:,curr_hum_ID)';
   
    
end

figure()

for i = 1:length(human_ID);
    
    night_ID = strcmp(dates_exp(1,i), all_dates);
    per_night = percentage(:,night_ID);
    
    humans_night = cell2mat(data_ALL(5,night_ID));
    
    [humans_sort, sort_ID] = sort(humans_night);
    
    per_sort = per_night(sort_ID);
    plot(humans_sort,per_sort, '-o')

    hold on 
    
end



scatter(human_ID, means_per,'*')

hold on

e = errorbar(human_ID, means_per,dev_per(1,:))

e.LineStyle = 'none';

ylabel('Landings percentage')
xlabel('human')

xlim([0.5, (length(human_ID)+0.5)]);

box off

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16)
set(gca,'TickDir','out');

legend(dates_exp)


%% scatter total

landings_all = cell2mat(data_ALL(2,:));

figure()

scatter(human, landings_all)

means_landings = nan(1, length(human_ID));
medians_landings = nan(1, length(human_ID));
dev_landings = nan(2, length(human_ID));

landings_humans = {};

for i=1:length(human_ID);
    
    curr_hum_ID = human ==human_ID(:,i);
    
    means_landings(1,i) = mean(landings_all(:,curr_hum_ID));
    
    dev_landings(1,i) = ste(landings_all(:,curr_hum_ID),0,2); % get the standard error of the mean
    
    dev_landings(2,i) = std(landings_all(:,curr_hum_ID),0,2); % get the standard error of the mean
    
    medians_landings(1,i) = median(landings_all(:,curr_hum_ID));
    
    landings_humans{1,i} = landings_all(:,curr_hum_ID)';
    
end

figure()

for i = 1:length(human_ID);
    
    night_ID = strcmp(dates_exp(1,i), all_dates);
    landings_night = landings_all(:,night_ID);
    
    humans_night = cell2mat(data_ALL(5,night_ID));
    
    [humans_sort, sort_ID] = sort(humans_night);
    
    landings_sort = landings_night(sort_ID);
        
    plot(humans_sort,landings_sort, '-o')

    hold on 
    
end



scatter(human_ID, means_landings,'*')

hold on

e = errorbar(human_ID, means_landings,dev_landings(1,:));

e.LineStyle = 'none';

ylabel('Total landings')
xlabel('human')

xlim([0.5, (length(human_ID)+0.5)]);


box off

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16)
set(gca,'TickDir','out');

legend(dates_exp)

%% 2022 weather dat analysis

weather = load("weather_2022.mat");

weather = weather.hobo_per_night; % load all hobo data

%% get temperature

temperature = weather(4,:);

exp_temp = cellfun(@(x) x(7:end), temperature, 'UniformOutput',false);

mean_temp = cellfun(@(x) mean(cell2mat(x),1), exp_temp, 'UniformOutput',false);

landing_total=nan(1,length(dates_exp));

for i=1:length(dates_exp);
    
    date_ID=cell2mat(cellfun(@(x) strcmp(x,dates_exp(1,i)), data_ALL(4,:),'UniformOutput',false));
   
    landings_date = landings_all(date_ID);
    landing_total(1,i)=sum(landings_date);

end


mean_t=cell2mat(mean_temp);


%% calculate R2 temp


p=polyfit(mean_t,landing_total,1);

yfit = polyval(p,mean_t);
yresid = landing_total - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(landing_total)-1) * var(landing_total);
rsq = 1 - SSresid/SStotal;

%% plot mean temeprature per night vs landings 
close 

figure()
scatter(mean_t, landing_total)
ylabel("Total landings")
xlabel("Mean temperature (°C)")
title('Temeprature VS Landings')
h=lsline;
set(h,"LineWidth",3)
txt=['R2 = ' num2str(rsq)];

xlim([16 22])

t=text(min(mean_t)+0.2,max(landing_total)-100,txt);
t.FontSize=16;

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16)
set(gca,'TickDir','out'); % The only other option is 'in'




%% get relative humidity

RH = weather(5,:);

exp_RH=cellfun(@(x) x(7:end), RH, 'UniformOutput',false);

mean_RH=cellfun(@(x) mean(cell2mat(x)), exp_RH, 'UniformOutput',false);
mean_RH=cell2mat(mean_RH);

%% calculate R2 RH

p=polyfit(mean_RH,landing_total,1);

yfit = polyval(p,mean_RH);
yresid = landing_total - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(landing_total)-1) * var(landing_total);
rsq = 1 - SSresid/SStotal;

%% plot mean RH per night vs landings 
close 

figure()
scatter(mean_RH, landing_total)
ylabel("Total landings")
xlabel("Mean relative humidity (%)")
title("RH VS Landings")
h=lsline;
set(h,"LineWidth",3)

xlim([89 98])

txt=['R2 = ' num2str(rsq)];
t=text(min(mean_RH)+0.5,max(landing_total)-100,txt);
t.FontSize=16;

a = get(gca,'XTickLabel') 

set(gca,'XTickLabel',a,'fontsize',16)


%% load lab data

clc
clear all

data = load('data_lab.mat');
data = data.data_AllC;


%% select lab data data

exp_type = data(6,:);

experiment_list={"bugdorm"};


selection= listdlg("ListString",experiment_list,"SelectionMode", "single");
experiment_type=experiment_list{1,selection};


type_ID=cell2mat(cellfun(@(x) strcmp(x,experiment_type), exp_type ,'UniformOutput',false));

data_type=data(:,type_ID);

type=unique([data_type{8,:}]);

%% landings

groupNumbers=nan(1,size(data_type,2));
exp_pos=data_type(8,:);
landings=cell2mat(data_type(1,:));
mean_landings=nan(1,length(type));
median_landings=nan(1,length(type));

for i=1:length(type);
    
    exp_ID=cell2mat(cellfun(@(x) strcmp(x,type(1,i)), exp_pos,'UniformOutput',false));
    groupNumbers(exp_ID)=i;
    
    mean_landings(1,i)=mean(landings(exp_ID));
    median_landings(1,i)=median(landings(exp_ID));    
    
end

landings_cell={};


for i=1:length(unique(groupNumbers));
    
    curr_type_ID=groupNumbers==i;
    landings_cell{1,i}=landings(curr_type_ID)';
end

std_groups=cellfun(@(x) std(x), landings_cell,'UniformOutput',false);
ste_groups=cellfun(@(x) ste(x,1,1), landings_cell,'UniformOutput',false);

%% plot mean landings
close

figure()
scatter(groupNumbers,landings)
hold on 
plot((1:1:length(type)),mean_landings,'k*', 'MarkerSize', 8, 'LineWidth', 2)

errorbar(mean_landings,cell2mat(ste_groups),'ko')

xticks([1:1:5])
xticklabels(type)
xlim([0.5 length(type)+0.5])
ylabel('Number of landings')
title('Landings with mean')

a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16)

% ylim([-10 max(landings)+10])

% breakyaxis([60 300])

hold off



%% Landing permutation test

landing_stats_mean={};
gn=num2cell(unique(groupNumbers));
 

 [pVals, pVals_C, H] = EST_multipleTests(landings_cell,'title',experiment_type,...
    'repNb',20000,'correctionType','benjami_P','fileOutput','./mean_percentage_stats.txt',...
    'alpha',0.05,'groupNames',type,...
    'testFunc','EST_permTest_meanDiff');










