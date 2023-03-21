%% import data

% This section imports the full.txt file from ivTrace. Information on frame
% rate, experiment type, treatment, position in cage and date have to be provided
% for the current trajectory being analyzed

clc
clear all

answer = inputdlg('Maximum number of regions'); % this is necessary for matlab to fill the empty spaces with NaN. The maximum number of regions = (maximum number of columns in .txt file - 1)/5
answer = cell2mat(answer);
answer = str2num(answer);
fps = cell2mat(inputdlg('Recording frame rate')); % movies were recorded at 10fps
fps = str2num(fps);

number_of_columns=(answer*5)+1;
columns=[];

date_exp = cell2mat(inputdlg('Date of the experiment')); % input the date of the experiment

experiment_list = {"CO2",  "1-human", "2-human", "Lab"}; % select the typer of experiment


selection= listdlg("ListString",experiment_list,"SelectionMode", "single");
experiment_type=experiment_list{1,selection};

if selection ==1
    
    position_list={"CO2", "1", "2", "3", "4", "5", "6", "7"}; % list treatment options (the numbers are the control positions relative to CO2);
    pos_selection= listdlg('PromptString','Select a treatment',...
        "ListString",position_list,"SelectionMode", "single"); % select the treatment
    position=position_list{1,pos_selection};
    
    pos_list={"1", "2", "3", "4", "5", "6","7","8"}; % list of positions relative to the cage
    position_selection= listdlg('PromptString','Select position in the cage',"ListString",pos_list,"SelectionMode", "single"); % select the position relative to the cage
    real_pos=pos_list{1,position_selection};
end




if selection ==2
    
    position_list={"Human 1", "1", "2", "3", "4", "5", "6","7"}; %list treatment options (the numbers are the control positions relative to human1);
    pos_selection= listdlg('PromptString','Select a treatment',"ListString",position_list,"SelectionMode", "single"); % select treatment
    position=position_list{1,pos_selection};
    
    pos_list={"1", "2", "3", "4", "5", "6","7","8"}; % list of positions relative to the cage
    position_selection= listdlg('PromptString','Select position in the cage',"ListString",pos_list,"SelectionMode", "single");% select the position relative to the cage
    real_pos=pos_list{1,position_selection};
    
end

if selection ==3 
    
    position_list={"Human 1", "Human 2", "1", "2", "3", "4", "5", "6","7"};  %list treatment options (the numbers are the control positions relative to human1);
    pos_selection= listdlg('PromptString','Select a treatment',"ListString",position_list,"SelectionMode", "single");  % select the treatment
    position=position_list{1,pos_selection};
    
    pos_list={"1", "2", "3", "4", "5", "6","7","8"}; % list of positions relative to the cage
    position_selection= listdlg('PromptString','Select position in the cage',"ListString",pos_list,"SelectionMode", "single");  % select the position relative to the cage
    real_pos=pos_list{1,position_selection};
    
end

if selection ==4
    
    pos_list={"CO2_heat", "CO2_alone", "Air_Heat", "Heat", "Air_alone", "no_stimulus"};
    position_selection= listdlg('PromptString','Select a treatment',"ListString",pos_list,"SelectionMode", "single");
    real_pos=pos_list{1,position_selection};
    
end


if answer==0
    
    full_data=load('full.txt');
    
else
    
    for i=1:number_of_columns;
        columns=[columns, '%f '];
        
    end
    
    full_data=readtable('full.txt','ReadVariableNames',false,'HeaderLines',0,'Format',columns);
    full_data=table2array(full_data);
end


%% Count number of landings

if answer==0
    
    total_landings=0;
    
else
    
    nan_IDX=nan(length(full_data), answer);
    
    count=2;
    
    for i=1:answer;
        
        nan_IDX(:,i)=~isnan(full_data(:,count)); % from the x position columns, identify frames with a region detected (i.e. not NaN)
        count=count+5;
        
    end
    
    occupationID_diff=diff(nan_IDX); % Find frames where a region appears (i.e. change from NaN to a number)
    
    occupation_startID=occupationID_diff==1; % find where a new region appears
    occupation_stopID=occupationID_diff==-1; % find when a reagion dissapears
    
    total_landings=sum(sum(occupation_startID)); % the total landings is the sum of all the instances where a new region appears
        
    
end

total_landings
%% target occupancy

nan_IDX = nan(length(full_data), answer);

count=2;

for i = 1:answer; % loop through the 7 columns containing the x position
    
    nan_IDX(:,i)=~isnan(full_data(:,count));  % in each column find the frames that had a region detected
    count=count+5;

end

total_per_time = sum(nan_IDX,2); % calculate the total number of regions per frame

%% plot target occupancy

close
figure(1)
plot(total_per_time)
xticks(0:length(total_per_time)/12:length(total_per_time));
xlim([0,length(total_per_time)+900]);
ylim([0,max(total_per_time)+1]);
xticklabels(["22:00", "22:30", "23:00", "23:30", "00:00", "00:30", "01:00", "01:30", "02:00", "02:30", "03:00", "03:30", "4:00"])
title("Target occupancy")
ylabel("Number of mosquitoes on target")
xlabel("Time of day")


%% Extract individual traces 

% individual trajectories are extracted by finding frames where a new
% region is detected and frames where a region dissapears. 
% each individual trajectory is stored in the cell traces_cell

frames=full_data(:,1);
start_frames=frames(occupation_startID(:,1));
stop_frames=frames(occupation_stopID(:,1));
duration_cell={};
traces_cell={};


for i=1:length(start_frames);
    
    if i>length(stop_frames)==1
        
        interval=(start_frames(i,1)+1):frames(end,1);
        nan_interval=nan_IDX(interval,:);
        total_interval=sum(nan_interval);
        sum_ID=sum(total_interval>0);
    
    else
        interval=(start_frames(i,1)+1):stop_frames(i,1);
        nan_interval=nan_IDX(interval,:);
        total_interval=sum(nan_interval);
        sum_ID=sum(total_interval>0);
    end
    
    if sum_ID==1
        
        real_duration=total_interval(1,1)/fps;
        curr_data=full_data(interval+1,:);
        
        duration_segment={};
        traces_segment={};
        traces_segment={curr_data(:,2:6)};
        traces_segment{1,1}(:,6)=curr_data(:,1);
        duration_segment(1,:)= {real_duration};
    
    elseif sum_ID==0 & size(interval,1)==1
        
        total_landings=total_landings-1;
    
    else
       
        curr_data=full_data(interval+1,:); % get section to analyze
        
        interval_diff=diff(nan_interval);
        stop_interval=interval_diff<0;  % find where animal leaves
        stop_interval(end,1)=1; % define the last line as a leave for the first column
        stop_all=logical(sum(stop_interval,2));
        start_interval=interval_diff>0; % find when animal arrives
        start_all=logical(sum(start_interval,2));
        stop_current_pos=find(stop_all==1); % index where stop events happen
        start_current_pos=find(start_all==1); % index where start events happen
        start_current_frames=find(start_all==1);    % index where start events happen
        total_animals=sum(start_all);   % total animals in the section
        
        
        for x=2:sum(sum(start_interval)>0); % set the new start ID

            new_start=start_interval(:,1:x-1)+stop_interval(:,x);
            start_interval(:,1:x-1)=new_start;
            
        end
        
        
        start_ID_vector=[];
        startint_sum=sum(start_interval,2);
        start_ident=find(startint_sum==1);
        
        [start_indent,start_ID_vector]=find(start_interval==1,length(startint_sum));
        
        [start_repeated,start_sorted]=sort(start_indent);
        start_ID_vector=start_ID_vector(start_sorted);
        
        
        for y=2:sum(sum(stop_interval)>0); % set the new stop ID

            new_stop=stop_interval(:,1:y-1)+stop_interval(:,y);
            stop_interval(:,1:y-1)=new_stop;
            
        end
        
        stop_repeated=[]; % create a stop vector with the new IDs
        
        for a=1:sum(sum(stop_interval)>0);
            
            stop_currID=find(stop_interval(:,a)==1);
            stop_repeated=[stop_repeated; stop_currID];
            
        end
        
        stop_repeated=stop_repeated-1;
        
        curr_data_sections={};
        
        column_count=2; %make it a variable
        
        selec_data=curr_data(:,2:end);
        
        stop_repeated=sort(stop_repeated);
        
        for n=1:length(stop_repeated); % separate segments
            
            start_section=start_repeated(n,1);
            stop_section=stop_repeated(n,1);
            frames_segment=curr_data(start_section:stop_section,1);
            
            select_column=(((start_ID_vector(n,1))*5)-4):(((start_ID_vector(n,1))*5));
            curr_data_sections{1,n}=selec_data(start_section:stop_section,select_column);
            curr_data_sections{1,n}(:,6)=frames_segment;
            
            if n==1 % define if it segment is a continuation (1), or not (0)
                curr_data_sections{2,n}=0;
            else
                start_prev=start_section-1;
                if isnan(selec_data(start_prev,select_column))==1;
                    curr_data_sections{2,n}=0;
                    
                else
                    curr_data_sections{2,n}=1;
                end
            end
            
            curr_data_sections{3,n}=curr_data_sections{1,n}(1,1); %start x position
            curr_data_sections{4,n}=curr_data_sections{1,n}(end,1);%end x position
            
        end
                
        final_segment={};
        
        
        for z=1:size(curr_data_sections,2);
            
            if curr_data_sections{2,z}==0
                
                final_pos=size(final_segment,2)+1;
                final_segment{1,final_pos}=curr_data_sections{1,z};
                final_segment{2,final_pos}=curr_data_sections{3,z};
                final_segment{3,final_pos}=curr_data_sections{4,z};
               
            
            else
                
                xpos_diff=abs((cell2mat(final_segment(3,:)))-curr_data_sections{3,z});
                [~,min_ID]=min(xpos_diff);
                segment_trace=final_segment{1,min_ID};
                new_segment=[segment_trace;curr_data_sections{1,z}];
                final_segment{1,min_ID}=new_segment;
                final_segment{3,min_ID}=curr_data_sections{4,z};
                
            
            end
            
           
            
        end
        
         duration_segment={};
         traces_segment={};
         traces_segment=final_segment(1,:);
         duration_segment(1,:)= cellfun(@(x) (length(x))/10,final_segment(1,:),'UniformOutput',false);
            
            
    
    end
    
    duration_cell = [duration_cell, duration_segment];
    traces_cell = [traces_cell, traces_segment];
    
end

%% Make final cell


% Save the data in a cell file "trace_data.mat"

trace_data={};

trace_data{1,1}=total_landings; % store the total number of landings

trace_data{1,2}='total landings';

if total_landings ==0 % if there are no landings
    
    for i=2:5;
        trace_data{i,1}=nan;
    end
    
   
    trace_data{2,2}='individual traces'; 
    trace_data{3,2}='raw data';
    trace_data{4,2}='occupation through time';
    trace_data{5,1}=experiment_type; % define the experiment type
    trace_data{5,2}='Experiment type';
    
    
    trace_data{6,1}=position; % define the position relative to the stimulus
    trace_data{6,2}='treatment';
    
    trace_data{7,1}=real_pos;
    trace_data{7,2}='Position';
    
    trace_data{8,1}=date_exp;
    trace_data{8,2}='Date';
    
    
elseif selection == 4 % if the experiment was a lab experiment
    
    trace_data{2,1}=total_per_time; % store occupation through time
    trace_data{2,2}='occupation through time';
   
    trace_data{3,1}=full_data; % store the entire raw data
    trace_data{3,2}='raw data';
    
    trace_data{4,1}=experiment_type; % define the experiment type
    trace_data{4,2}='Experiment type';

    trace_data{5,1}=date_exp;
    trace_data{5,2}='Date';
    
    trace_data{6,1}=real_pos;
    trace_data{6,2}='Treatment';

    
else % if it the experiment was a semi-field experiment and one or more region were detected
    
    
    trace_data{2,1}=traces_cell; % store each individual trace
    trace_data{2,2}='individual traces';
    
    
    trace_data{3,1}=full_data; % store the entire raw data
    trace_data{3,2}='raw data';
    
    trace_data{4,1}=total_per_time; % store occupation through time
    trace_data{4,2}='occupation through time';
    
    trace_data{5,1}=experiment_type; % define the experiment type
    trace_data{5,2}='Experiment type';
    
    
    trace_data{6,1}=position; % define the position relative to the stimulus
    trace_data{6,2}='treatment';
    
    trace_data{7,1}=real_pos;
    trace_data{7,2}='Position';
    
    trace_data{8,1}=date_exp;
    trace_data{8,2}='Date';
end




save('trace_data.mat','trace_data')  % save the cell

% the trace_data.mat file is a cell that contains the data extracted from
% each trajectory. 
% The first line contains the total landings
% The second line contains the individual trajectories
% The third line contains the occupancy through time 
% The fifth line contains the experiment type information
% The sixth line contains the treatment (e.g. CO2, Human1, control)
% The seventh line contains the position relative to the cage
% the eigths line contains the date of the experiment










