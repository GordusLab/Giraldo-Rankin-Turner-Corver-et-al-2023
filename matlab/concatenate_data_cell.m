%% Concatenate semi-field data analyzed with ivTrace

% this script will find the trace_data.mat which contains the data for each
% trajectory file obtained from ivTrace. These files can be created with
% the script proccess_trajectories.m. 

% Each trace_data.mat file contains a cell with the following information:

% The first line contains the total landings
% The second line contains the individual trajectories
% The third line contains the occupancy through time 
% The fifth line contains the experiment type information
% The sixth line contains the treatment (e.g. CO2, Human1, control)
% The seventh line contains the position relative to the cage
% the eigths line contains the date of the experiment

% The processed data from each trajctory file is then concatenated to a
% single cell (data_AllC) for further analysis. 


Folder=cd;
FileList = dir(fullfile(Folder, '**/**/**/', 'trace_data.mat')); % find the trace_data.mat files and list thieir directories

paths={FileList.folder}'; % obtain the paths for each file
real_paths={};
count=1;

for i=1:3:length(paths);
    
    real_paths{count,1}=paths{i,1};
    count=count+1;
end

%% load data onto cell

data_AllC={}; % create an empty cell to store the data

for i=1:length(real_paths); % loop theough each folder, load the trace.mat cell and add it to the data_AllC cell
    
    curr_data=load([real_paths{i,1} '/trace_data.mat']);
    data_add=curr_data.trace_data;
    data_AllC(1:length(data_add),i)=data_add(1:end,1);
    
end

%%

save('data_AllC.mat', 'data_AllC','-v7.3') % save the cell