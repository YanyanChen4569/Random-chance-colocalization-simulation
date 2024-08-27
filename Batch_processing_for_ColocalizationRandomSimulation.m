%%% Batch simulation for random chance colocalization 
%%% Yanyan Chen (08/27/2024; Columbia University; yc4569@cumc.columbia.edu)
close all; 
clear all; clc;
%%
% Import the Summary table of the Coordinates from two channels

% file_path='/Users/chenyanyan/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Confocal Facility YC/Image analysis/Image analysis for Alondra Schweizer Burguete/Colocolization_Simulation/CS29iALS_ISO images/';
% T = readtable([file_path '081524_CS29iALS_ISO peaks merged Results.csv']);

file_path='/Users/chenyanyan/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Confocal Facility YC/Image analysis/Image analysis for Alondra Schweizer Burguete/Colocolization_Simulation/323-2 images/';
T = readtable([file_path '082224_323-2 Images 1,2,3,4,6 merged Results.csv']);


% Acquire the mask file name from the table
% Change to the mask file name
MaskFileName=regexprep(T.filename,'.tif_Results.csv','_Mask.tif');
% Aquire number of masks
Num_Masks=length(MaskFileName);

% Extract the MaskName,Num_Ch1,Num_Ch2 from the table 
% And feed these parameters to the random simulation
% Ch1: FUS; Ch2: DNAJB6

Array_FUSNum=T.NumberOfFUSMarkers;
Array_DNAJB6Num=T.NumberOfDNAJB6Markers;

Num_iteration=20;  % Num_iteration is the round of the simulations in each mask
for i=1:1:Num_Masks
    MaskName=MaskFileName{i};
    Num_Ch1=Array_FUSNum(i);
    Num_Ch2=Array_DNAJB6Num(i);
    [Num_Colo_multiple_simulation]=ColocolizationAnalysis(MaskName,Num_Ch1,Num_Ch2,Num_iteration,i);
    % Store the number of colocalization from multiple simulation in to an
    % array
    Num_Colo_multiple_simulation_array{i}=Num_Colo_multiple_simulation;
    for iter=1:1:Num_iteration
        AvergNumColoWithIter(iter)=mean(Num_Colo_multiple_simulation([1:iter]));
    end
    AvergNumColoWithIter_array{i}=AvergNumColoWithIter;
    PrecentColo=AvergNumColoWithIter/Num_Ch1;
    PrecentColo_array{i}=PrecentColo;

    AvergNumColoWithIter20(i)=AvergNumColoWithIter(end);
    PrecentColo_withIter20(i)=PrecentColo(end);

end 

%% Add the simulated results to the orignial table and export

T.PrecentColo_experiment=T.NumberOfFUSWithin0_071UmOfDNAJB6./T.NumberOfFUSMarkers;
T.NumberofFUScolocalized_simulation=Num_Colo_multiple_simulation_array';
T.AvergNumColoWithIter20_simulation=AvergNumColoWithIter20';
T.PrecentColo_withIter20_simulation=PrecentColo_withIter20';

% Save the table
 writetable(T,'082224_323-2 Images 1,2,3,4,6 merged Results_withSimulation_complete dataset.csv')
%% Plot the precent vs iteration

figure()
for mm=1:1:Num_Masks
plot(PrecentColo_array{1, mm},'LineWidth',2)
hold on
end
set(gca, 'YScale', 'lin', 'YLim', [0 0.2], 'xlim',[1 100], 'FontSize',16);
xlabel('Number of simulations');
ylabel('Precent of FUS with Colocalization');
box off
%% Plot the precent of colocalization between experiment and simulation

x1=ones(Num_Masks,1);
x2=2*ones(Num_Masks,1);
figure()
scatter(x1,T.PrecentColo_experiment)
hold on
scatter(x2,T.PrecentColo_withIter20_simulation)
set(gca, 'YScale', 'lin', 'YLim', [0 1], 'xlim',[0 3], 'FontSize',16);
xticklabels({'' ,'', 'Experiment', '','Simulation',''}); set(gca,'TickLabelInterpreter', 'tex')

%% Totoal precent 

total_precent_exp=sum(T.NumberOfFUSWithin0_071UmOfDNAJB6)/sum(T.NumberOfFUSMarkers)
total_precent_simul=sum(T.AvergNumColoWithIter20_simulation)/sum(T.NumberOfFUSMarkers)

%% Statistical difference

file_path1='/Users/chenyanyan/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Confocal Facility YC/Image analysis/Image analysis for Alondra Schweizer Burguete/Colocolization_Simulation/';
stack_data_table = readtable([file_path 'Num_colocalization_per_stack_combined dataset.csv']);
experiment_stack=stack_data_table.Num_colocalization_exp_per_stack;
simulation_stack=stack_data_table.Num_colocalization_simul_per_stack;
[h3,p3] = ttest2(experiment_stack,simulation_stack,'Tail','right')

[Num_stacks,Num_colum]=size(stack_data_table)
x1=ones(Num_stacks,1);
x2=2*ones(Num_stacks,1);
figure()
scatter(x1,experiment_stack)
hold on
scatter(x2,simulation_stack)
set(gca, 'YScale', 'lin','XLim',[0 3],'FontSize',16);
xticklabels({'' ,'', 'Experiment', '','Simulation',''}); set(gca,'TickLabelInterpreter', 'tex')
ylabel('Number of FUS with Colocalization');

