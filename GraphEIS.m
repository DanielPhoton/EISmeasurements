%% clear
clear all
clc
%% import the CV, EIS, Voc data
% note: CV is cyclic voltametery and is just a JV curve
% 
start_path = '~/Documents/MATLAB/NanoTechFiles/EIS measurements/';

dialog_title = 'Select EIS data folder to graph.';
EIS_folder = uigetdir(start_path,dialog_title); % Asks for the EIS folder

cd(EIS_folder); % change to that directory

folderparts = regexp(pwd, '/', 'split'); % get folder names

listing = dir('*.DTA'); % get list of all the files that end with .dta

numfiles = length(listing);
tag = cell(1,numfiles); % measurement type (CV, EIS, or Voc)
a=0;b=0;c=0;% counters

for K = 1:numfiles
    [tag{K}] = importEIStag(listing(K).name); %imports the measurement tag
    
    if strcmp(tag{K},'EISPOT') % if the tag matches 'EISPOT'
        a=a+1;
        [Freq{a},Zreal{a},Zimag{a},Zsig{a},Zmod{a},Zphz{a},Idc{a},Vdc{a}] = importEISfile(listing(K).name); % Imports all data
    
    elseif strcmp(tag{K},'CV')
        b=b+1;
        [CVT{b},CVV{b},CVJ{b},Vu{b},Sig{b},CVAch{b}] = importEISCV(listing(K).name);% Imports all data
    
    elseif strcmp(tag{K},'CORPOT')
        c=c+1;
        [VocT{c},VocV{c},Vm{c},VocAch{c}] = importEISVoc(listing(K).name);% Imports all data
    
    else
        disp('Error: Unknown tag');% Displays error if measurement type is new or unknown
    
    end
    
end
%% Cyclical Voltametry parameter (CV Parameter)
% Gets the area of device
defaultarea=0.084; % cm^2
choicearea = menu('How do you want to enter the Area?','Load Matlab Data in folder','Enter Area manually');

if choicearea == 2 % 2 means enter manually
    areacell=inputdlg({['Please enter the area of the device in cm^{2}. ']}, 'Area' , [1 60]);
    area=str2num(areacell{1});

elseif choicearea==1 % 1 means Load Data
    try
        load('area.mat');
    catch
        area=defaultarea;
    end
else
    area=defaultarea;
end

save('area.mat','area'); % save area to file in the folder
% get SC parameters by calling the SCparameters function
if b ~=0 % with CV data
    scanerate = zeros(2,length(CVV));
    for K = 1:length(CVV)
        [Voc{K} Jsc{K} FF{K} Eff{K} scanrate(:,K)] = SCparameters(CVV{K},CVT{K},CVJ{K}.*1000./area); % note that .*1000./area is to get current density in mA/cm^2
    end
    % Get Medians for the CV SC parameters
    
    for K = 1:length(CVV)
        Vocmedian(1,K)=median(Voc{K}(1,:));Vocmedian(2,K)=median(Voc{K}(2,:));
        Jscmedian(1,K)=median(Jsc{K}(1,:));Jscmedian(2,K)=median(Jsc{K}(2,:));
        FFmedian(1,K)=median(FF{K}(1,:));FFmedian(2,K)=median(FF{K}(2,:));
        Effmedian(1,K)=median(Eff{K}(1,:));Effmedian(2,K)=median(Eff{K}(2,:));
    end
else % no CV data 
    for K = 1:2
        Vocmedian(1,K)=0;Vocmedian(2,K)=0;
        Jscmedian(1,K)=0;Jscmedian(2,K)=0;
        FFmedian(1,K)=0;FFmedian(2,K)=0;
        Effmedian(1,K)=0;Effmedian(2,K)=0;
    end
end
%% Plot EIS 

colors = {'b','r','k',[0 .6 0],'m',[.5 .6 .7],[.3 .3 .3]}; % Define colours 
choicescale = menu('Do you want to autoscale JV graph?','Yes','No');

if choicescale==0 || choicescale==2 %no
    choicescale=0;
    prompt = {'Enter low x limits:','Enter upper x limits:','Enter lower y limit:','Enter upper y limit:'};
    dlg_title = 'Axis limits';
    num_lines = 1;
    defaultans = {'-0.5','1.3','-26','20'};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    
    xlow = str2num(answer{1}); 
    xhigh = str2num(answer{2});
    ylow = str2num(answer{3}); 
    yhigh = str2num(answer{4});
elseif choicescale==1 %yes
    choicescale=1;
end

h=figure('Position',[50, 50, 1440, 900]);
suptitle([folderparts(length(folderparts)-1) folderparts(length(folderparts))]); %make a main title  that is based on the folder names

set(gca, 'Position', [1.5 1.5 .4 .4]);% move blank graph off screen
ha = tight_subplot(3,2,[.1 .05],[.1 .1],[.05 .01]) % make the subplot graphs at good resolution

axes(ha(1));set(gca, 'Position', [1.5 1.5 .4 .4]);% Move blank graph off screen
axes(ha(5));set(gca, 'Position', [1.5 1.5 .4 .4]);% Move blank graph off screen
axes(ha(3));set(gca, 'Position', [.05 .4 .445 .5]);

if b ~=0 % if we have JV data 
    for i =1:length(CVV)
        plot(CVV{i},CVJ{i}.*1000./area,'color',colors{i});hold on % Plots JV data
    end
end

for i =1:length(Vdc)
    scatter(median(Vdc{i}),median(Idc{i}.*1000./area),90,colors{i},'filled'); % Plots dots at the EIS measurement points.
end

plot([-100 100],[0 0],'color','k');plot([0 0],[-100 100],'color','k'); % Plot x and y axis
title('JV Curve'); grid on; ylim([ylow yhigh]); xlim([xlow xhigh]);
ylabel('Current Density (mA/cm^{2})'); xlabel('Voltage (V)');
axes(ha(2)); %switch to the next subplot

for i =1:length(Zphz)
    scatter(Freq{i},Zphz{i},36,colors{i},'filled');hold on % plot Phase part of the Bode graph
end

title('Bode Graph'); grid on;set(gca,'XScale','log');
ylabel('Phase (degrees)'); xlabel('Frequency (Hz)'); xlim([.4 1e6]);ylim([-90 0]);
axes(ha(4));

for i =1:length(Zmod)
    scatter(Freq{i},Zmod{i},36,colors{i},'filled');hold on % plot modulus part of the bode graph
end

title('Bode Graph'); grid on;set(gca,'XScale','log');set(gca,'YScale','log');
ylabel('Modulus (Ohms)'); xlabel('Frequency (Hz)'); xlim([.4 1e6]);ylim([10 10e4]);
axes(ha(6));

for i =1:length(Zmod)
    scatter(Zreal{i},-Zimag{i},36,colors{i},'filled');hold on % plot Nyquist data as Circles
    scatter(-Zreal{i},Zimag{i},36,colors{i},'filled','s'); % plot Negative data as Qquares  
end

title('Nyquist Chart'); grid on;set(gca,'XScale','log');set(gca,'YScale','log');
ylabel('-Z Imag (Ohms)'); xlabel('Z Real (Ohms)');  ylim([10e-3 10e4]);%xlim([1 1e6]);
axes(ha(3));

t = uitable; % make table for solar cell values
t.Data = {Vocmedian(1,1),Jscmedian(1,1),FFmedian(1,1),Effmedian(1,1);...
            Vocmedian(2,1),Jscmedian(2,1),FFmedian(2,1),Effmedian(2,1);...
            Vocmedian(1,2),Jscmedian(1,2),FFmedian(1,2),Effmedian(1,2);...
            Vocmedian(2,2),Jscmedian(2,2),FFmedian(2,2),Effmedian(2,2)};
t.Position = [122 204 477 95];
t.ColumnName = {'Voc (V)','Jsc (mA/cm^2)','FF','Eff (%)'};
t.RowName = {'Before EIS up','Before EIS down','After EIS up','After EIS down'};