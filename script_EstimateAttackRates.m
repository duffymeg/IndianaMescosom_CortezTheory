% Script for estimating attack rates on clones
% AUTHORS:  Michael Cortez, Florida State University, mcortez@fsu.edu
%           Adapted code provided by Spencer Hall, Indiana University, sprhall@indiana.edu
% DATE:     June 2021
% PROJECT:  Predation trials with Chaoborus
% PURPOSE:  Estimate attack rate parameter for Type 1 functional response

% INPUTS:   "Data_clone" -- column 1 is Daphnia number provided; 
%                           column 2 is number eaten

% METHOD:   Step 1: Estimate attack rate (ahat) by minimizing the negative  
%                   log likelihood assuming binomially distributed data
%                   Value is analytically computed and accounts for
%                   prey density decreasing over time
%           Step 2: Bootstrap data -- this code preserves the number of
%                   replicates in each treatment
%           Step 3: Estimate percent confidence intervals (pCI)
%           Step 4: Plot results

% OUTPUTS:  all outputs are to screen
%           ahat -- estimated attack rates (1/hr/predator)
%           pCI_lower -- lower bound for percent confidence intervals
%           pCI_upper -- upper bound for percent confidence intervals
%           Plot Panel A -- Data, estimated functional response, and bcCI
%           Plot Panel B -- histogram of bootstrapped attack rates

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Experiment Parameters

num_clone = 9;                 % Number of clones
T = 16;                         % Length of experiment is 16 hours
y = 1;                          % number of predators in all treatments
Dens = [1 2 5 10];              % prey numbers in different treatments


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import Data

Data_BD0542=readmatrix('PredationData/BD05-42.xlsx','NumHeaderLines',1);
Data_BD0846=readmatrix('PredationData/BD08-46.xlsx','NumHeaderLines',1);
Data_BD1964=readmatrix('PredationData/BD19-64.xlsx','NumHeaderLines',1);
Data_CB2468=readmatrix('PredationData/CB24-68.xlsx','NumHeaderLines',1);
Data_DW2258=readmatrix('PredationData/DW22-58.xlsx','NumHeaderLines',1);
Data_DW2975=readmatrix('PredationData/DW29-75.xlsx','NumHeaderLines',1);
Data_IL1443=readmatrix('PredationData/IL14-43.xlsx','NumHeaderLines',1);
Data_ML3082=readmatrix('PredationData/ML30-82.xlsx','NumHeaderLines',1);
Data_ML3284=readmatrix('PredationData/ML32-84.xlsx','NumHeaderLines',1);

Names = ['Data_BD0542'; 'Data_BD0846'; 'Data_BD1964'; 'Data_CB2468';...
    'Data_DW2258'; 'Data_DW2975'; 'Data_IL1443'; 'Data_ML3082'; 'Data_ML3284'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute best fitting attack rate and confidence interval for each clone

ahat = zeros(num_clone,1);      % estimated attack rate 1/predator/day
pCI_alpha = 0.05;               % alpha for percent confidence interval (pCI)
pCI_lower = zeros(num_clone,1); % lower bound for pCI
pCI_upper = zeros(num_clone,1); % lower bound for pCI

for ii = 1:num_clone
   
    % Step 1: Compute ML estimate
    
        data = eval(Names(ii,:));     % Data for clone ii
        SumN = sum(data(:,1));  % total number of each clone provided
        SumE = sum(data(:,2));  % total number of each clone eaten
        ahat(ii) = -log(1-SumE/SumN)/T/y; % attack rate 


    % Step 2: Create Bootstrapped Data
    
    % A. Separate data into different treatment groups
        cat1 = find(data(:,1)==Dens(1));
        cat2 = find(data(:,1)==Dens(2));
        cat3 = find(data(:,1)==Dens(3));
        cat4 = find(data(:,1)==Dens(4));
    
    % B. Bootstrap data
        
        BS = 1000;				% bootstrap number
        ncat1=length(cat1);		% number of replicates for treatment 1
        ncat2=length(cat2);		% number of replicates for treatment 2
        ncat3=length(cat3);		% number of replicates for treatment 3
        ncat4=length(cat4);     % number of replicates for treatment 4

        Bparsb = zeros(BS,1);	% initalize holder of all of the bootstrapped CIs

        for kk=1:1:BS			% from 1 to number of boostraps
                                            
           order=round((ncat1-1)*rand(ncat1,1)+1);	% create bootstrapped sample
           data_cat1=data(cat1(order),:);      % create bootstrapped data for treatment 1

           order=round((ncat2-1)*rand(ncat2,1)+1);	% create bootstrapped sample
           data_cat2=data(cat2(order),:);		% create bootstrapped data for treatment 2

           order=round((ncat3-1)*rand(ncat3,1)+1);	% create bootstrapped sample
           data_cat3=data(cat3(order),:);		% create bootstrapped data for treatment 3

           order=round((ncat4-1)*rand(ncat4,1)+1);	% create bootstrapped sample
           data_cat4=data(cat4(order),:);		% create bootstrapped data for treatment 4
           
           rdata=[data_cat1; data_cat2; data_cat3; data_cat4];	    % Combined bootstrapped data
           
           % Compute ML estimate using bootstrapped data
           SumN = sum(rdata(:,1));  % total number of each clone provided
           SumE = sum(rdata(:,2));  % total number of each clone eaten
           Bparsb(kk) = -log(1-SumE/SumN)/T/y;
        end
        
        
        
    % Step 3: Sort bootstrapped data and compute percent confidence intervals
        
        Bparsb = sort(Bparsb);
        pCI_lower(ii) = Bparsb(floor(pCI_alpha/2*BS));
        pCI_upper(ii) = Bparsb(ceil((1-pCI_alpha/2)*BS));
        
        
    % Step 4: Plot Results
    
        figure('Name',Names(ii,:));
        
        subplot(1,2,1) % Plot of data and best fitting Type 1 
        hold on
        plot(data(:,1),data(:,2),'o')   % Data
        xvals = 0:10;
        % Plot best fitting Type 1 FR
        yvals = xvals.*(1-exp(-ahat(ii)*y*T));
        plot(xvals,yvals,'-k')
        % Plot Type 1 with lower bcCI
        yvals = xvals.*(1-exp(-pCI_lower(ii)*y*T));
        plot(xvals,yvals,'--k')  % Type 1 with lower bcCI
        % Plot Type 1 with upper bcCI
        yvals = xvals.*(1-exp(-pCI_upper(ii)*y*T));
        plot(xvals,yvals,'--k')  % Type 1 with upper bcCI

        xlabel('Prey Available'); ylabel('Prey Eaten');
        axis([0 11 0 (max(data(:,2))+2) ]);

        
        subplot(1,2,2) % histogram of the bootstrapped values
        hold on
        hist(Bparsb, 20)    % Bootstrapped values
                            % estimated value is denoted in red
        plot([ahat(ii) ahat(ii)],[0 max(hist(Bparsb, 20))],'-r')    
        xlabel('attack rates') 
end
