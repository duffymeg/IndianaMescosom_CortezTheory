% Script for applying Geber Method to mesocosm experiment data

% AUTHOR:   Michael Cortez, Florida State University, mcortez@fsu.edu
% DATE:     March 2022
% PURPOSE:  Use Geber Method (Ellner et al. 2005, Ecology Letters) to
%           estimate effect sizes of ecology and evolution on 
%           reproduction number (R)

% INPUTS:   genotype.frequencies_MHC_converted.csv -- data table of clonal
%           frequencies in weeks 0, 2, 6, 9
%           meso.expt.densities_MHC_converted.csv -- data table of total
%           host densities and infected densities in weeks 1-9.
%           Order of clones is (1) BD05-42, (2) BD08-46	(3) BD19-64	
%           (4) CB24-68	(5) DW22-58	(6) DW29-75	(7) IL14-43	(8) ML30-82
%           (9) ML32-84

% METHOD:   Step 1: Define parameters estimated from other experiments 
%           Step 2: Import data
%           Step 3: Extract densities and frequencies for weeks 0 & 2
%           Step 4: Compute R values at each week and eco/evo effect sizes
%           Step 5: Plot estimated R values
%           Step 6: Plot eco and evo effect sizes
%           Step 7: Plot ratio of eco and evo effect sizes

% OUTPUTS:  R_msred -- Estimated R values at weeks 0 & 2
%           DeltaEco_1a -- Estimated effect sizes of ecology between weeks
%           0 and 2
%           DeltaEvo_1a -- Estimated effect sizes of evolution between weeks
%           0 and 2

% IMPORTANT: Only data for weeks 0 and 2 plotted/analyzed because predator
%           densities declined after.

clear all


%% Step 1: Define parameters 

    
    % per spore probability of infection
    p_i = [5.54 5.64 0 5.25 6.88 5.77 11.7 8.58 4.93]*10^(-4);
    
    % uptake rates in presence of spores (L/hr)
    fSi = [1.44 1.74 0.994 1.55 2.46 2.65 1.54 0.58 2.78]*10^(-4);
    fIi = 0.5*fSi;
    
    % attack rates on susceptible individuals(1/predator/hr)
    a_i = [0.0197 0.0134 0.0116 0.0062 0.0144 0.0082 0.016 0.0157 0.0096];
    
    % increase in attack rate on infected individuals
    omega = 1;
    
    % spore yield (spores/indiv.) 
    chi_i = [78231 36990 83035 70578 89909 146128 108163 41667 92610];
    
    % Disease-induced morality rate (1/hr)
    mu_i = 1/20/24;
    
    % spore degradation rate (1/day)
    delta = 0.0083;
    
    % R0_i(infinity) Competence for host
    R0i = p_i.*chi_i;
    
    % Reduction in spore burst size for P = 0, 0.1, 0.5, and 1
    xi_0pred = [1 1 1 1 1 1 1 1 1];
    xi_01pred = [0.2104 0.4099 0.4954 0.7929 0.3681 0.6828 0.3098 0.3200 0.6032];
    xi_05pred = [0.0088 0.0163 0.0213 0.0796 0.0144 0.0431 0.0121 0.0125 0.0310];
    xi_1pred=[0.0040 0.0059 0.0071 0.0188 0.0055 0.0116 0.0049 0.0050 0.0091];

    

%% Step 2: Import Data

    % A. Import total host density data 
        % Order of columns is 
        % (1) Parasites present: 1 = yes, 0 = no
        % (2) Predator density (1/L)
        % (3) tank
        % (4) Week
        % (5) Julian Date
        % (6) Host density (1/L)
        % (7) Infected host density (1/L)
        DenData=readmatrix(...
            'MesocosmData/meso.expt.densities_MHC_converted.csv',...
            'NumHeaderLines',1);

    % B. Import clone frequencies data 
        % Order of columns is 
        % (1) Parasites present: 1 = yes, 0 = no
        % (2) Predator density (1/L)
        % (3) tank
        % (4) Week
        % (5) Julian Date
        % (6) number of host sampled
        % (7-15) Freqs of (7) BD05-42, (8) BD08-46, (9) BD19-64, (10) CB24-68, 
        %     (11) DW22-58, (12) DW29-75, (13) IL14-43, (14) ML30-82, (15) ML32-84
        FreqData=readmatrix(...
            'MesocosmData/genotype.frequencies_MHC_converted.csv',...
            'NumHeaderLines',1);
        

%% Step 4: Extract Host densities and frequencies at weeks 0, 2, 6, 9

    % A. Extract all densities, sampling days, and predator densities
        %row is tank, column is time
        SampleDay = zeros(52,9)+NaN;
        HostDen = zeros(52,9)+NaN; 
        InfectDen = zeros(52,9)+NaN; 
        PredDen = zeros(52,1)+NaN;

        for ii = [1:24 27:32 35:52]
            % Find rows for tank ii 
            tankrow = find(DenData(:,3)==ii);
            % Sampling Day 
            SampleDay(ii,1:9) = DenData(tankrow,5)';
            % Host densities for tank ii
            HostDen(ii,1:9) = DenData(tankrow,6)';
            InfectDen(ii,1:9) = DenData(tankrow,7);
            % Predator densities for tank ii
            PredDen(ii,1) = DenData(tankrow(1),2);
        end

    % B. Find clonal frequencies for each replicate on weeks 2, 6, 9
        % row is tank, column is freq of each clone
        Freq_week0 = zeros(52,9)+NaN;
        Freq_week2 = zeros(52,9)+NaN; 
        Freq_week6 = zeros(52,9)+NaN; 
        Freq_week9 = zeros(52,9)+NaN; 

        for ii = [1:24 27:32 35:52]
            % Find rows for tank ii 
            tankrow = find(FreqData(:,3)==ii);
            % Host frequencies for tank ii
            Freq_temp = FreqData(tankrow,7:15);

            % Find rows for weeks 2, 6, 9
            % Note that weeks are recorded as -1, 1, 5, 8 in data set
            weekcol_entry1 = find(FreqData(tankrow,4)==1);
            weekcol_entry5 = find(FreqData(tankrow,4)==5);
            weekcol_entry8 = find(FreqData(tankrow,4)==8);

            % Set host frequencies for tank ii at week -1
            Freq_week0(ii,1:9) = zeros(1,9)+1/9;

            % Host frequencies for tank ii at weeks 1, 5, 8
            if isempty(weekcol_entry1) 
            else
                Freq_week2(ii,1:9) = Freq_temp(weekcol_entry1,:);
            end
            if isempty(weekcol_entry5) 
            else
                Freq_week6(ii,1:9) = Freq_temp(weekcol_entry5,:);
            end
            if isempty(weekcol_entry8) 
            else
                Freq_week9(ii,1:9) = Freq_temp(weekcol_entry8,:);
            end
        end

        
        
%% Step 5: Apply price equation 

    %A. Values to be reported 
        R_msred = zeros(52,4)+NaN;      % Estimated R values at each week
        DeltaEco_1a = zeros(52,3)+NaN;  % Effect of Ecology
        DeltaEvo_1a = zeros(52,3)+NaN;  % Effect of Evolution
   
    %B. Compute real and hypothetical values of R, then compute effect sizes
        % R_X_Y is R value computed using densities at time X and 
        % frequencies at time Y
    for ii = [1:24 27:32 35:52]
       
       % i. Set value of xi based on level of predation 
            if PredDen(ii) == 0
                xi_i = xi_0pred;
            elseif PredDen(ii) == 0.1
                xi_i = xi_01pred;
            elseif PredDen(ii) == 0.5
                xi_i = xi_05pred;
            elseif PredDen(ii) == 1
                xi_i = xi_1pred;
            else
                xi_i = Error;
            end
        
       % ii. Compute "real" R(N,q) values
            % Compute R(N(0),q(0))
            % NOTE: Predator density is zero for all replicates in week 0
            % Shared demoninator is m + sum(fSi*(N-I)*qi + fIi*I*qi)
            denom = delta + sum(...
               fSi.*(Freq_week0(ii,1:9)*(HostDen(ii,1)-InfectDen(ii,1)))+...
               fIi.*(Freq_week0(ii,1:9)*InfectDen(ii,1)));
            % Sum of numerators 
            numer = sum(p_i.*fSi.*(Freq_week0(ii,1:9)*...
               (HostDen(ii,1)-InfectDen(ii,1))).*chi_i...
               .*(mu_i+omega*xi_i.*a_i.*PredDen(ii))./(mu_i+omega*a_i*PredDen(ii)));
            R_0_0 = numer/denom;
            R_msred(ii,1) = R_0_0;
        
            % Compute R(N(2),q(2))
            % Shared demoninator is m + sum(fSi*(N-I)*qi + fIi*I*qi)
            denom = delta + sum(...
               fSi.*(Freq_week2(ii,1:9)*(HostDen(ii,2)-InfectDen(ii,2)))+...
               fIi.*(Freq_week2(ii,1:9)*InfectDen(ii,2)));
            % Sum of numerators 
            numer = sum(p_i.*fSi.*(Freq_week2(ii,1:9)*...
               (HostDen(ii,2)-InfectDen(ii,2))).*chi_i...
               .*(mu_i+omega*xi_i.*a_i.*PredDen(ii))./(mu_i+omega*a_i*PredDen(ii)));
            R_2_2 = numer/denom;
            R_msred(ii,2) = R_2_2;
       
            % Compute R(N(6),q(6))
            % Shared demoninator is m + sum(fSi*(N-I)*qi + fIi*I*qi)
            denom = delta + sum(...
               fSi.*(Freq_week6(ii,1:9)*(HostDen(ii,6)-InfectDen(ii,6)))+...
               fIi.*(Freq_week6(ii,1:9)*InfectDen(ii,6)));
            % Sum of numerators
            numer = sum(p_i.*fSi.*(Freq_week6(ii,1:9)*...
               (HostDen(ii,6)-InfectDen(ii,6))).*chi_i...
                .*(mu_i+omega*xi_i.*a_i.*PredDen(ii))./(mu_i+omega*a_i*PredDen(ii)));
            R_6_6 = numer/denom;
            R_msred(ii,3) = R_6_6;
        
            % Compute R(N(9),q(9))
            % Shared demoninator is m + sum(fSi*(N-I)*qi + fIi*I*qi)
            denom = delta + sum(...
               fSi.*(Freq_week9(ii,1:9)*(HostDen(ii,9)-InfectDen(ii,9)))+...
               fIi.*(Freq_week9(ii,1:9)*InfectDen(ii,9)));
            % Sum of mumerators
            numer = sum(p_i.*fSi.*(Freq_week9(ii,1:9)*...
               (HostDen(ii,9)-InfectDen(ii,9))).*chi_i...
                .*(mu_i+omega*xi_i.*a_i.*PredDen(ii))./(mu_i+omega*a_i*PredDen(ii)));
            R_9_9 = numer/denom;
            R_msred(ii,4) = R_9_9;
       
       % iii. Compute hypothetical R(N(t),q(t+1)) values
            % Compute R(N(0),q(2))
            % NOTE: Predator density is zero for all replicates in week -1
            % Shared demoninator is m + sum(fSi*(N-I)*qi + fIi*I*qi)
            denom = delta + sum(...
               fSi.*(Freq_week2(ii,1:9)*(HostDen(ii,1)-InfectDen(ii,1)))+...
               fIi.*(Freq_week2(ii,1:9)*InfectDen(ii,1)));
            % Sum of numerator
            numer = sum(p_i.*fSi.*(Freq_week2(ii,1:9)*...
               (HostDen(ii,1)-InfectDen(ii,1))).*chi_i...
               .*(mu_i+omega*xi_i.*a_i.*PredDen(ii))./(mu_i+omega*a_i*PredDen(ii)));
            R_0_2 = numer/denom;
       
            % Compute R(N(2),q(6))
            % Shared demoninator is m + sum(fSi*(N-I)*qi + fIi*I*qi)
            denom = delta + sum(...
               fSi.*(Freq_week6(ii,1:9)*(HostDen(ii,2)-InfectDen(ii,2)))+...
               fIi.*(Freq_week6(ii,1:9)*InfectDen(ii,2)));
            % Sum of numerators
            numer = sum(p_i.*fSi.*(Freq_week6(ii,1:9)*...
               (HostDen(ii,2)-InfectDen(ii,2))).*chi_i...
               .*(mu_i+omega*xi_i.*a_i.*PredDen(ii))./(mu_i+omega*a_i*PredDen(ii)));
            R_2_6 = numer/denom;
       
            % Compute R(N(6),q(9))
            % Shared demoninator is m + sum(fSi*(N-I)*qi + fIi*I*qi)
            denom = delta + sum(...
               fSi.*(Freq_week9(ii,1:9)*(HostDen(ii,6)-InfectDen(ii,6)))+...
               fIi.*(Freq_week9(ii,1:9)*InfectDen(ii,6)));
            % Sum of numerators
            numer = sum(p_i.*fSi.*(Freq_week9(ii,1:9)*...
               (HostDen(ii,6)-InfectDen(ii,6))).*chi_i...
               .*(mu_i+omega*xi_i.*a_i.*PredDen(ii))./(mu_i+omega*a_i*PredDen(ii)));
            R_6_9 = numer/denom;
       
       % iv. Compute hypothetical R(N(t+1),q(t)) values
            % Compute R(N(2),q(0))
            % Shared demoninator is m + sum(fSi*(N-I)*qi + fIi*I*qi)
            denom = delta + sum(...
               fSi.*(Freq_week0(ii,1:9)*(HostDen(ii,2)-InfectDen(ii,2)))+...
               fIi.*(Freq_week0(ii,1:9)*InfectDen(ii,2)));
            % Sum of numerators
            numer = sum(p_i.*fSi.*(Freq_week0(ii,1:9)*...
               (HostDen(ii,2)-InfectDen(ii,2))).*chi_i...
               .*(mu_i+omega*xi_i.*a_i.*PredDen(ii))./(mu_i+omega*a_i*PredDen(ii)));
            R_2_0 = numer/denom;
       
            % Compute R(N(6),q(2))
            % Shared demoninator is m + sum(fSi*(N-I)*qi + fIi*I*qi)
            denom = delta + sum(...
               fSi.*(Freq_week2(ii,1:9)*(HostDen(ii,6)-InfectDen(ii,6)))+...
               fIi.*(Freq_week2(ii,1:9)*InfectDen(ii,6)));
            % Sum of numerators
            numer = sum(p_i.*fSi.*(Freq_week2(ii,1:9)*...
               (HostDen(ii,6)-InfectDen(ii,6))).*chi_i...
               .*(mu_i+omega*xi_i.*a_i.*PredDen(ii))./(mu_i+omega*a_i*PredDen(ii)));
            R_6_2 = numer/denom;
       
            % Compute R(N(9),q(6))
            % Shared demoninator is m + sum(fSi*(N-I)*qi + fIi*I*qi)
            denom = delta + sum(...
               fSi.*(Freq_week6(ii,1:9)*(HostDen(ii,8)-InfectDen(ii,8)))+...
               fIi.*(Freq_week6(ii,1:9)*InfectDen(ii,8)));
            % Sum of numerators
            numer = sum(p_i.*fSi.*(Freq_week6(ii,1:9)*...
               (HostDen(ii,8)-InfectDen(ii,8))).*chi_i...
               .*(mu_i+omega*xi_i.*a_i.*PredDen(ii))./(mu_i+omega*a_i*PredDen(ii)));
            R_9_6 = numer/denom;
       
       % v. Compute effect sizes
            % Compute ecology effect sizes
            DeltaEco_1a(ii,1) = (R_2_2-R_0_2 + R_2_0-R_0_0)/2 ...
                /(SampleDay(ii,2) -SampleDay(ii,1));
            DeltaEco_1a(ii,2) = (R_6_6-R_2_6 + R_6_2-R_2_2)/2 ...
                /(SampleDay(ii,6) -SampleDay(ii,2));
            DeltaEco_1a(ii,3) = (R_9_9-R_6_9 + R_9_6-R_6_6)/2 ...
                /(SampleDay(ii,9) -SampleDay(ii,6));
   
            % Compute evolution effect sizes
            DeltaEvo_1a(ii,1) = (R_2_2-R_2_0 + R_0_2-R_0_0)/2 ...
                /(SampleDay(ii,2) -SampleDay(ii,1));
            DeltaEvo_1a(ii,2) = (R_6_6-R_6_2 + R_2_6-R_2_2)/2 ...
                /(SampleDay(ii,6) -SampleDay(ii,2));
            DeltaEvo_1a(ii,3) = (R_9_9-R_9_6 + R_6_9-R_6_6)/2 ...
                /(SampleDay(ii,9) -SampleDay(ii,6));
        
    end

    
    
%% Step 6: Plot R values
   
    figure('Name','Calculated R Values');

    % Yes Disease Treatments
        subplot(1,2,2)
        hold on
            for ii =[7, 8, 9, 44, 45, 46 ]
                h0=plot([0,2],log10(R_msred(ii,1:2)),'-o', 'Color', '#6baed6','LineWidth',1.5,...
                    'MarkerFaceColor', '#6baed6');
            end
            for ii =[1, 23, 24, 30, 31, 32 ]
                h01=plot([0,2],log10(R_msred(ii,1:2)),'-o', 'Color', '#3182bd','LineWidth',1.5,...
                    'MarkerFaceColor','#3182bd');
            end
            for ii =[16, 17, 18, 38, 39, 40 ]
                h05=plot([0,2],log10(R_msred(ii,1:2)),'-o', 'Color', '#08519c','LineWidth',1.5,...
                    'MarkerFaceColor','#08519c');
            end
            for ii =[10, 11, 12, 50, 51, 52]
                h1=plot([0,2],log10(R_msred(ii,1:2)),'-ok','LineWidth',1.5,...
                    'MarkerFaceColor','k');
            end
            plot([0 9],[0 0],'--r','Color', '#69369e','LineWidth',2)
            
            legend([h0 h01 h05 h1],'0 predator/L','0.1 predator/L','0.5 predator/L',...
                '1 predator/L', 'FontSize', 16,'location','northwest')
            
        set(gca,'fontsize',16)
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.04, 0.025], 'LineWidth', 1)
        rectangle('Position',[0,-1.2,2,3.7],'EdgeColor','k','LineWidth', 1 )
        ylim([-1.2 2.5]); yticks([-1 0 1 2])
        xlim([0 2]); xticks([0 1 2])
        title('Parasites', 'FontSize', 20)
        ylabel('log($R$)', 'FontSize', 20,'interpreter','Latex')
        xlabel('Time (week)', 'FontSize', 20)
        text(-0.2,2.4,'b','FontSize',24,'FontWeight','Bold')

    % No Disease Treatments
        subplot(1,2,1)
        hold on
            for ii =[13, 14, 15, 35, 36, 37 ]
                h0=plot([0,2],log10(R_msred(ii,1:2)),'-o', 'Color', '#6baed6','LineWidth',1.5,...
                    'MarkerFaceColor', '#6baed6');
            end
            for ii =[2, 3, 22, 47, 48, 49]
                h01=plot([0,2],log10(R_msred(ii,1:2)),'-o', 'Color', '#3182bd','LineWidth',1.5,...
                    'MarkerFaceColor','#3182bd');
            end
            for ii =[19, 20, 21, 27, 28, 29]
                h05=plot([0,2],log10(R_msred(ii,1:2)),'-o', 'Color', '#08519c','LineWidth',1.5,...
                    'MarkerFaceColor','#08519c');
            end
            for ii =[4, 5, 6, 41, 42, 43]
                h1=plot([0,2],log10(R_msred(ii,1:2)),'-ok','LineWidth',1.5,...
                    'MarkerFaceColor','k');
            end
            plot([0 9],[0 0],'--r','Color', '#69369e','LineWidth',2)
            
            legend([h0 h01 h05 h1],'0 predator/L','0.1 predator/L','0.5 predator/L',...
                '1 predator/L', 'FontSize', 16,'location','northwest')
        
        set(gca,'fontsize',16)
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.04, 0.025], 'LineWidth', 1)
        rectangle('Position',[0,-1.2,2,3.7],'EdgeColor','k','LineWidth', 1 )
        ylim([-1.2 2.5]); yticks([-1 0 1 2])
        xlim([0 2]); xticks([0 1 2])
        title('No parasites', 'FontSize', 20)
        ylabel('log$(R_0)$', 'FontSize', 20,'interpreter','Latex')
        xlabel('Time (week)', 'FontSize', 20)    
        text(-0.2,2.4,'a','FontSize',24,'FontWeight','Bold')

    

%% Step 7: Plot effect sizes 
%     figure('Name','Effect Sizes - Exact Frequencies');
% 
%     % Yes Disease Treatments
%         subplot(4,4,1)
%         hold on
%             for ii =[7, 8, 9, 44, 45, 46 ]
%                 plot([2,6,9],DeltaEco_1a(ii,:),'-ok')
%             end
%         title('Yes Disease, \Delta_{eco}')
%         ylabel({'{\bf 0 Pred Density}'; '\Delta_{eco}'})
%         ylim([-1,4])
% 
%         subplot(4,4,2)
%         hold on
%             for ii =[7, 8, 9, 44, 45, 46 ]
%                 plot([2,6,9],DeltaEvo_1a(ii,:),'-ok')
%             end
%         title('Yes Disease, \Delta_{evo}')
%         ylabel({''; '\Delta_{evo}'})
%         ylim([-1,4])
%         
%         subplot(4,4,5)
%         hold on
%             for ii =[1, 23, 24, 30, 31, 32 ]
%                 plot([2,6,9],DeltaEco_1a(ii,:),'-ok')
%             end
%         ylabel({'{\bf 0.1 Pred Density}'; '\Delta_{eco}'})
%         ylim([-1,3])
%         
%         subplot(4,4,6)
%         hold on
%             for ii =[1, 23, 24, 30, 31, 32 ]
%                 plot([2,6,9],DeltaEvo_1a(ii,:),'-ok')
%             end
%         ylabel({''; '\Delta_{evo}'})
%         ylim([-1,3])
%         
%         
%         subplot(4,4,9)
%         hold on
%             for ii =[16, 17, 18, 38, 39, 40 ]
%                 plot([2,6,9],DeltaEco_1a(ii,:),'-ok')
%             end
%         ylabel({'{\bf 0.5 Pred Density}'; '\Delta_{eco}'})
%         ylim([-1,1])
%         
%         subplot(4,4,10)
%         hold on
%             for ii =[16, 17, 18, 38, 39, 40 ]
%                 plot([2,6,9],DeltaEvo_1a(ii,:),'-ok')
%             end
%         ylabel({''; '\Delta_{evo}'})
%         ylim([-1,1])
% 
%         subplot(4,4,13)
%         hold on
%             for ii =[10, 11, 12, 50, 51, 52]
%                 plot([2,6,9],DeltaEco_1a(ii,:),'-ok')
%             end
%         ylabel({'{\bf 1 Pred Density}'; '\Delta_{eco}'})
%         ylim([-1.5,1])
%         xlabel('Week')
%         
%         subplot(4,4,14)
%         hold on
%             for ii =[10, 11, 12, 50, 51, 52]
%                 plot([2,6,9],DeltaEvo_1a(ii,:),'-ok')
%             end
%         ylabel({''; '\Delta_{evo}'})
%         ylim([-1.5,1])
%         xlabel('Week')
%         
% 
%     % No Disease Treatments
%         subplot(4,4,3)
%         hold on
%             for ii =[13, 14, 15, 35, 36, 37 ]
%                 plot([2,6,9],DeltaEco_1a(ii,:),'-ok')
%             end
%         title('No Disease, \Delta_{eco}')
%         ylabel({''; '\Delta_{eco}'})
%         ylim([-1,4])
%         
%         subplot(4,4,4)
%         hold on
%             for ii =[13, 14, 15, 35, 36, 37 ]
%                 plot([2,6,9],DeltaEvo_1a(ii,:),'-ok')
%             end
%         title('No Disease, \Delta_{evo}')
%         ylabel({''; '\Delta_{evo}'})
%         ylim([-1,4])
% 
%         subplot(4,4,7)
%         hold on
%             for ii =[2, 3, 22, 47, 48, 49]
%                 plot([2,6,9],DeltaEco_1a(ii,:),'-ok')
%             end
%         ylabel({''; '\Delta_{eco}'})
%         ylim([-1,3])
%         
%         subplot(4,4,8)
%         hold on
%             for ii =[2, 3, 22, 47, 48, 49]
%                 plot([2,6,9],DeltaEvo_1a(ii,:),'-ok')
%             end
%         ylabel({''; '\Delta_{evo}'})
%         ylim([-1,3])
% 
%         subplot(4,4,11)
%         hold on
%             for ii =[19, 20, 21, 27, 28, 29]
%                 plot([2,6,9],DeltaEco_1a(ii,:),'-ok')
%             end
%         ylabel({''; '\Delta_{eco}'})
%         ylim([-1,1])
%         
%         subplot(4,4,12)
%         hold on
%             for ii =[19, 20, 21, 27, 28, 29]
%                 plot([2,6,9],DeltaEvo_1a(ii,:),'-ok')
%             end
%         ylabel({''; '\Delta_{evo}'})
%         ylim([-1,1])
%         
% 
%         subplot(4,4,15)
%         hold on
%             for ii =[4, 5, 6, 41, 42, 43]
%                 plot([2,6,9],DeltaEco_1a(ii,:),'-ok')
%             end
%         ylabel({''; '\Delta_{eco}'})
%         ylim([-1.5,1])
%         xlabel('Week')
%         
%         subplot(4,4,16)
%         hold on
%             for ii =[4, 5, 6, 41, 42, 43]
%                 plot([2,6,9],DeltaEvo_1a(ii,:),'-ok')
%             end
%         ylabel({''; '\Delta_{evo}'})
%         ylim([-1.5,1])
%         xlabel('Week')
%         
             
        
%% Step 8: Plot ratios of effect sizes 
%     figure('Name','Ratios of Effect Sizes - Exact Frequencies');
% 
%     % Yes Disease Treatments
%         subplot(4,2,1)
%         hold on
%             for ii =[7, 8, 9, 44, 45, 46 ]
%                 plot([2,6,9],log10(abs(DeltaEco_1a(ii,:)./DeltaEvo_1a(ii,:))),'-ok')
%             end
%         title('Yes Disease, R')
%         plot([-1 5],[0 0],'-r')
%         ylabel({'{\bf 0 Pred Density}'; 'log10(\Delta_{eco}/\Delta_{evo})'})
%         ylim([-2,2])
% 
%         subplot(4,2,3)
%         hold on
%             for ii =[1, 23, 24, 30, 31, 32 ]
%                 plot([2,6,9],log10(abs(DeltaEco_1a(ii,:)./DeltaEvo_1a(ii,:))),'-ok')
%             end
%         ylabel({'{\bf 0.1 Pred Density}'; 'log10(\Delta_{eco}/\Delta_{evo})'})
%         plot([-1 5],[0 0],'-r')
%         ylim([-2,2])
%         
%         subplot(4,2,5)
%         hold on
%             for ii =[16, 17, 18, 38, 39, 40 ]
%                 plot([2,6,9],log10(abs(DeltaEco_1a(ii,:)./DeltaEvo_1a(ii,:))),'-ok')
%             end
%         ylabel({'{\bf 0.5 Pred Density}'; 'log10(\Delta_{eco}/\Delta_{evo})'})
%         plot([-1 5],[0 0],'-r')
%         ylim([-2,2])
%         
%         subplot(4,2,7)
%         hold on
%             for ii =[10, 11, 12, 50, 51, 52]
%                 plot([2,6,9],log10(abs(DeltaEco_1a(ii,:)./DeltaEvo_1a(ii,:))),'-ok')
%             end
%         ylabel({'{\bf 1 Pred Density}'; 'log10(\Delta_{eco}/\Delta_{evo})'})
%         plot([-1 5],[0 0],'-r')
%         ylim([-2,2])
%         xlabel('Week')
%         
%        
%     % No Disease Treatments
%         subplot(4,2,2)
%         hold on
%             for ii =[13, 14, 15, 35, 36, 37 ]
%                 plot([2,6,9],log10(abs(DeltaEco_1a(ii,:)./DeltaEvo_1a(ii,:))),'-ok')
%             end
%         title('No Disease, R0')
%         ylabel({''; 'log10(\Delta_{eco}/\Delta_{evo})'})
%         plot([-1 5],[0 0],'-r')
%         ylim([-2,2])
%         
%        
%         subplot(4,2,4)
%         hold on
%             for ii =[2, 3, 22, 47, 48, 49]
%                 plot([2,6,9],log10(abs(DeltaEco_1a(ii,:)./DeltaEvo_1a(ii,:))),'-ok')
%             end
%         ylabel({''; 'log10(\Delta_{eco}/\Delta_{evo})'})
%         plot([-1 5],[0 0],'-r')
%         ylim([-2,2])
%         
%         subplot(4,2,6)
%         hold on
%             for ii =[19, 20, 21, 27, 28, 29]
%                 plot([2,6,9],log10(abs(DeltaEco_1a(ii,:)./DeltaEvo_1a(ii,:))),'-ok')
%             end
%         ylabel({''; 'log10(\Delta_{eco}/\Delta_{evo})'})
%         plot([-1 5],[0 0],'-r')
%         ylim([-2,2])
%         
%         subplot(4,2,8)
%         hold on
%             for ii =[4, 5, 6, 41, 42, 43]
%                 plot([2,6,9],log10(abs(DeltaEco_1a(ii,:)./DeltaEvo_1a(ii,:))),'-ok')
%             end
%         ylabel({''; 'log10(\Delta_{eco}/\Delta_{evo})'})
%         plot([-1 5],[0 0],'-r')
%         ylim([-2,2])
%         xlabel('Week')
        
              

%%%% Plot combining panels with ecology, evolution, and ratio

% Single plot 
figure('Name','Geber Method Results');

    % Ecology
    % Yes Disease Treatments
        subplot(3,2,2)
        hold on
            for ii =[7, 8, 9, 44, 45, 46 ]
                h0=plot([0],DeltaEco_1a(ii,1),'-o', 'Color', '#6baed6','LineWidth',1.5,...
                    'MarkerFaceColor','#6baed6');
            end
            for ii =[1, 23, 24, 30, 31, 32 ]
                h01=plot([0.1],DeltaEco_1a(ii,1),'-o', 'Color', '#3182bd','LineWidth',1.5,...
                    'MarkerFaceColor','#3182bd');
            end
            for ii =[16, 17, 18, 38, 39, 40 ]
                h05=plot([0.5],DeltaEco_1a(ii,1),'-o', 'Color', '#08519c','LineWidth',1.5,...
                    'MarkerFaceColor','#08519c');
            end
            for ii =[10, 11, 12, 50, 51, 52]
                h1=plot([1],DeltaEco_1a(ii,1),'-ok','LineWidth',1.5,...
                    'MarkerFaceColor','k');
            end
            
            legend([h0 h01 h05 h1],'0 predator/L','0.1 predator/L','0.5 predator/L',...
                '1 predator/L', 'FontSize', 12)
            
        set(gca,'fontsize',16)
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.04, 0.025], 'LineWidth', 1)
        rectangle('Position',[0,-1.5,1,5.5],'EdgeColor','k','LineWidth', 1 )
        ylim([-1.5 4]); yticks([-1 0 1 2 3])
        xlim([0 1]); xticks([0 0.1 0.5 1])
        ylabel('$\Delta_{eco}$', 'Fontsize', 20,'interpreter','Latex','color','#FFFFFF')
        xlabel('Predator density', 'Fontsize', 20,'color','#FFFFFF')
        text(-0.2,4,'b','FontSize',24,'FontWeight','Bold')
        text(0.2,5.5,'Parasites{\it (R)}','Fontsize',24,'FontWeight','Bold');
        
        
    % Ecology
    % No Disease Treatments
        subplot(3,2,1)
        hold on
            for ii =[13, 14, 15, 35, 36, 37 ]
                h0=plot([0],DeltaEco_1a(ii,1),'-o', 'Color', '#6baed6','LineWidth',1.5,...
                    'MarkerFaceColor','#6baed6');
            end
            for ii =[2, 3, 22, 47, 48, 49]
                h01=plot([0.1],DeltaEco_1a(ii,1),'-o', 'Color', '#3182bd','LineWidth',1.5,...
                    'MarkerFaceColor','#3182bd');
            end
            for ii =[19, 20, 21, 27, 28, 29]
                h05=plot([0.5],DeltaEco_1a(ii,1),'-o', 'Color', '#08519c','LineWidth',1.5,...
                    'MarkerFaceColor','#08519c');
            end
            for ii =[4, 5, 6, 41, 42, 43]
                h1=plot([1],DeltaEco_1a(ii,1),'-ok','LineWidth',1.5,...
                    'MarkerFaceColor','k');
            end
            
            legend([h0 h01 h05 h1],'0 predator/L','0.1 predator/L','0.5 predator/L',...
                '1 predator/L', 'FontSize', 12)
            
        set(gca,'fontsize',16)
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.04, 0.025], 'LineWidth', 1)
        rectangle('Position',[0,-1.5,1,5.5],'EdgeColor','k','LineWidth', 1 )
        ylim([-1.5 4]); yticks([-1 0 1 2 3])
        xlim([0 1]); xticks([0 0.1 0.5 1])
        ylabel('$\Delta_{eco}$', 'Fontsize', 20,'interpreter','Latex')
        xlabel('Predator density', 'Fontsize', 20,'color','#FFFFFF')  
        text(-0.2,4,'a','FontSize',24,'FontWeight','Bold')
        text(-0.3,-1,'Ecology','Fontsize',24,'FontWeight','Bold','rotation',90);
        text(0.1,5.5,'No parasites{\it (R_0)}','Fontsize',24,'FontWeight','Bold');

    % Evolution
    % Yes Disease Treatments
        subplot(3,2,4)
        hold on
            for ii =[7, 8, 9, 44, 45, 46 ]
                h0=plot([0],DeltaEvo_1a(ii,1),'-o', 'Color', '#6baed6','LineWidth',1.5,...
                    'MarkerFaceColor','#6baed6');
            end
            for ii =[1, 23, 24, 30, 31, 32 ]
                h01=plot([0.1],DeltaEvo_1a(ii,1),'-o', 'Color', '#3182bd','LineWidth',1.5,...
                    'MarkerFaceColor','#3182bd');
            end
            for ii =[16, 17, 18, 38, 39, 40 ]
                h05=plot([0.5],DeltaEvo_1a(ii,1),'-o', 'Color', '#08519c','LineWidth',1.5,...
                    'MarkerFaceColor','#08519c');
            end
            for ii =[10, 11, 12, 50, 51, 52]
                h1=plot([1],DeltaEvo_1a(ii,1),'-ok','LineWidth',1.5,...
                    'MarkerFaceColor','k');
            end
            
  %          legend([h0 h01 h05 h1],'0 pred','0.1 pred','0.5 pred','1 pred')
            
        set(gca,'fontsize',16)
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.04, 0.025], 'LineWidth', 1)
        rectangle('Position',[0,-1,1,2.5],'EdgeColor','k','LineWidth', 1 )
        ylim([-1 1.5]); yticks([-1 0 1])
        xlim([0 1]); xticks([0 0.1 0.5 1])
        ylabel('$\Delta_{evo}$', 'Fontsize', 20,'interpreter','Latex','color','#FFFFFF')
        xlabel('Predator density', 'Fontsize', 20,'color','#FFFFFF')    
        text(-0.2,1.5,'d','FontSize',24,'FontWeight','Bold')

    % Evolution
    % No Disease Treatments
        subplot(3,2,3)
        hold on
            for ii =[13, 14, 15, 35, 36, 37 ]
                h0=plot([0],DeltaEvo_1a(ii,1),'-o', 'Color', '#6baed6','LineWidth',1.5,...
                    'MarkerFaceColor','#6baed6');
            end
            for ii =[2, 3, 22, 47, 48, 49]
                h01=plot([0.1],DeltaEvo_1a(ii,1),'-o', 'Color', '#3182bd','LineWidth',1.5,...
                    'MarkerFaceColor','#3182bd');
            end
            for ii =[19, 20, 21, 27, 28, 29]
                h05=plot([0.5],DeltaEvo_1a(ii,1),'-o', 'Color', '#08519c','LineWidth',1.5,...
                    'MarkerFaceColor','#08519c');
            end
            for ii =[4, 5, 6, 41, 42, 43]
                h1=plot([1],DeltaEvo_1a(ii,1),'-ok','LineWidth',1.5,...
                    'MarkerFaceColor','k');
            end
            
  %          legend([h0 h01 h05 h1],'0 pred','0.1 pred','0.5 pred','1 pred')
            
        set(gca,'fontsize',16)
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.04, 0.025], 'LineWidth', 1)
        rectangle('Position',[0,-1,1,2.5],'EdgeColor','k','LineWidth', 1 )
        ylim([-1 1.5]); yticks([-1 0 1])
        xlim([0 1]); xticks([0 0.1 0.5 1])
        ylabel('$\Delta_{evo}$', 'Fontsize', 20,'interpreter','Latex')
        xlabel('Predator density', 'Fontsize', 20,'color','#FFFFFF') 
        text(-0.3,-1,'Evolution','Fontsize',24,'FontWeight','Bold','rotation',90);
        text(-0.2,1.5,'c','FontSize',24,'FontWeight','Bold')

    % Ratio of eco:evo
    % Yes Disease Treatments
        subplot(3,2,6)
        hold on
            for ii =[7, 8, 9, 44, 45, 46 ]
                h0=plot([0],log10(abs(DeltaEco_1a(ii,1)./DeltaEvo_1a(ii,1))),'-o', 'Color', '#6baed6','LineWidth',1.5,...
                    'MarkerFaceColor','#6baed6');
            end
            for ii =[1, 23, 24, 30, 31, 32 ]
                h01=plot([0.1],log10(abs(DeltaEco_1a(ii,1)./DeltaEvo_1a(ii,1))),'-o', 'Color', '#3182bd','LineWidth',1.5,...
                    'MarkerFaceColor','#3182bd');
            end
            for ii =[16, 17, 18, 38, 39, 40 ]
                h05=plot([0.5],log10(abs(DeltaEco_1a(ii,1)./DeltaEvo_1a(ii,1))),'-o', 'Color', '#08519c','LineWidth',1.5,...
                    'MarkerFaceColor','#08519c');
            end
            for ii =[10, 11, 12, 50, 51, 52]
                h1=plot([1],log10(abs(DeltaEco_1a(ii,1)./DeltaEvo_1a(ii,1))),'-ok','LineWidth',1.5,...
                    'MarkerFaceColor','k');
            end
            plot([0 8],[0 0],'--r','Color', '#69369e','LineWidth',2)
            
    %        legend([h0 h01 h05 h1],'0 pred','0.1 pred','0.5 pred','1 pred')
            
        set(gca,'fontsize',16)
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.04, 0.025], 'LineWidth', 1)
        rectangle('Position',[0,-1.5,1,4.5],'EdgeColor','k','LineWidth', 1 )
        ylim([-1.5 3]); yticks([-1 0 1 2 3])
        xlim([0 1]); xticks([0 0.1 0.5 1])
        %title('Eco:Evo', 'Fontsize', 20)
        ylabel('log$(|\Delta_{eco}/\Delta_{evo}|)$', 'Fontsize', 20,'interpreter','Latex','color','#FFFFFF')
        xlabel('Predator density', 'Fontsize', 20)    
        text(-0.2,4,'f','FontSize',24,'FontWeight','Bold')
        
    % Ratio of eco:evo
    % No Disease Treatments
        subplot(3,2,5)
        hold on
            for ii =[13, 14, 15, 35, 36, 37 ]
                h0=plot([0],log10(abs(DeltaEco_1a(ii,1)./DeltaEvo_1a(ii,1))),'-o', 'Color', '#6baed6','LineWidth',1.5,...
                    'MarkerFaceColor','#6baed6');
            end
            for ii =[2, 3, 22, 47, 48, 49]
                h01=plot([0.1],log10(abs(DeltaEco_1a(ii,1)./DeltaEvo_1a(ii,1))),'-o', 'Color', '#3182bd','LineWidth',1.5,...
                    'MarkerFaceColor','#3182bd');
            end
            for ii =[19, 20, 21, 27, 28, 29]
                h05=plot([0.5],log10(abs(DeltaEco_1a(ii,1)./DeltaEvo_1a(ii,1))),'-o', 'Color', '#08519c','LineWidth',1.5,...
                    'MarkerFaceColor','#08519c');
            end
            for ii =[4, 5, 6, 41, 42, 43]
                h1=plot([1],log10(abs(DeltaEco_1a(ii,1)./DeltaEvo_1a(ii,1))),'-ok','LineWidth',1.5,...
                    'MarkerFaceColor','k');
            end
            plot([0 8],[0 0],'--r','Color', '#69369e','LineWidth',2)
            
        %    legend([h0 h01 h05 h1],'0 pred','0.1 pred','0.5 pred','1 pred')
            
        set(gca,'fontsize',16)
        set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.04, 0.025], 'LineWidth', 1)
        rectangle('Position',[0,-1.5,1,4.5],'EdgeColor','k','LineWidth', 1 )
        ylim([-1.5 3]); yticks([-1 0 1 2 3])
        xlim([0 1]); xticks([0 0.1 0.5 1])
        %title('Eco:evo', 'Fontsize', 20)
        ylabel('log$(|\Delta_{eco}/\Delta_{evo}|)\:\:\:\:\:\:\:\:\:\:$', 'Fontsize', 20,'interpreter','Latex')
        xlabel('Predator density', 'Fontsize', 20)    
        text(-0.3,-2,'eco:evo','Fontsize',24,'FontWeight','Bold','rotation',90);
        text(-0.2,4,'e','FontSize',24,'FontWeight','Bold')
    
    
   
