% Script for evaluating equations related to healthy herd hypothesis

% AUTHOR:   Michael Cortez, Florida State University, mcortez@fsu.edu
% DATE:     March 2022
% PURPOSE:  Evaluate sensitivity equations defining how increased predator
%           density affects total host density

% INPUTS:   genotype.frequencies_MHC_converted.csv -- data table of clonal
%           frequencies in weeks 0, 2, 6, and 9
%           meso.expt.densities_MHC_converted.csv -- data table of total
%           host densities and infected densities in weeks 1-9.
%           Order of clones is (1) BD05-42, (2) BD08-46	(3) BD19-64	
%           (4) CB24-68	(5) DW22-58	(6) DW29-75	(7) IL14-43	(8) ML30-82
%           (9) ML32-84

% METHOD:   Step 1: Define parameters estimated from other experiments 
%           Step 2: Import data
%           Step 3: Extract densities and frequencies for weeks 0,2,6,9
%           Step 4: Compute sign of sensitivity

% OUTPUTS:  Signs of \partial N^*/\partial P (i.e., negative of numerator) when
%           (i) all clones present at equal frequencies
%           (ii) clone frequencies at week 2 of "parasites" treatments
%           (iii) clone frequencies at week 2 of "no parasites" treatments

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
    omega = 2;
    
    % spore yield (spores/indiv.) 
    chi_i = [78231 36990 83035 70578 89909 146128 108163 41667 92610];
    
    % Disease-induced morality rate (1/hr)
    mu_i = 1/20/24;
    
    % spore degradation rate (1/day)
    delta = 0.0083;
    
    % R0_i(infinity) Competence for host
    R0i = p_i.*chi_i;
    
    % initial density of spores (spores)
    z0_i = [21   57  128  109  139  226  167   64  143];
    
    % maximum burst size (spores)
    sigma_i = 1.0e+05*[0.8059 0.3811 0.8554 0.7271 0.9262 1.5054 1.1143 0.4293 0.9541];
    
    % maximum exponential growth rate of spores within individuals (1/hr)
    % Taken from Auld et al. 2014
    r = 0.0208;
    
    % Reduction in spore burst size for P = 0, 0.1, 0.5, and 1
    x_i_0pred = [1 1 1 1 1 1 1 1 1];
    x_i_01pred = [0.2104 0.4099 0.4954 0.7929 0.3681 0.6828 0.3098 0.3200 0.6032];
    x_i_05pred = [0.0088 0.0163 0.0213 0.0796 0.0144 0.0431 0.0121 0.0125 0.0310];
    x_i_1pred=[0.0040 0.0059 0.0071 0.0188 0.0055 0.0116 0.0049 0.0050 0.0091];

    % Derivative dx_i/dP at P = 0, 0.1, 0.5, and 1
    P=0;
    dx_idP_0pred = -x_i_0pred.*r.*omega.*a_i.*(sigma_i-z0_i)./(sigma_i-z0_i+...
        z0_i.*exp(r./(omega.*a_i.*P+mu_i)))/(P.*omega.*a_i+mu_i).^2;
    P=0.1;
    dx_idP_01pred = -x_i_01pred.*r.*omega.*a_i.*(sigma_i-z0_i)./(sigma_i-z0_i+...
        z0_i.*exp(r./(omega.*a_i.*P+mu_i)))/(P.*omega.*a_i+mu_i).^2;
    P=0.5;
    dx_idP_05pred = -x_i_05pred.*r.*omega.*a_i.*(sigma_i-z0_i)./(sigma_i-z0_i+...
        z0_i.*exp(r./(omega.*a_i.*P+mu_i)))/(P.*omega.*a_i+mu_i).^2;
    P=1;
    dx_idP_1pred = -x_i_1pred.*r.*omega.*a_i.*(sigma_i-z0_i)./(sigma_i-z0_i+...
        z0_i.*exp(r./(omega.*a_i.*P+mu_i)))/(P.*omega.*a_i+mu_i).^2;

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
        

%% Step 3: Extract Host densities and frequencies at weeks 0, 2, 6, 9

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

            % Find rows for weeks 1, 5, 8
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

        
        
%% Step 4: Compute sensitivities

    Yvals= 0:0.1:1;
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 scrsz(4) scrsz(3)*0.4 scrsz(4)*.3]);
    
% Case 1: All genotypes present with equal frequencies
    % Definte average parameters
    subplot(1,3,1)
    hold on
    xlabel('Equilibrium Prevalence (Y)','Fontsize', 14)
    ylabel('Sign of dN/dP', 'Fontsize', 14)
    %xticks([0 0.5 1])
    xlim([0 1])
    yticks([-0.03 -0.02 -0.01 0 ])
    ylim([-0.035 0.005])
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.04, 0.025], 'LineWidth', 1,'fontsize',14)
        rectangle('Position',[0,-0.035,1,0.04],'EdgeColor','k','LineWidth', 1 )
    title({'Equal frequencies';''}, 'Fontsize', 14)
    plot([0 1],[0 0],'--r','Color', '#69369e','LineWidth',1.5)
    
    Freqs = [1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9 1/9];
    a_S_bar = sum(Freqs.*a_i);
    a_I_bar = sum(Freqs.*a_i.*omega);
    chimu_bar = sum(Freqs.*chi_i.*mu_i);
    chim_bar = 0; % m_i assumed to be negligible
    chixai_bar_0pred = sum(Freqs.*chi_i.*x_i_0pred.*a_i.*omega);
    chixai_bar_01pred = sum(Freqs.*chi_i.*x_i_01pred.*a_i.*omega);
    chixai_bar_05pred = sum(Freqs.*chi_i.*x_i_05pred.*a_i.*omega);
    chixai_bar_1pred = sum(Freqs.*chi_i.*x_i_1pred.*a_i.*omega);
    m_bar = 0; %m_i assumed to be negligible
    mu_bar = sum(Freqs.*mu_i);
    pf_bar = sum(Freqs.*p_i.*fSi); % S and I assumed to have equal uptake
    dxdP_0pred = sum(Freqs.*dx_idP_0pred); % derivative dx/dP at P=0
    dxdP_01pred = sum(Freqs.*dx_idP_01pred); % derivative dx/dP at P=0
    dxdP_05pred = sum(Freqs.*dx_idP_05pred); % derivative dx/dP at P=0
    dxdP_1pred = sum(Freqs.*dx_idP_1pred); % derivative dx/dP at P=0
    
    
    P=0;
    Z = (mu_bar+m_bar+a_I_bar*P).*Yvals./pf_bar./(1-Yvals);
    factorZ=(mu_bar-a_S_bar*P+a_I_bar*P)./(mu_bar+m_bar+a_I_bar*P+pf_bar*Z);
    dxdP = dxdP_0pred;
    dNdP_0pred = -(a_S_bar.*(1-Yvals)+a_I_bar.*Yvals + ...
        -a_I_bar.*(chimu_bar+chim_bar)./(chimu_bar+chim_bar+chixai_bar_0pred.*P)*factorZ+...
        chixai_bar_0pred.*(mu_bar+m_bar)./(chimu_bar+chim_bar+chixai_bar_0pred.*P)*factorZ+...
        (mu_bar+m_bar+a_I_bar*P)./(chimu_bar+chim_bar+chixai_bar_0pred.*P).*factorZ.*dxdP*P);
    
    P=0.1;
    Z = (mu_bar+m_bar+a_I_bar*P).*Yvals./pf_bar./(1-Yvals);
    dxdP = dxdP_01pred;
    factorZ=(mu_bar-a_S_bar*P+a_I_bar*P)./(mu_bar+m_bar+a_I_bar*P+pf_bar*Z);
    dNdP_01pred = -(a_S_bar.*(1-Yvals)+a_I_bar.*Yvals + ...
        -a_I_bar.*(chimu_bar+chim_bar)./(chimu_bar+chim_bar+chixai_bar_01pred.*P)*factorZ+...
        chixai_bar_01pred.*(mu_bar+m_bar)./(chimu_bar+chim_bar+chixai_bar_01pred.*P)*factorZ+...
        (mu_bar+m_bar+a_I_bar*P)./(chimu_bar+chim_bar+chixai_bar_01pred.*P).*factorZ.*dxdP*P);
    
    P=0.5;
    Z = (mu_bar+m_bar+a_I_bar*P).*Yvals./pf_bar./(1-Yvals);
    factorZ=(mu_bar-a_S_bar*P+a_I_bar*P)./(mu_bar+m_bar+a_I_bar*P+pf_bar*Z);
    dxdP = dxdP_05pred;
    dNdP_05pred = -(a_S_bar.*(1-Yvals)+a_I_bar.*Yvals + ...
        -a_I_bar.*(chimu_bar+chim_bar)./(chimu_bar+chim_bar+chixai_bar_05pred.*P)*factorZ+...
        chixai_bar_05pred.*(mu_bar+m_bar)./(chimu_bar+chim_bar+chixai_bar_05pred.*P)*factorZ+...
        (mu_bar+m_bar+a_I_bar*P)./(chimu_bar+chim_bar+chixai_bar_05pred.*P).*factorZ.*dxdP*P);
    P=1;
    Z = (mu_bar+m_bar+a_I_bar*P).*Yvals./pf_bar./(1-Yvals);
    factorZ=(mu_bar-a_S_bar*P+a_I_bar*P)./(mu_bar+m_bar+a_I_bar*P+pf_bar*Z);
    dxdP = dxdP_1pred;
    dNdP_1pred = -(a_S_bar.*(1-Yvals)+a_I_bar.*Yvals + ...
        -a_I_bar.*(chimu_bar+chim_bar)./(chimu_bar+chim_bar+chixai_bar_1pred.*P)*factorZ+...
        chixai_bar_1pred.*(mu_bar+m_bar)./(chimu_bar+chim_bar+chixai_bar_1pred.*P)*factorZ+...
        (mu_bar+m_bar+a_I_bar*P)./(chimu_bar+chim_bar+chixai_bar_1pred.*P).*factorZ.*dxdP*P);

    h0=plot(Yvals,dNdP_0pred,'-','Color', '#6baed6','LineWidth',1.5);
    h01=plot(Yvals,dNdP_01pred,'-','Color', '#3182bd','LineWidth',1.5);
    h05=plot(Yvals,dNdP_05pred,'-','Color', '#08519c','LineWidth',1.5);
    h1=plot(Yvals,dNdP_1pred,'-','Color', 'k','LineWidth',1.5);
    
    legend([h0 h01 h05 h1],'0 predator/L','0.1 predator/L','0.5 predator/L','1 predator/L',...
        'FontSize', 12,'location','southwest')
    text(-0.2,0.005,'a','FontSize',16,'FontWeight','Bold')
    
    
% Case 2: Week 2 frequencies for With parasites treatments
    subplot(1,3,3)
    hold on
    xlabel('Equilibrium Prevalence (Y)', 'Fontsize', 14)
    ylabel('Sign of dN/dP','FontWeight','bold','color','#FFFFFF')
    %xticks([0 0.5 1])
    xlim([0 1])
    yticks([-0.03 -0.02 -0.01 0 ])
    ylim([-0.035 0.005])
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.04, 0.025], 'LineWidth', 1,'fontsize',14)
        rectangle('Position',[0,-0.035,1,0.04],'EdgeColor','k','LineWidth', 1 )
    title({'Week 2 frequencies for'; 'parasite treatments'}, 'Fontsize', 14) 
    plot([0 1],[0 0],'--r','Color', '#69369e','LineWidth',1.5)
    
    % Replicates with predator density = 0
    for ii = [7, 8, 9, 44, 45, 46 ]
        
        P=0;
        Freqs = Freq_week2(ii,1:9);
        a_S_bar = sum(Freqs.*a_i);
        a_I_bar = sum(Freqs.*a_i.*omega);
        chimu_bar = sum(Freqs.*chi_i.*mu_i);
        chim_bar = 0; % m_i assumed to be negligible
        chixai_bar = sum(Freqs.*chi_i.*x_i_0pred.*a_i.*omega);
        m_bar = 0; %m_i assumed to be negligible
        mu_bar = sum(Freqs.*mu_i);
        pf_bar = sum(Freqs.*p_i.*fSi); % S and I assumed to have equal uptake
        dxdP = sum(Freqs.*dx_idP_0pred);
        
        Z = (mu_bar+m_bar+a_I_bar*P).*Yvals./pf_bar./(1-Yvals);
        factorZ=(mu_bar-a_S_bar*P+a_I_bar*P)./(mu_bar+m_bar+a_I_bar*P+pf_bar*Z);
        dNdP = -(a_S_bar.*(1-Yvals)+a_I_bar.*Yvals + ...
            -a_I_bar.*(chimu_bar+chim_bar)./(chimu_bar+chim_bar+chixai_bar.*P)*factorZ+...
            chixai_bar.*(mu_bar+m_bar)./(chimu_bar+chim_bar+chixai_bar.*P)*factorZ+...
            (mu_bar+m_bar+a_I_bar*P)./(chimu_bar+chim_bar+chixai_bar.*P).*factorZ.*dxdP*P);
        
        plot(Yvals,dNdP,'-','Color', '#6baed6','LineWidth',1.5)
    end

    % Replicates with predator density = 0.1
    for ii = [1, 23, 24, 30, 31, 32 ]
        
        P=0.1;
        Freqs = Freq_week2(ii,1:9);
        a_S_bar = sum(Freqs.*a_i);
        a_I_bar = sum(Freqs.*a_i.*omega);
        chimu_bar = sum(Freqs.*chi_i.*mu_i);
        chim_bar = 0; % m_i assumed to be negligible
        chixai_bar = sum(Freqs.*chi_i.*x_i_01pred.*a_i.*omega);
        m_bar = 0; %m_i assumed to be negligible
        mu_bar = sum(Freqs.*mu_i);
        pf_bar = sum(Freqs.*p_i.*fSi); % S and I assumed to have equal uptake
        dxdP = sum(Freqs.*dx_idP_01pred);
        
        Z = (mu_bar+m_bar+a_I_bar*P).*Yvals./pf_bar./(1-Yvals);
        factorZ=(mu_bar-a_S_bar*P+a_I_bar*P)./(mu_bar+m_bar+a_I_bar*P+pf_bar*Z);
        dNdP = -(a_S_bar.*(1-Yvals)+a_I_bar.*Yvals + ...
            -a_I_bar.*(chimu_bar+chim_bar)./(chimu_bar+chim_bar+chixai_bar.*P)*factorZ+...
            chixai_bar.*(mu_bar+m_bar)./(chimu_bar+chim_bar+chixai_bar.*P)*factorZ+...
            (mu_bar+m_bar+a_I_bar*P)./(chimu_bar+chim_bar+chixai_bar.*P).*factorZ.*dxdP*P);
        
        plot(Yvals,dNdP,'-','Color', '#3182bd','LineWidth',1.5)
    end


    % Replicates with predator density = 0.5
    for ii = [16, 17, 18, 38, 39, 40 ]
        
        P=0.5;
        Freqs = Freq_week2(ii,1:9);
        a_S_bar = sum(Freqs.*a_i);
        a_I_bar = sum(Freqs.*a_i.*omega);
        chimu_bar = sum(Freqs.*chi_i.*mu_i);
        chim_bar = 0; % m_i assumed to be negligible
        chixai_bar = sum(Freqs.*chi_i.*x_i_05pred.*a_i.*omega);
        m_bar = 0; %m_i assumed to be negligible
        mu_bar = sum(Freqs.*mu_i);
        pf_bar = sum(Freqs.*p_i.*fSi); % S and I assumed to have equal uptake
        dxdP = sum(Freqs.*dx_idP_05pred);
        
        Z = (mu_bar+m_bar+a_I_bar*P).*Yvals./pf_bar./(1-Yvals);
        factorZ=(mu_bar-a_S_bar*P+a_I_bar*P)./(mu_bar+m_bar+a_I_bar*P+pf_bar*Z);
        dNdP = -(a_S_bar.*(1-Yvals)+a_I_bar.*Yvals + ...
            -a_I_bar.*(chimu_bar+chim_bar)./(chimu_bar+chim_bar+chixai_bar.*P)*factorZ+...
            chixai_bar.*(mu_bar+m_bar)./(chimu_bar+chim_bar+chixai_bar.*P)*factorZ+...
            (mu_bar+m_bar+a_I_bar*P)./(chimu_bar+chim_bar+chixai_bar.*P).*factorZ.*dxdP*P);
        
        plot(Yvals,dNdP,'-','Color', '#08519c','LineWidth',1.5)
    end
    
    % Replicates with predator density = 1
    for ii =[10, 11, 12, 50, 51, 52]
        
        P=1;
        Freqs = Freq_week2(ii,1:9);
        a_S_bar = sum(Freqs.*a_i);
        a_I_bar = sum(Freqs.*a_i.*omega);
        chimu_bar = sum(Freqs.*chi_i.*mu_i);
        chim_bar = 0; % m_i assumed to be negligible
        chixai_bar = sum(Freqs.*chi_i.*x_i_1pred.*a_i.*omega);
        m_bar = 0; %m_i assumed to be negligible
        mu_bar = sum(Freqs.*mu_i);
        pf_bar = sum(Freqs.*p_i.*fSi); % S and I assumed to have equal uptake
        dxdP = sum(Freqs.*dx_idP_1pred);
        
        Z = (mu_bar+m_bar+a_I_bar*P).*Yvals./pf_bar./(1-Yvals);
        factorZ=(mu_bar-a_S_bar*P+a_I_bar*P)./(mu_bar+m_bar+a_I_bar*P+pf_bar*Z);
        dNdP = -(a_S_bar.*(1-Yvals)+a_I_bar.*Yvals + ...
            -a_I_bar.*(chimu_bar+chim_bar)./(chimu_bar+chim_bar+chixai_bar.*P)*factorZ+...
            chixai_bar.*(mu_bar+m_bar)./(chimu_bar+chim_bar+chixai_bar.*P)*factorZ+...
            (mu_bar+m_bar+a_I_bar*P)./(chimu_bar+chim_bar+chixai_bar.*P).*factorZ.*dxdP*P);
        
        plot(Yvals,dNdP,'-','Color', 'k','LineWidth',1.5)
    end

    text(-0.2,0.005,'c','FontSize',16,'FontWeight','Bold')
    
    
% Case 3: Final frequencies for No parasites treatments
    subplot(1,3,2)
    hold on
    xlabel('Equilibrium Prevalence (Y)', 'Fontsize', 14)
    ylabel('Sign of dN/dP','FontWeight','bold','color','#FFFFFF')
    %xticks([0 0.5 1])
    xlim([0 1])
    yticks([-0.03 -0.02 -0.01 0 ])
    ylim([-0.035 0.005])
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.04, 0.025], 'LineWidth', 1,'fontsize',14)
        rectangle('Position',[0,-0.035,1,0.04],'EdgeColor','k','LineWidth', 1 )
    title({'Week 2 frequencies for'; 'no parasite treatments' }, 'Fontsize', 14) 
    plot([0 1],[0 0],'--r','Color', '#69369e','LineWidth',1.5)
    
    % Replicates with predator density = 0
    for ii = [13, 14, 15, 35, 36, 37 ]
        
        P=0;
        Freqs = Freq_week2(ii,1:9);
        a_S_bar = sum(Freqs.*a_i);
        a_I_bar = sum(Freqs.*a_i.*omega);
        chimu_bar = sum(Freqs.*chi_i.*mu_i);
        chim_bar = 0; % m_i assumed to be negligible
        chixai_bar = sum(Freqs.*chi_i.*x_i_0pred.*a_i.*omega);
        m_bar = 0; %m_i assumed to be negligible
        mu_bar = sum(Freqs.*mu_i);
        pf_bar = sum(Freqs.*p_i.*fSi); % S and I assumed to have equal uptake
        dxdP = sum(Freqs.*dx_idP_0pred);
        
        Z = (mu_bar+m_bar+a_I_bar*P).*Yvals./pf_bar./(1-Yvals);
        factorZ=(mu_bar-a_S_bar*P+a_I_bar*P)./(mu_bar+m_bar+a_I_bar*P+pf_bar*Z);
        dNdP = -(a_S_bar.*(1-Yvals)+a_I_bar.*Yvals + ...
            -a_I_bar.*(chimu_bar+chim_bar)./(chimu_bar+chim_bar+chixai_bar.*P)*factorZ+...
            chixai_bar.*(mu_bar+m_bar)./(chimu_bar+chim_bar+chixai_bar.*P)*factorZ+...
            (mu_bar+m_bar+a_I_bar*P)./(chimu_bar+chim_bar+chixai_bar.*P).*factorZ.*dxdP*P);
        
        plot(Yvals,dNdP,'-','Color', '#6baed6','LineWidth',1.5)
    end

    % Replicates with predator density = 0.1
    for ii = [2, 3, 22, 47, 48, 49]
        
        P=0.1;
        Freqs = Freq_week2(ii,1:9);
        a_S_bar = sum(Freqs.*a_i);
        a_I_bar = sum(Freqs.*a_i.*omega);
        chimu_bar = sum(Freqs.*chi_i.*mu_i);
        chim_bar = 0; % m_i assumed to be negligible
        chixai_bar = sum(Freqs.*chi_i.*x_i_01pred.*a_i.*omega);
        m_bar = 0; %m_i assumed to be negligible
        mu_bar = sum(Freqs.*mu_i);
        pf_bar = sum(Freqs.*p_i.*fSi); % S and I assumed to have equal uptake
        dxdP = sum(Freqs.*dx_idP_01pred);
        
        Z = (mu_bar+m_bar+a_I_bar*P).*Yvals./pf_bar./(1-Yvals);
        factorZ=(mu_bar-a_S_bar*P+a_I_bar*P)./(mu_bar+m_bar+a_I_bar*P+pf_bar*Z);
        dNdP = -(a_S_bar.*(1-Yvals)+a_I_bar.*Yvals + ...
            -a_I_bar.*(chimu_bar+chim_bar)./(chimu_bar+chim_bar+chixai_bar.*P)*factorZ+...
            chixai_bar.*(mu_bar+m_bar)./(chimu_bar+chim_bar+chixai_bar.*P)*factorZ+...
            (mu_bar+m_bar+a_I_bar*P)./(chimu_bar+chim_bar+chixai_bar.*P).*factorZ.*dxdP*P);
        
        plot(Yvals,dNdP,'-','Color', '#3182bd','LineWidth',1.5)
    end


    % Replicates with predator density = 0.5
    for ii = [19, 20, 21, 27, 28, 29]
        
        P=0.5;
        Freqs = Freq_week2(ii,1:9);
        a_S_bar = sum(Freqs.*a_i);
        a_I_bar = sum(Freqs.*a_i.*omega);
        chimu_bar = sum(Freqs.*chi_i.*mu_i);
        chim_bar = 0; % m_i assumed to be negligible
        chixai_bar = sum(Freqs.*chi_i.*x_i_05pred.*a_i.*omega);
        m_bar = 0; %m_i assumed to be negligible
        mu_bar = sum(Freqs.*mu_i);
        pf_bar = sum(Freqs.*p_i.*fSi); % S and I assumed to have equal uptake
        dxdP = sum(Freqs.*dx_idP_05pred);
        
        Z = (mu_bar+m_bar+a_I_bar*P).*Yvals./pf_bar./(1-Yvals);
        factorZ=(mu_bar-a_S_bar*P+a_I_bar*P)./(mu_bar+m_bar+a_I_bar*P+pf_bar*Z);
        dNdP = -(a_S_bar.*(1-Yvals)+a_I_bar.*Yvals + ...
            -a_I_bar.*(chimu_bar+chim_bar)./(chimu_bar+chim_bar+chixai_bar.*P)*factorZ+...
            chixai_bar.*(mu_bar+m_bar)./(chimu_bar+chim_bar+chixai_bar.*P)*factorZ+...
            (mu_bar+m_bar+a_I_bar*P)./(chimu_bar+chim_bar+chixai_bar.*P).*factorZ.*dxdP*P);
        
        plot(Yvals,dNdP,'-','Color', '#08519c','LineWidth',1.5)
    end
    
    % Replicates with predator density = 1
    for ii =[4, 5, 6, 41, 42, 43]
        
        P=1;
        Freqs = Freq_week2(ii,1:9);
        a_S_bar = sum(Freqs.*a_i);
        a_I_bar = sum(Freqs.*a_i.*omega);
        chimu_bar = sum(Freqs.*chi_i.*mu_i);
        chim_bar = 0; % m_i assumed to be negligible
        chixai_bar = sum(Freqs.*chi_i.*x_i_1pred.*a_i.*omega);
        m_bar = 0; %m_i assumed to be negligible
        mu_bar = sum(Freqs.*mu_i);
        pf_bar = sum(Freqs.*p_i.*fSi); % S and I assumed to have equal uptake
        dxdP = sum(Freqs.*dx_idP_1pred);
        
        Z = (mu_bar+m_bar+a_I_bar*P).*Yvals./pf_bar./(1-Yvals);
        factorZ=(mu_bar-a_S_bar*P+a_I_bar*P)./(mu_bar+m_bar+a_I_bar*P+pf_bar*Z);
        dNdP = -(a_S_bar.*(1-Yvals)+a_I_bar.*Yvals + ...
            -a_I_bar.*(chimu_bar+chim_bar)./(chimu_bar+chim_bar+chixai_bar.*P)*factorZ+...
            chixai_bar.*(mu_bar+m_bar)./(chimu_bar+chim_bar+chixai_bar.*P)*factorZ+...
            (mu_bar+m_bar+a_I_bar*P)./(chimu_bar+chim_bar+chixai_bar.*P).*factorZ.*dxdP*P);
        
        plot(Yvals,dNdP,'-','Color', 'k','LineWidth',1.5)
    end
    
    text(-0.2,0.005,'b','FontSize',16,'FontWeight','Bold')
