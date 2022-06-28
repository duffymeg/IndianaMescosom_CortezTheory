% Script for computing parameter values

% AUTHOR:   Michael Cortez, Florida State University, mcortez@fsu.edu
% DATE:     June 2021
% PURPOSE:  Compute parameter values 

% INPUTS:   Parameters for different host clones. Order of parameters is:
%           (1) BD05-42     (2) BD08-46	(3) BD19-64	(4) CB24-68	
%           (5) DW22-58(84) (6) DW29-75	(7) IL14-43	(8) ML30-82 (9) ML32-84
%           Fraction of infected individuals from exposure trials
%           Host filtering rates
%           Estimated predator attack rates
%           Estimated host mortality rates
%           Estimated burst sizes for individuals not consumed

% METHOD:   Step 1: Compute per spore probabilities of infection 
%           Step 2: Compute reduction in spore burst size for infected
%           individuals that are consumed by predators, for each level of
%           predator density


% OUTPUTS:  SporesConsumed -- number of spores consumed during exposure exp
%           p_i -- per spore probability of infection
%           x_i_Ypred -- reduction in spore burst size at pred density = Y

clear all

format long

%% Measured Parameters

% Order of clones is 
% Fraction of clones infected in exposure experiments, unitless
    FracInfect = [0.3684 0.45 0 0.375 0.7619 0.6842  0.8333 0.2353 0.6111];
% Filtering rates, L/hr/indiv
    fSi = [1.44 1.74 0.994 1.55 2.46 2.65 1.54 0.58 2.78]*10^(-4);
% Disease-induced mortality rate, 1/hr
    mu = 1/20/24; 
% Predator attack rates on susceptible individuals, 1/hr/indiv
% Increase in attack rate on infected individuals
    a_i = [0.0197 0.0134 0.0116 0.0062 0.0144 0.0082 0.016 0.0157 0.0096];
    omega = 2; 
% Burst sizes
    chi_i= [78231 36990  NaN 70578 89909 146128 108163 41667 92610];
    chi_i(3) = mean(chi_i([1:2,4:end]));
    
%% Step 1: Estimating per spore probabilities
    
    % numbers of spores consumed by a single individual in 0.045 mL of water
    % in 24 hours with initial spore density of 200000 spores/L 
    SporesConsumed = 200000*(1-exp(-fSi*24/0.045))*0.045; 
    p_i = FracInfect./SporesConsumed;
    
    
%% Step 2: Reduction in spore burst size for consumed individuals
    % Let z be the density of spores in a single host
    % Model growth using the logistic equation dz/dt = r(1-z/sigma) with
    % initial condition z(0) = z0
    
    % Based on Ault et al. 2014, assume r = 0.5/day = 0.0208/hr
        r =0.5/24;
    
    % Step A: Compute values for z0 and sigma assuming assuming
    %         z(13 days) = 0.5*sigma; z(20 days) = chi_i 
        z0_i = zeros(1,length(chi_i)); % initial conditions
        sigma_i = zeros(1,length(chi_i)); % maximum burst sizes
        t1 = 13*24; % 13 days * 24 hours
        t2 = 20*24; % 20 days * 24 hours

        for ii=1:length(chi_i)
            z0_i(ii) = chi_i(ii)*(exp(r*t1)+exp(r*t2))/exp(r*t2)/(exp(r*t1)+1);
            sigma_i(ii) = chi_i(ii)*(exp(r*t1)+exp(r*t2))/exp(r*t2);
        end
    
    % Step B: Calculate average time till death after exposure in the
    %         presence of predators
        time0 = 1/mu; % time if predator density is 0/L
        time1 = 1./(mu+omega*a_i*0.1); % time if predator density is 0.1/L
        time2 = 1./(mu+omega*a_i*0.5); % time if predator density is 0.5/L
        time3 = 1./(mu+omega*a_i*1); % time if predator density is 1/L
        
        
    % Step C: Calculate reduction in spore burst size
        % solution to ODE is z(t) = sigma*z0*exp(rt)/(sigma-z0+z0*exp(rt))
        
        % Check for no reduction if predator density is 0/L
        x_i_0pred = sigma_i.*z0_i.*exp(r*time0)./(sigma_i-z0_i+z0_i.*exp(r*time0))...
            ./chi_i;
        % Reduction if predator density is 0.1/L
        x_i_01pred = sigma_i.*z0_i.*exp(r*time1)./(sigma_i-z0_i+z0_i.*exp(r*time1))...
            ./chi_i;
        % Reduction if predator density is 0.5/L
        x_i_05pred = sigma_i.*z0_i.*exp(r*time2)./(sigma_i-z0_i+z0_i.*exp(r*time2))...
            ./chi_i;
        % Reduction if predator density is 1/L
        x_i_1pred = sigma_i.*z0_i.*exp(r*time3)./(sigma_i-z0_i+z0_i.*exp(r*time3))...
            ./chi_i;
        
        

        
        
        
