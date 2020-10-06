%% ENPH 455 Hybrid Rocket Motor Design Project (Iterative Code)
%Written by: Cole Mero 
%October 29, 2019
%
%function [regRate, m_fuel, G_tot] = RegRate(m_ox, m_fuel, currentD, rho, Lp, dt) 
%
%Takes in the oxidizer mass flow rate, fuel mass flow rate, the current
%diameter, the density of the solid fuel, the length of the combustion
%chamber, and the differential time step. Returns the instantaneous
%regression rate, fuel mass flow rate, and propellant mass flux. 

%Note that references to SPAD are the textbook Space Propuslsion Analysis
%and Design by Humble. 

%% Function to Calculate Instantaneous Regression Rate + Fuel Mass Flow Rate

function [regRate, m_fuel, G_tot] = RegRate(m_ox, currentD, rho, Lp, AmplificationFactor)

% a = 2.066*10^(-5);   %[kg^(-n)m^(1+2n-m)s^(n-1)] (reference value found on pg 385 of SPAD) 
% n = 0.75;              %(reference value found on pg 385 of SPAD) 
% m = -0.15;             %(reference value found on pg 385 of SPAD) 
% 
% G_tot = (m_ox + m_fuel)/currentA; %[kg/m^2*s] Calculate the instantaneous total mass flux 
% 
% 
% regRate = a*(G_tot^n)*Lp^m; %[m/s] Calcualte the instantaneous regression rate
% m_fuel = regRate*rho*pi*currentD*Lp*(dt*regRate); %[kg/s] Calculate the instantaneous fuel mass flow rate

currentA = pi*(currentD/2)^2;  %[m^2]      Calculate the current combustion port cross sectional area

%The below constants are derived for rather specific and stupid units for
%the terms in the regression rate law, example values for these constants
%can be found in SPAD Table 7.2 Pg 389. Typically these values are found
%from experimental data, thus I used a research paper and played with the
%values to emulate their results. Research paper:
%https://pdfs.semanticscholar.org/b004/f9574f56c810273f5bbf7e494ad3598abdae.pdf

Ao = 0.000007; 
No = 0.8; 
Mo = -0.20;

Gox = (m_ox/currentA); %[kg/m^2*s] Oxidizer flux in units for the regression rate calc;
regRate = Ao*(Gox^No)*(Lp^Mo) * AmplificationFactor; %[m/s] Regression rate calculation (eqn 7.36 from SPAD pg 385)

m_fuel = regRate*rho*(pi*currentD*Lp); %[kg/s] Instantaneous fuel mass flow rate
Gfuel = m_fuel/currentA; %[kg/m^2*s] Instantaneous fuel mass flux 

G_tot = Gox + Gfuel; %Calculates the instantaneous total mass flux 

end 