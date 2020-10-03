%% NOT A PIPE BOMB rocket design 
%Written by: Thomas Kulin 
%October 2, 2020
%
%function [regRate, m_fuel, G_tot] = HelixRegRate(m_ox, m_fuel, currentD, rho, Lp, dt) 
%
%Takes in the oxidizer mass flow rate, fuel mass flow rate, the current
%diameter, the density of the solid fuel, the length of the combustion
%chamber, and the differential time step. Returns the instantaneous
%regression rate, fuel mass flow rate, and propellant mass flux. 

%Note that references are to Engineering Model for Hybrid Fuel Regression
%Rate Amplification Using Helical Ports by Stephen A whitmore


%% Function to Calculate Instantaneous Regression Rate + Fuel Mass Flow Rate

function [regRate, m_fuel, m_total] = HelixRegRate(time, m_ox, r_L)


% A_throat : choked nozzle throat area [cm^2]
% P_0 : combustion chamber pressure [kPa]
% gamma : ratio of specific heats (from RPA)
% R_g : gas constant for combustion products (from RPA) [J/(kg*K)]
% T_0 : combustion flame temperature (from RPA) [K]
% rho_fuel : solid fuel grain material density [kg/(m^3)]

m_total = A_throat * P_0 * sqrt(gamma / (R_g * T_0) * (2/(gamma+1))^((gamma+1)/(gamma-1)));  %nozzle exit mass flow rate [kg/s]
m_fuel = m_total - m_ox;  %fuel mass flow rate [kg/s]

regRate = m_fuel/(2*pi*rho_fuel*r_L*L);  %mean longitudinal regression rate [cm/s]




end 