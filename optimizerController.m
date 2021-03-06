
%---<INNITIAL PARAMETERS>---%
%The innitial parameters represents the innitial rocket design for the
%optimizer to start with, this design should produce a rocket that can fly
%to a non-trivial altitude and is structurally sound 
innitialParams = [0.5, 0.008, 0.5, 0.002, 0.002, 0.5, 0.002, 0.37, 0.4, 0.4, 3, 0.0082, 0.0172, 0.017];

%---<UPPER BOUNDS>---%
%The upper bounds represent the maximum values any parameters are allowed
%to have, many of these considerations are made from a manufacturability
%perspective.
upperBounds = [1, 0.02, 1, 0.002, 0.002, 3, 0.02, 1, 1, 3, 6, 0.02, 0.04, 0.05];

%---<LOWER BOUNDS>---%
%The lower bounds represent the lowest design values the rocket is allowed
%to have. Many of these considerations are made from a manufacturing
%perspective but the most important characteristic is that no values can be
%zero, which may return NaN results from the simulation.
lowerBounds = [0.1, 0.002, 0.1, 0.002, 0.002, 0.1, 0.002, 0.01, 0.02, 0.01, 0.1, 0.002, 0.001, 0.002];

%Innitiate optimization
x = fmincon(@rocketModel, innitialParams, [], [], [], [], lowerBounds, upperBounds);