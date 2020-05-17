function [score] = rocketModel(input) 
%ROCKET MODEL -> Hybrid rocket simulation software
%   This function takes various parameters and models the flight of a
%   hybrid rocket from them. The input parameters are as follows.
    input
    tD = input(1); %Outer diameter of the rocket 
    tT = input(2); %Tube thickness of the rocker
    t1L = input(3); %Length of tube 1, (Avionics and recovery)
    c1T = input(4); %Coupler 1 thickness 
    c2T = input(5); %Coupler 2 thickness
    t3L = input(6); %Tube 3 length
    bT = input(7); %Oxidizer bulkhead thickness
    fuelCore = input(8); %Fuel core diameter 
    fuelDia = input(9); %Fuel core max diameter
    fuelLength = input(10); %Length of ABS fuel core
    m_ox = input(11); %Oxidizer mass flow rate
    nozzleThroat = input(12); %Nozzle throat diameter
    nozzleExit = input(13); %Nozzle exit diameter
    chamberT = input(14); %Fuel chamber thickness 
    

    dt = 0.5; %[s] integration timestep
    altitude = zeros(1,2000);
    velocity = zeros(1,2000);
    acceleration = zeros(1,2000);
    liftoff = true;
    
    %INNITIALIZE FLIGHT STRESS MONITORS
    stressXt1 = zeros(1,2000); %Stress in tube 1 pointing radially
    stressYt1 = zeros(1,2000); %Stress in tube 1 pointing axially
    
    stressXt2 = zeros(1,2000); %Stress in tube 2 pointing radially
    stressYt2 = zeros(1,2000); %Stress in tube 2 pointing axaially
    
    stressXt3 = zeros(1,2000); %Stress in tube 3 pointing radially
    stressYt3 = zeros(1,2000); %Stress in tube 3 pointing axially
    %END OF FLIGHT STRESS MONITOR INNITIALIZATION

    L = 1; %[m] nosecone length

    %CALCULATE MOTOR THRUST CURVE AND OPTIMAL OXIDIZER AMOUNT ALSO T2 LEN
    thrustC = motorThrust(fuelCore, fuelDia, fuelLength, m_ox, nozzleThroat, nozzleExit, dt);
    burnTime = max(thrustC(:,1))

    oxMass = burnTime*m_ox; %[kg] optimal mass of oxidiser


    parachuteM = 2; %[kg]
    avionicsM = 2; %[kg]
    payloadM = 0; %[kg]
    nozzleM = 0.5; %[kg]

    rhoAl =1300; %[kg/m^3]
    yeildAl = 1300*10^6; %[N/m^2] Carbon yeild stress UPDATE, SWAPPED TO CARBON
    sf = 0.99; %safety factor
    yeildAl = yeildAl*sf;
    rhoABS = 1052;%[kg/m^3]

    %OXIDIZER TANK CONSIDERATIONS

    pNox_20c = 6315 * 10^3; %[N/m^2]

    rhoNoxL_20c = 743.9; %[kg/m^3] density before tank drain ~20c liquid
    rhoNoxG_20c = 190.0; %[kg/m^3] density before tank drain ~20c vapour

    %Assume u = h since tank is a constant volume system
    uNoxL_20c = -241; %[kj/kg] specific internal energy. 
    uNoxG_20c = -94.4; %[kj/kg] specific internal energy. 

    x = @(Utot, mTot, uVap, uLiq) ((Utot/mTot)-uLiq) / (uVap - uLiq);
    tankVol = @(mass, xCur, rhoL, rhoG) mass * ((1-xCur)/rhoL + xCur/rhoL);

    % Assuming ~100% tank fill with liquid at startup Utot = uLiq*mTot
    xT0 = x(oxMass*uNoxL_20c, oxMass, uNoxG_20c, uNoxL_20c);
    vol = tankVol(oxMass, xT0, rhoNoxL_20c, rhoNoxG_20c);

    volBh = 2*(4/24)*pi*(tD-2*bT)^3; %Volume added by the two bulkheads
    lOx = ((vol-volBh) / ((tD-2*tT)^2 * (pi/4))); %[m] length of the ox tank in meters
    t2L = lOx + 2*tD; %[m] length of entire tube 2 (accounting for coupling)

    %END OF OXIDIZER TANK CONSIDERATIONS


    %AERODYNAMICS CONSIDERATIONS
    cD = 0.5; %Drag coeffecient, to be replaced later with better estimates
    drag = @(D, rho, v) rho*v^2*(1/2)*(pi/4)*D^2;
    airDensity = 1.2; %[kg/m^3] Air density (starting at alt = 0)
    externalTemp = 30; %[C] Air temp (holding variable)
    groundTemp = 30; %[C] Ground temp ~30 Spaceport
    %END OF AERODYNAMICS CONSIDERATIONS

    %MASS CONSIDERATIONS
    tubeMass = @(D, T, L) (1/4)*pi*(D^2 - (D-2*T)^2)*L *rhoAl;
    bulkheadMass = (4/24)*pi*(tD^3 - (tD-2*bT)^3)*rhoAl;
    noseconeMass = ((pi*L)/12) *(tD^2 - (tD-2*tT)^2)*rhoAl;

    t1M = tubeMass(tD, tT, t1L);
    t2M = tubeMass(tD, tT, t2L);
    t3M = tubeMass(tD, tT, t3L);

    c1M = tubeMass(tD-2*tT, c1T, 2*tD);
    c2M = tubeMass(tD-2*tT, c2T, 2*tD);

    chamberM = tubeMass(fuelDia, chamberT, fuelLength);
    fuelM = (1/4)*pi*(fuelDia^2 - fuelCore^2)*fuelLength*rhoABS;

    dryMass = t1M + t2M + t3M + c1M + c2M + chamberM + noseconeMass + payloadM...
        + avionicsM + parachuteM + bulkheadMass*2 + nozzleM; %calculate mechanical weight

    wetMass = dryMass + oxMass + fuelM; %calculate the fully loaded weight of the craft
    
    %BEGIN TIME LOOP%
    i = 2;
    while(i<500)

        
        if(dt*i > burnTime)
            motorForce = 0;
            mCur = dryMass;
        else
            motorForce = thrustC(i,2);
            mCur = wetMass - m_ox*dt*i - thrustC(i,3); %estimate the current mass of the rocket where t = dt*i
        end
        
        %---<UPDATE POSITION>---%
        acceleration(i) = (motorForce - drag(tD, airDensity, velocity(i-1)) - mCur*9.81)/mCur;

        
        if(acceleration(i)>0)
            liftoff = false;
        end

        velocity(i) = velocity(i-1) + acceleration(i)*dt;
        altitude(i) = altitude(i-1) + velocity(i)*dt;
        
        if(altitude(i) < 11000) %TROPOSPHERE
            externalTemp = 15.04 - 0.00649*altitude(i);
            externalPres = 101.29 * ((externalTemp + 273.1)/288.08)^5.256;
            airDensity = externalPres/(0.28705*(externalTemp + 273));
            
        elseif(altitude(i) >= 11000 && altitude(i) < 25000) %STRATOSPHERE (LOWER)
            externalTemp = -56.46;
            externalPres = exp(1.73 - 0.000157*altitude(i));
            airDensity = externalPres/(0.28705*(externalTemp + 273));
            
        else
            externalTemp = -131.21 + 0.00299*altitude(i); %STRATOSPHERE (UPPER)
            externalPres = 2.488 * ( (externalTemp+273.1)/216.6)^(-11.388);
            airDensity = externalPres/(0.28705*(externalTemp + 273));
        end
        
        %---<CALCULATE STRESSES>---%
        stressXt1(i) = (((noseconeMass + parachuteM + t1M + payloadM)*(acceleration(i))) + drag(tD, 1.2, velocity(i))) / ((pi/4)*(tD^2 - (tD-tT)^2)); %calculate stress on tube due to mass acceleration and nosecone drag
        stressYt1(i) = 0;

        stressXt2(i) = ((m_ox*dt*i + 2*bulkheadMass + t2M + avionicsM)*acceleration(i) + stressXt1(i)) / ((pi/4)*(tD^2 - (tD-tT)^2)); %Calculate stress X stress on oxidizer tube 
        stressYt2(i) = (pNox_20c * (tD-2*tT))/(2*tT); %hoop stress on tank walls from pressurized oxidizer

        stressXt3(i) = ((motorForce)*acceleration(i)) / ((pi/4)*(tD^2 - (tD-tT)^2)); % axial stress on 
        stressYt3(i) = (thrustC(i,4)* (fuelDia-2*chamberT))/(2*chamberT);
        
         i = i + 1;
        if(i>2000)
            break
        end


        
    end
    %END TIME LOOP%
    
    %---<COMPUTE SCORE>---%
    verbose = [acceleration', velocity', altitude', stressXt1', stressYt1', stressXt2', stressYt2', stressXt3', stressYt3'];
    safetyFactors = [max(stressXt1)/yeildAl, max(stressYt1)/yeildAl, max(stressXt2)/yeildAl, max(stressYt2)/yeildAl, max(stressXt3)/yeildAl, max(stressYt3)/yeildAl];
    
    
    if(max(safetyFactors > 1))
        score = abs((100000 - max(altitude)) + mean(safetyFactors) * 1000);
        "failure expected"
    else
        score = abs(100000 - max(altitude));
        "failure not expected"
    end
    max(altitude)
    plot(acceleration(1:500))
    %plot(thrustC(1:500,2))
    %printf("MAX ALTITUDE - > %d m", max(altitude));
    %---<SCORE COMPUTED>---%

end

