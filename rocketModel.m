function [score] = rocketModel(tD, tT, t1L, c1T, c2T, t3L, bT, fuelCore, fuelDia, fuelLength, m_ox, nozzleThroat, nozzleExit, chamberT) 
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    dt = 0.4; %[s] integration timestep
    altitude = zeros(1,2000);
    velocity = zeros(1,2000);
    acceleration = zeros(1,2000);
    liftoff = true;

    L = 1; %[m] nosecone length

    %CALCULATE MOTOR THRUST CURVE AND OPTIMAL OXIDIZER AMOUNT ALSO T2 LEN
    thrustC = motorThrust(fuelCore, fuelDia, fuelLength, m_ox, nozzleThroat, nozzleExit)
    burnTime = max(thrustC(:,1));

    oxMass = burnTime*m_ox; %[kg] optimal mass of oxidiser


    parachuteM = 2; %[kg]
    avionicsM = 2; %[kg]
    payloadM = 0; %[kg]
    nozzleM = 0.5; %[kg]

    rhoAl = 2720; %[kg/m^3]
    rhoABS = 1052;%[kg/m^3]

    %OXIDIZER TANK CONSIDERATIONS

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
        + avionicsM + parachuteM + bulkheadMass*2 + nozzleM %calculate mechanical weight

    wetMass = dryMass + oxMass + fuelM %calculate the fully loaded weight of the craft
    
    %BEGIN TIME LOOP%
    i = 2;
    while(velocity(i)>=0 && liftoff==true)

        
        if(dt*i > burnTime)
            motorForce = 0;
            mCur = dryMass;
        else
            motorForce = interp1(thrustC(:,1), thrustC(:,2), dt*i);
            mCur = wetMass - m_ox*dt*i - interp1(thrustC(:,1), thrustC(:,3), dt*i) %estimate the current mass of the rocket where t = dt*i
        end

        %---<UPDATE POSITION>---%
        acceleration(i) = (motorForce - drag(tD, 1.2, velocity(i-1)) - mCur*9.81)/mCur;
        if(acceleration(i)>0)
            liftoff = false;
        end
        velocity(i) = velocity(i-1) + acceleration(i)*dt;
        altitude(i) = altitude(i-1) + velocity(i)*dt;
        acceleration(i);
        i = i + 1;
        if(i>2000)
            break
        end
    end
    %END TIME LOOP%
    
    score = acceleration;

end

