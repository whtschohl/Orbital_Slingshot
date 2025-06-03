clc, clear, close all;


%% initialising values
%==========================================================================

%Years
%--------------------------------------------------------------------------
earthYear = 2020;
marsYear = 1;                                                               %Count of mars years starting at earth year 2020 

%the Earth
%--------------------------------------------------------------------------
earthAph = kmToAU(152097701);                                               %Earth's Aphelion
earthPer = kmToAU(147098074);                                               %Earth's Perihelion
earthE = 0.0166967;                                                         %Earth's Eccentricity
earthN = 0.9855796;                                                         %Earth's Daily Motion
earthL = 328.40353;                                                         %Earth's Mean Longitude 
earthP = 102.8514;                                                          %Earth's Longitude at Perihelion

%Mars
%--------------------------------------------------------------------------
marsAph = kmToAU(249*10^6);                                                 %Mars' Aphelion
marsPer = kmToAU(206*10^6);                                                 %Mars' Perihelion
marsE = 0.0934231;                                                          %Mars' Eccentricity
marsM = 6.39*10^23;                                                         %Mars' Mass
marsRad = kmToAU(3389.5);                                                   %Mars' Radius
marsN = 0.5240613;                                                          %Mars' Daily Motion
marsL = 262.42784;                                                          %Mars' Mean Longitude
marsP = 3360.882;                                                           %Mars' Longitude at Perihelion

%the Moon
%--------------------------------------------------------------------------
moonAph = kmToAU(405696);                                                   %Moon's Aphelion
moonPer = kmToAU(363104);                                                   %Moon's Perihelion
moonE = 0.0549;                                                             %Moon's Eccentricity
moonM = 0.07346*10^24;                                                      %Moon's Mass
moonRad = kmToAU(1737.4);                                                   %Moon's Volumetric Radius
moonN = 13.17619501;                                                        %Moon's Daily Motion
moonL = (moonPer + moonAph) / 2;                                            %Moon's Mean Longitude

%% finding Axes
%==========================================================================
earthAxes = findMajorAndMinorAxis(earthAph, earthPer, earthE);              %Getting Earth's Axes
marsAxes = findMajorAndMinorAxis(marsAph, marsPer, marsE);                  %Getting Mars' Axes
moonAxes = findMajorAndMinorAxis(moonAph, moonPer, moonE);                  %Getting Moon's Axes

%% finding orbital paths of the Earth and Mars
%==========================================================================
earthOrbit = orbitalPathFinder(earthAxes(1), earthAxes(2), (earthAxes(1)... %getting earths orbital path
    - earthPer), 0);
marsOrbit = orbitalPathFinder(marsAxes(1), marsAxes(2), (marsAxes(1) -...   %getting mars's orbital path
    marsPer), 0);

%% Position of Celestial bodies
%==========================================================================

%Earth Position
%--------------------------------------------------------------------------
zerosVector = zeros([1, 2]);                                                %initialising Vector to decrease RAM usage
if floor(rem(earthYear, 4)) ~= 0                                            %checking if year is a leap year
    earthTrueAnomaly = zeros([365, 2]);                                     %initialising matrix to decrease RAM usage
    earthPos = zeros([365, 2]);                                             %initialising Vector to decrease RAM usage

    for i = 1:365                                                           %for every day in a non-leap year
        earthTrueAnomaly(i, :) = getTrueAnomaly(earthN, i, earthL, ...      %getting distance from the san and the true anomaly (angle between earth(x,y) and sun(0, 0)
            earthP, earthE, earthAxes(1));
        earthPos(i, :) = getPlanetPosition(earthTrueAnomaly(i, 1),...       %assinging x y co-ordinates for all earth positions
            earthTrueAnomaly(i, 2));
    end
                                                                            %changing indexes of positions so that day 1 = index
    for i = 1:365                                                           %
        for j = (i+1):365                                                   %
            if earthTrueAnomaly(i, 2) > earthTrueAnomaly(j, 2)              %finding smaller true anomaly (smallest true annomaly is the first day of the year)
                zerosVector(:) = earthTrueAnomaly(i, :) ;                   %assigning largest value to a temperary vector
                earthTrueAnomaly(i, :) = earthTrueAnomaly(j, :);            %assigning smallest value to the smallest index 
                earthTrueAnomaly(j, :) = zerosVector(:);                    %assigning largest value to the largest index
                                                                            %same process for the position
                zerosVector(:) = earthPos(i, :) ;
                earthPos(i, :) = earthPos(j, :);
                earthPos(j, :) = zerosVector(:);
            end
        end
    end
else                                                                        %same process as above for if year is a leap year
    earthTrueAnomaly = zeros([366, 2]);                                     
    earthPos = zeros([366, 2]);

    for i = 1:366
        earthTrueAnomaly(i, :) = getTrueAnomaly(earthN, i, earthL,...
            earthP, earthE, earthAxes(1));
        earthPos(i, :) = getPlanetPosition(earthTrueAnomaly(i, 1),...
            earthTrueAnomaly(i, 2));
    end
    for i = 1:366
        for j = (i+1):366
            if earthTrueAnomaly(i, 2) > earthTrueAnomaly(j, 2)
                zerosVector(:) = earthTrueAnomaly(i, :) ;
                earthTrueAnomaly(i, :) = earthTrueAnomaly(j, :);
                earthTrueAnomaly(j, :) = zerosVector(:);
                
                zerosVector(:) = earthPos(i, :) ;
                earthPos(i, :) = earthPos(j, :);
                earthPos(j, :) = zerosVector(:);
            end
        end
    end
end

%Mars Position 
%--------------------------------------------------------------------------
marsTrueAnomaly = zeros([687, 2]);                                          %initialising matrix to decrease RAM usage
marsPos = zeros([687, 2]);                                                  %initialising matrix to decrease RAM usage

for i = 1:687                                                               %for i goes through each one of mars's days
    marsTrueAnomaly(i, :) = getTrueAnomaly(marsN, i, marsL, marsP, ...      %same process as above for assigning values to the matrixes
                                                    marsE, marsAxes(1));
    marsPos(i, :) = getPlanetPosition(marsTrueAnomaly(i, 1),...
                                    marsTrueAnomaly(i, 2));
end
for i = 1:687                                                               %same process as above for switching the order of values 
        for j = (i+1):687
            if marsTrueAnomaly(i, 2) > marsTrueAnomaly(j, 2)
                zerosVector(:) = marsTrueAnomaly(i, :) ;
                marsTrueAnomaly(i, :) = marsTrueAnomaly(j, :);
                marsTrueAnomaly(j, :) = zerosVector(:);
                
                zerosVector(:) = marsPos(i, :) ;
                marsPos(i, :) = marsPos(j, :);
                marsPos(j, :) = zerosVector(:);
            end
        end
end

%% Days
%==========================================================================
days = getDay(130, earthYear, marsYear);                                      %returns a 4x3 martix with values of the days, launch days, landing days, and years for earth, mars, and the moon

%% Moon Orbit
%==========================================================================
moonOrbit = orbitalPathFinder(moonAxes(1), moonAxes(2), ...                 %getting the moons orbit
    earthPos(days(1, 2), 1), - earthPos(days(1, 2), 2));

%% trajectory of the Rocket to Mars
%==========================================================================

%avg speed of the rocket in km/s
%--------------------------------------------------------------------------
rocketSpeed = 100/36;                                                      

%Find mars day where mars true anomaly v is = earths true anomaly + 180 
%--------------------------------------------------------------------------
[minValue, days(2, 3)] = min(abs(marsTrueAnomaly(:, 2) - ...
    (earthTrueAnomaly(days(1, 2), 2) + 180).'));

%Get Crusial Values to fin the orbital path of the trajectory
%--------------------------------------------------------------------------
t = earthTrueAnomaly(days(1, 2), 2) * pi / 180;                             %get true anomaly in radians on launch day
trajeR = earthTrueAnomaly(days(1, 2), 1);                                   %get total distance between earth and sun on launch day
trajmR = marsTrueAnomaly(days(2, 3), 1);                                    %get total distance between mars and teh sun on landing day
trajA = (trajeR + trajmR)/2;                                                %get the semi-major axis of the trajectory
trajC = trajA - trajeR;                                                     %get the focaldistance of the trajectory
trajB = sqrt(trajA^2 - trajC^2);                                            %get the semi minor axis of the trajectory
trajCX = trajC * cos(t);                                                    %get the x co-ordinate of the centre of the trajectory
trajCY = trajC * sin(t);                                                    %get the y co-ordinate of the centre of the trajectory
v = earthTrueAnomaly(days(1, 2), 2) * pi/180;

%Use Keplars third law to find landing days and travel time
%--------------------------------------------------------------------------
periodOfJourneyToMars = ceil(ceil(sqrt(trajA^3) * 365) / 2)                 %get total travel time of the rocket using keplars thrid law
idealMarsDay = days(2, 3) - periodOfJourneyToMars;                          %find the ideal position of mars for a perfect hermann transfer

%trajectory Orbit and its rotation
%--------------------------------------------------------------------------
projOrbit = orbitalPathFinder(trajA, trajB, trajCX, trajCY, v);             %get the orbit of the trajectory



%% Slingshot
%==========================================================================
%Get Position of the moon
moonTrueAnomaly = zeros([27, 2]);                                           %initialising matrix to decrease RAM usage
moonPos = zeros([27, 2]);                                                   %initialising matrix to decrease RAM usage

for i = 1:27
    moonTrueAnomaly(i, :) = getTrueAnomaly(moonN, i, moonL, moonPer,...     %same process of getting values as for mars and earth
        moonE, moonAxes(1));
    moonPos(i, :) = getPlanetPosition(moonTrueAnomaly(i, 1),...
                                    moonTrueAnomaly(i, 2));
end
for i = 1:27                                                                %same process of assigning each value to the correct day as above
        for j = (i+1):27
            if moonTrueAnomaly(i, 2) > moonTrueAnomaly(j, 2)
                zerosVector(:) = moonTrueAnomaly(i, :) ;
                moonTrueAnomaly(i, :) = moonTrueAnomaly(j, :);
                moonTrueAnomaly(j, :) = zerosVector(:);
                
                zerosVector(:) = moonPos(i, :) ;
                moonPos(i, :) = moonPos(j, :);
                moonPos(j, :) = zerosVector(:);
            end
        end
end

% intersection of moon orbit and projectory
%--------------------------------------------------------------------------
f = 10.^-2;                                                                 %the decimal that we need rounded to
MO = round(f*moonOrbit(:))/f;                                               %rounding the moons orbit values to the nearest f
PO = round(f*projOrbit(:))/f;                                               %rounding the trajectory orbit values to the nearest f

[k, m] = find(MO == PO);                                                    %finding the position where the moon is closest to the trajectory orbit

%initialising values
%--------------------------------------------------------------------------
t = moonTrueAnomaly(days(3, 1), 2) * pi / 180;                              %same method for initialising the values for the projectory 
a = (moonTrueAnomaly(days(3, 1), 1) + k(2)) / 2; 
c = a - moonPer;
b = sqrt(a^2 - c^2);
slingCX =  c * cos(t);
slingCY =  c * sin(t); 

%Find sling shot orbital path
%--------------------------------------------------------------------------
slingOrbit = orbitalPathFinder(a, b, slingCX, slingCY, t);                  % finding slingshot orbit

% Get time to get around the earth
%--------------------------------------------------------------------------
periodOfSlingShot = ceil(ceil(sqrt(a^3) * 27) / 2)                          %time to travel around the earth
totalPeriodOfJourney = periodOfJourneyToMars + periodOfSlingShot

%% Plotting the Orbit and Planets
%==========================================================================
marsSecFocus = marsE*marsAxes(1);                                           %getting second focus of the mars orbit   

xAxes = axes('Xlim', [-2 2], 'XTick', -2:0.2:2, 'NextPlot', 'add');         %setting limits for the x and y axes in AU
yAxes = axes('Ylim', [-2 2], 'YTick', -2:0.2:2, 'NextPlot', 'add');

figure(1)
hold on                                                                     %to plot all orbits and plannets on one graph

%Plot the celetial bodies
%--------------------------------------------------------------------------
plot(0, 0, 'y*');                                                           %Plotting the sun
plot(earthPos(days(1, 2), 1), earthPos(days(1, 2), 2), 'ko')             	%Plotting the earth
plot(marsPos(days(2, 2), 1), marsPos(days(2, 2), 2), 'mo')                  %Plotting mars's ideal position
plot(marsPos(idealMarsDay, 1), marsPos(idealMarsDay, 2), 'ro')            	%Plotting Mars
plot(moonPer + earthPos(days(1, 2), 1), earthPos(days(1, 2), 2), 'c*')      %Plotting the Moon
plot(marsPos(days(2, 3), 1), marsPos(days(2, 3), 2), 'go')                  %Plotting Mars on landing day

%Plot Orbits
%--------------------------------------------------------------------------
plot(-earthOrbit(1, :),earthOrbit(2, :), 'b')                               %Plotting the earth's orbit
plot(-marsOrbit(1, :)-marsSecFocus, marsOrbit(2,:), 'r')                    %Plotting mars's orbit
plot(-moonOrbit(1, :), moonOrbit(2, :), 'm')                                %Plotting the moon's orbit
plot(projOrbit(1, :),projOrbit(2,:), 'g')                                   %Plotting the trajectories orbit

%Plots for testing all posible positions of earth and mars
%--------------------------------------------------------------------------
% for i = 1:365
%     plot(earthPos(i, 1), earthPos(i, 2), 'm*');
% end
% 
% for i = 1:687
%     plot(marsPos(i, 1), marsPos(i, 2), 'g*');
% end

hold off

%% Calculating Required Fuel
%==========================================================================
%Initialising Rocket Values
%--------------------------------------------------------------------------
m0 = 285228;                                                                %in kilograms
m = 1070;                                                                   %in kilograms

earthVe = 11.2;                                                             %velocity to get to earths orbit in km/s
moonVe = 2.38;                                                              %velocity to get to the moons orbit in km/s

%Using the rocket equation
%--------------------------------------------------------------------------
earthDeltaV = earthVe * log(m0 / m)                                         %calculating the deltaV (fuel needed) to escape earths atmosphere
moonDeltaV = moonVe * log(m0 / m)                                           %calculating the deltaV (fuel needed) to escape moons atmosphere

%% Finishing Touches
%==========================================================================
%add legends
%--------------------------------------------------------------------------
legend('sun', 'earth', 'mars', 'ideal mars position', 'moon', ...
        'mars on landing day', ...
        'location', 'northwest', ...
        'Box', 'off')
    
title('slingshot of a rocket around the earth towards mars')