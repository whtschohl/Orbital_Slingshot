%% Function to calculate orital path of a planet (elipsis)
function out = orbitalPathFinder(a, b, x0, y0, theata)
    
    % a is the semi-major axis of the orbital path
    % b is the semi-minor axis of the orbital path
    % e is the eccentricity of the orbital path    
    % x0 is the distance between the centre of the orbital path and the 
    %       earth on the x axis only. That means that it is equal to the 
    %       semi-major axis - the perhelion on 01.01    
    % y0 is the distance between the centre of the orbital path and the
    %       earth on the the y axis only. That means that it is equal to 
    %       0 on 01.01 
    
    t=-pi:0.01:pi;
    
    if ~exist('theata', 'var')
        x = a*cos(t) - x0;
        y = b*sin(t) - y0;
        out = [x; y];
    else
        x = a*cos(t + theata) - x0;
        y = b*sin(t + theata) - y0;
        out = [x; y];
    end
end