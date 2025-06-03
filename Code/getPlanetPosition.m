function position = getPlanetPosition(r, v)

    x = r * cos(v * (pi/180));
    y = r * sin(v * (pi/180));
    
    position = [x, y];  
end