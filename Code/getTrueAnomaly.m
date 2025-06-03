function vr = getTrueAnomaly(n, d, l, p, e, a)
    
    M = n * d + l - p;                                                      %Calculate Mean Anomaly
    M = M - floor(M/360) * 360;                                             %Ensure that the Mean Anomaly doesn't exceed 360 degrees
    
    v = M + 180/pi * ((2 *e - e^3 / 4) * sin(M) + ...                       %Calculate True Anomaly
        5 / 4 * e^2 * sin(2*M) + 13 / 12 * e^2 * sin(3*M));
    
    r = a * (1 - e^2) / (1 + e * cos(v));                                   %Calculating the Distance from the first focus (commonly the sun)
    
    if v > 360
        v = v - 360;
    end
    
    vr = [r, v];                                                            %output the distance and the true anomaly
    
end