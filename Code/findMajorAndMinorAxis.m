function answer = findMajorAndMinorAxis(aphilion, perhilion, e)
    a = (aphilion + perhilion)/2;
    b = (sqrt(a^2 * (1 - e^2)));
    answer = [a b];
end