function day = getDay(earthDay, earthYear, marsYear, earthLandingDay, marsLandingDay)

marsYear1 = marsYear;                                                       %Mars's initial year 
earthYear1 = earthYear;                                                     %Earth's initial year

earthYear = earthYear + floor(earthDay / 365);                              %get true year

if floor(rem(earthYear, 4)) ~= 0                                            % Checking if the year is a leap year
    earthLaunchDay = earthDay - (earthYear - earthYear1)*365;                    %launchDay is not allowed to be bigger than the days in the year
else
    earthLaunchDay = earthDay - (earthYear - earthYear1)*365 + ...          %launchDay on a leap year
                                                 (earthYear-earthYear1)/4;
end

marsYear = marsYear + floor(earthDay / 687);                                % get true mars year
marsDay = (earthDay - 5);                                                   % - 5 because masr is 5 days behind on its orbit


if (marsDay - (marsYear - marsYear1)*687) < 0                                                        % ensuring that mars day is a whole number
    marsDay = marsDay + 687;                                                % adding days from mars's previous year
    marsYear = marsYear - 1;                                                % ensuring that mars year is correct
end

if marsYear ~= 0
    marsLaunchDay = marsDay - (marsYear - marsYear1)*687;                       % get mars launch day
else
    marsLaunchDay = marsDay;
end

if ~exist('earthLandingDay', 'var')
    earthLandingDay = 0;
end

if ~exist('marsLandingDay', 'var')
    marsLandingDay = 0;
end

moonRev = floor(earthDay / 27);
moonDay = earthLaunchDay - 27 * moonRev;

day = [earthDay, earthLaunchDay, earthLandingDay ;marsDay, marsLaunchDay,...
    marsLandingDay; moonDay, 0, 0; earthYear, marsYear, 0];

end
