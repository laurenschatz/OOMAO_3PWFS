function radian = arcsec(val)
if ischar(val)
    val = str2double(val);
end
r=1;
radian = val*constants.arcsec2radian;
end
