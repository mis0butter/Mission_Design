%% Problem 2 

clear 
clc 

JD=2449141.61285;
D=JD-2451545;

% n=5477.5+7442;

% mean longitude of sun 
L=280.46+0.9856474*D;
Ll=floor(L/360);
L=L-360*Ll;

% Mean anomaly 
g=357.528+0.9856003*D;
gg=floor(g/360);
g=g-360*gg;

% ecliptic longitude 
lambda=L+1.915*sind(g)+0.02*sind(2*g);

% ecliptic latitude 
B=0; 

% obliquity of ecliptic 
e=23.439-0.0000004*D; 

% right ascension 
a=atand(cosd(e)*tand(lambda)); 

% declination 
d=asind(sind(e)*sind(lambda)); 

% Earth orientation 
Thetag=18.697374458+24.06570982*(JD-2451545); 

% convert into degrees 
Thetag=Thetag*15; 
Tg=floor(Thetag/360); 
Thetag=Thetag-360*Tg; 

% convert into longitude and latitude 
Long=a-Thetag;
Lat=d;

if Long <-180
    Long=Long+180;
end

if Long >180
    Long=Long-180;
end