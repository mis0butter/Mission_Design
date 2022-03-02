clear
clc

JD=2457271.50000;
n=JD-2451545;

% n=5477.5+7442;

L=280.46+0.9856474*n;
Ll=floor(L/360);
L=L-360*Ll;
g=357.528+0.9856003*n;
gg=floor(g/360);
g=g-360*gg;

lambda=L+1.915*sind(g)+0.02*sind(2*g);
B=0;
e=23.439-0.0000004*n;
a=atand(cosd(e)*tand(lambda));

d=asind(sind(e)*sind(lambda));


Thetag=18.697374458+24.06570982*(JD-2451545);
Thetag=Thetag*15;
Tg=floor(Thetag/360);
Thetag=Thetag-360*Tg;

Long=a-Thetag;
Lat=d;

if Long <-180
    Long=Long+180;
end

if Long >180
    Long=Long-180;
end