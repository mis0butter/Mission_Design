clear
clc

%% Constants
mu=3.986004415e14;
ae=6378136.3;
we=7.292115e-5;
g=9.81;
J2=1.082e-3;
ws=2*pi/365.2422/24/60/60;

%% orbit
alt=350000;
a=alt+ae;
e=0;
I=35;
Long=+5.157;
n=sqrt(mu/a^3);

%% O Precession Calcs
Odot=-(3/2)*n*(ae/a)^2*J2*(1/(1-e^2)^(1/2))*cosd(I); % O precession
Odotday=Odot*3600*24; % for day
wsd=ws*3600*24; % Sun Precession for day
Cs=2*pi/(Odotday-wsd); % Sun Cycle period days
tc=20+43/60+47/3600; % clock time to decimal
LMTT=tc+(Long/15); % Local time at crossing TRMM
LMTR=22+20/60; % Local time Resurs
time_diff=(LMTT-LMTR); % Difference time
O_diff_change_day=24/Cs; %initial difference in time
Offset=time_diff/(O_diff_change_day); % Days since last cross

k=0:round(abs(365.2425/Cs))-1; % amount of crosses in one year
Jk=round(21-Offset)+fix(k*abs(Cs)); % offset by initial cross beforehand and set to cycle every Solar cycle
