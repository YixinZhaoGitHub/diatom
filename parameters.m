function p=parameters(n,NC)
%

p.vDp = 15; % advection velocity of non-dia phyto detritus m/d
p.vDd = 20;
p.diffml=100;%mixed layer diff     diffusivity m2/day 
p.diffdl=5;%deep layer diff
p.mld=50;% for constant simulation
p.w=1;% themocline width
% water column
p.Zm=300; % depth
p.dz = 2; % grid cell size
p.z = (0+p.dz/2):p.dz:(p.Zm-p.dz/2); %grid vector
p.xgrid=length(p.z); % grid cell number
%
p.Iin = 200; %incident light intensity   W/m2
p.k = 0.07; %light attenuation m2 (mmol N)-1   0.05
p.Kbg = 0.05;  % background turbidity m-1   0.0375
p.pmax=0.7; % maximal specific production rate   day-1
p.l = 0.03; % intrinsic loss rate d-1
p.Al= 0.12; % affinity of light [day-1]/[W/m2]  beckmann 2007  0.12
p.An= 1.67; % affinity of nutrient   d-1/[mmol N m-3]  beckmann 2007  1.67
p.tau= 0.1; % remineralization coefficient
p.gamma= 1.5; %loss by zooplankton  m3 (mmol N)-1 d-1
p.NB=400;

%----------------------------nutrient conditions---------------
  switch NC
      case 1 % Si:N 1
  p.NB=p.NB; 
  p.SB=p.NB*1; 
      case 2 % Si:N 0.7
  p.NB=p.NB; 
  p.SB=p.NB*0.7;
      case 3 % Si:N 0.4
  p.NB=p.NB;  
  p.SB=p.NB*0.4; 
      case 4 % Si:N 0.1
  p.NB=p.NB;  
  p.SB=p.NB*0.1;% 
  end

%diatom parameter
% n is group number ----- from 'defense specialist' to 'competitor specialist'
p.pmaxd=linspace(0.9,1.3,n);
p.gammad=linspace(0.5,1.2,n); % predation loss coefficient  gamma*Biomass^2
A=getAffinity;
p.Ald=A.light;%[0.1000 0.1500 0.2500 0.5000 1.5000]
p.And=A.nitrogen;% m3 mmol N-1 d-1
p.Asd=A.silicon;% m3 mmol Si-1 d-1
p.vd=linspace(1,0,n);%
p.rhoSN=1.05;%  molar ratio Si:N  mmol Si/mmol N

if n==0 %defense specialist
  p.pmaxd=1.1;
  p.gammad=0.9;  %predation loss coefficient  gamma*Biomass^2
  p.Ald=p.Al*1.3;
  p.And=p.An*1.2;
  p.vd=2;
  p.pmaxd(5)=0;
  p.gammad(5)=0; % predation loss coefficient  gamma*Biomass^2
  p.Ald(5)=0;
  p.And(5)=0;
  p.vd(5)=0;
  p.rhoNS(1:5)=1;%

elseif n>0 && n<5
  p.pmaxd(5)=0;
  p.gammad(5)=0; % predation loss coefficient  gamma*Biomass^2
  p.Ald(5)=0;
  p.And(5)=0;
  p.vd(5)=0;
  p.rhoNS(n+1:5)=1;%

else
end

end