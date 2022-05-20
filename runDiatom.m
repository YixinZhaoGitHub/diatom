%
% Use RunandPlot.m
%

%
function  sim=runDiatom(n,NC,agg)

arguments
n {mustBeInteger} = 5; %0=defense specialist 1=competitor specialist 0/1/2/3/4/5
NC {mustBeInteger} = 1 % N- Si- Nutrient ratio 1=1 2=0.7 3=0.4 4=0.1
agg logical =1;
end

tic
%------------------------load aggLibrary----------------------
sLibname=loadaggLibrary;
%------------------------load parameters----------------------
p=parameters(n,NC);
% ------------------------initial value------------------------
N= 0*p.z+5; % mmol N/m3
S= 0*p.z+5; % mmol Si/m3
P= 0*p.z; % mmol N/m3
P= 0*p.z+0.001; % mmol N/m3

D1=0*p.z;  D2=0*p.z; D3=0*p.z; D4=0*p.z; D5=0*p.z; %Diatom     mmol N/m3
switch n
    case 0
D1=D1+0.001; % 
    case 1
D1=D1+0.001; % 
    case 2
D1=D1+0.001; D2=D2+0.001; %
    case 3
D1=D1+0.001; D2=D2+0.001; D3=D3+0.001; %
    case 4
D1=D1+0.001; D2=D2+0.001; D3=D3+0.001; D4=D4+0.001; %
    case 5 
D1=D1+0.001;  D2=D2+0.001; D3=D3+0.001; D4=D4+0.001; D5=D5+0.001; 
end

Dp= 0*p.z; %detritus of non-diatomphytoplankton  mmol N/m3
Dd= 0*p.z; %detritus of diatom mmol N/m3

%-----------------------------run----------------------------------
t=365*20;
tspan=0:t;
[p.diffmat,p.zmld,p.zmldgrid]=getDiffusivity(t,p);% get diffusity and mixed layer depth
y=[N,S,P,D1,D2,D3,D4,D5,Dp,Dd];

[t,y]=ode45(@derivativeseasonal, tspan, y, [],p);

%---------------------------Result output-------------------------
sim.t=t;
sim.NC=NC;
sim.p=p;
sim.y=y;
sim.N=y(1:length(t),1:p.xgrid);
sim.S=y(1:length(t),p.xgrid+1:2*p.xgrid);
sim.P=y(1:length(t),2*p.xgrid+1:3*p.xgrid);
sim.D1=y(1:length(t),3*p.xgrid+1:4*p.xgrid);
sim.D2=y(1:length(t),4*p.xgrid+1:5*p.xgrid);
sim.D3=y(1:length(t),5*p.xgrid+1:6*p.xgrid);
sim.D4=y(1:length(t),6*p.xgrid+1:7*p.xgrid);
sim.D5=y(1:length(t),7*p.xgrid+1:8*p.xgrid);
sim.Dp=y(1:length(t),8*p.xgrid+1:9*p.xgrid);
sim.Dd=y(1:length(t),9*p.xgrid+1:10*p.xgrid);

%plotdiatom(t,sim,p,season,NC)
%--------------------------derivative function----------------------
    function dydt = derivativeseasonal(t, y, p)

              i=2:p.xgrid;
              N=y(1:p.xgrid);
              S=y(p.xgrid+1:2*p.xgrid);
              P=y(2*p.xgrid+1:3*p.xgrid);
              D1=y(3*p.xgrid+1:4*p.xgrid);
              D2=y(4*p.xgrid+1:5*p.xgrid);
              D3=y(5*p.xgrid+1:6*p.xgrid);
              D4=y(6*p.xgrid+1:7*p.xgrid);
              D5=y(7*p.xgrid+1:8*p.xgrid);
              Dp=y(8*p.xgrid+1:9*p.xgrid);
              Dd=y(9*p.xgrid+1:10*p.xgrid);

     %diff N
                 Jdn(i) = -p.diffmat(:,floor(t)+1).*(N(i)-N(i-1))/p.dz;
                 Jdn(1) = 0;             % boundary top
                 Jdn(p.xgrid+1) = -p.diffdl*(p.NB-N(p.xgrid))/p.dz;  % boundary bottom
     %diff S
                 Jds(i) = -p.diffmat(:,floor(t)+1).*(S(i)-S(i-1))/p.dz;
                 Jds(1) = 0;             % boundary top
                 Jds(p.xgrid+1) = -p.diffdl*(p.SB-S(p.xgrid))/p.dz;  % boundary bottom

     %diff P   
              Jdp(i) = -p.diffmat(:,floor(t)+1).*(P(i)-P(i-1))/p.dz;
              Jdp(1) = 0;             % boundary top
              Jdp(p.xgrid+1) = 0;     % boundary bottom  
    
     %adv Diatom1 Defence specialist
                 Jad1(i)=p.vd(1)*D1(i-1);
                 Jad1(1) = 0;  % boundary top
                 Jad1(p.xgrid+1) = 0; % bottom p.vd(1)*D1(p.xgrid)
     %diff Diatom1   
              Jdd1(i) = -p.diffmat(:,floor(t)+1).*(D1(i)-D1(i-1))/p.dz;
              Jdd1(1) = 0;             % boundary top
              Jdd1(p.xgrid+1) = 0;     % boundary bottom
     %adv Diatom2
                 Jad2(i)=p.vd(2)*D2(i-1);
                 Jad2(1) = 0;  % boundary top
                 Jad2(p.xgrid+1) = 0; % bottom p.vd(2)*D2(p.xgrid)
     %diff Diatom2   
              Jdd2(i) = -p.diffmat(:,floor(t)+1).*(D2(i)-D2(i-1))/p.dz;
              Jdd2(1) = 0;             % boundary top
              Jdd2(p.xgrid+1) = 0;     % boundary bottom
     %adv Diatom3
                 Jad3(i)=p.vd(3)*D3(i-1);
                 Jad3(1) = 0;  % boundary top
                 Jad3(p.xgrid+1) = 0; % bottom p.vd(3)*D3(p.xgrid)
     %diff Diatom3   
              Jdd3(i) = -p.diffmat(:,floor(t)+1).*(D3(i)-D3(i-1))/p.dz;
              Jdd3(1) = 0;             % boundary top
              Jdd3(p.xgrid+1) = 0;     % boundary bottom 
     %adv Diatom4
                 Jad4(i)=p.vd(4)*D4(i-1);
                 Jad4(1) = 0;  % boundary top
                 Jad4(p.xgrid+1) = 0; % bottom p.vd(4)*D4(p.xgrid)
     %diff Diatom4   
              Jdd4(i) = -p.diffmat(:,floor(t)+1).*(D4(i)-D4(i-1))/p.dz;
              Jdd4(1) = 0;             % boundary top
              Jdd4(p.xgrid+1) = 0;     % boundary bottom
     %adv Diatom5 Competitor specialist
                 Jad5(i)=p.vd(5)*D5(i-1);
                 Jad5(1) = 0; % boundary top
                 Jad5(p.xgrid+1) = 0; % bottom p.vd(5)*D5(p.xgrid)
     %diff Diatom5   
              Jdd5(i) = -p.diffmat(:,floor(t)+1).*(D5(i)-D5(i-1))/p.dz;
              Jdd5(1) = 0;             % boundary top
              Jdd5(p.xgrid+1) = 0;     % boundary bottom

     %adv Dp
                 JaDp(i)=p.vDp*Dp(i-1);
                 JaDp(1) = 0;  % boundary top
                 JaDp(p.xgrid+1) = p.vDp*Dp(p.xgrid); % bottom
     %diff Dp
                 JdDp(i)= -p.diffmat(:,floor(t)+1).*(Dp(i)-Dp(i-1))/p.dz;
                 JdDp(1)= 0; %top
                 JdDp(p.xgrid+1)= 0;   %bottom
     %adv Dd
                 JaDd(i)=p.vDd*Dd(i-1);
                 JaDd(1) = 0;  % boundary top
                 JaDd(p.xgrid+1) = p.vDd*Dd(p.xgrid); % bottom
     %diff Dd
                 JdDd(i)= -p.diffmat(:,floor(t)+1).*(Dd(i)-Dd(i-1))/p.dz;
                 JdDd(1)= 0; %top
                 JdDd(p.xgrid+1)= 0;   %bottom

% aggregation flux update
if agg==1
Dmld=sum(D1(1:p.zmldgrid)+...
         D2(1:p.zmldgrid)+...
         D3(1:p.zmldgrid)+...
         D4(1:p.zmldgrid)+...
         D5(1:p.zmldgrid)); %Diatom Biomass in mixed layer        
              if Dmld>45
 [Jad1(i),Jad2(i),Jad3(i),Jad4(i),Jad5(i),JaDd(i)]=calllib(sLibname, 'f_getagg'...
     ,p.xgrid-1,D1(1:end-1),D2(1:end-1),D3(1:end-1),D4(1:end-1),D5(1:end-1),Dd(1:end-1));
JaDd(p.xgrid+1) = 200*Dd(p.xgrid); 
              end
end
      %assemble  
              JN = Jdn;
              JS = Jds;
              JP = Jdp;
              JD1= Jad1+Jdd1;
              JD2= Jad2+Jdd2;
              JD3= Jad3+Jdd3;
              JD4= Jad4+Jdd4;
              JD5= Jad5+Jdd5;
              JDp = JaDp+JdDp;
              JDd = JaDd+JdDd;

            dNdt = -(JN(2:p.xgrid+1)-JN(1:p.xgrid))/p.dz;
            dSdt = -(JS(2:p.xgrid+1)-JS(1:p.xgrid))/p.dz;
            dPdt = -(JP(2:p.xgrid+1)-JP(1:p.xgrid))/p.dz;
            dD1dt= -(JD1(2:p.xgrid+1)-JD1(1:p.xgrid))/p.dz;
            dD2dt= -(JD2(2:p.xgrid+1)-JD2(1:p.xgrid))/p.dz; 
            dD3dt= -(JD3(2:p.xgrid+1)-JD3(1:p.xgrid))/p.dz; 
            dD4dt= -(JD4(2:p.xgrid+1)-JD4(1:p.xgrid))/p.dz; 
            dD5dt= -(JD5(2:p.xgrid+1)-JD5(1:p.xgrid))/p.dz;
            dDpdt = -(JDp(2:p.xgrid+1)-JDp(1:p.xgrid))/p.dz;           
            dDddt = -(JDd(2:p.xgrid+1)-JDd(1:p.xgrid))/p.dz;

   %reaction

           % light
           I = getLight(P,D1,D2,D3,D4,D5,Dp,Dd,t,p);

           % Liebig's law of Growth rate
           g= p.pmax* min( I.*p.Al./(p.pmax+I.*p.Al) , N.*p.An./(p.pmax+N.*p.An) );

           trans1=[I.*p.Ald(1)./(p.pmaxd(1)+I.*p.Ald(1)),  N.*p.And(1)./(p.pmaxd(1)+N.*p.And(1)),  S.*p.Asd(1)./(p.pmaxd(1)+S.*p.Asd(1)) ];
           gd1= p.pmaxd(1)* min(trans1,[],2);
           trans2=[I.*p.Ald(2)./(p.pmaxd(2)+I.*p.Ald(2)),  N.*p.And(2)./(p.pmaxd(2)+N.*p.And(2)),  S.*p.Asd(2)./(p.pmaxd(2)+S.*p.Asd(2)) ];
           gd2= p.pmaxd(2)* min(trans2,[],2);
           trans3=[I.*p.Ald(3)./(p.pmaxd(3)+I.*p.Ald(3)),  N.*p.And(3)./(p.pmaxd(3)+N.*p.And(3)),  S.*p.Asd(3)./(p.pmaxd(3)+S.*p.Asd(3)) ];
           gd3= p.pmaxd(3)* min(trans3,[],2);
           trans4=[I.*p.Ald(4)./(p.pmaxd(4)+I.*p.Ald(4)),  N.*p.And(4)./(p.pmaxd(4)+N.*p.And(4)),  S.*p.Asd(4)./(p.pmaxd(4)+S.*p.Asd(4)) ];
           gd4= p.pmaxd(4)* min(trans4,[],2);
           trans5=[I.*p.Ald(5)./(p.pmaxd(5)+I.*p.Ald(5)),  N.*p.And(5)./(p.pmaxd(5)+N.*p.And(5)),  S.*p.Asd(5)./(p.pmaxd(5)+S.*p.Asd(5)) ];
           gd5= p.pmaxd(5)* min(trans5,[],2);
  %---------update---------
           dNdt= -g.*P -gd1.*D1 -gd2.*D2 -gd3.*D3 -gd4.*D4 -gd5.*D5 ...
               + p.tau.*(Dp+Dd)...
               + dNdt'; % -uptake+remineralization -J
           dSdt= -gd1.*D1*p.rhoSN -gd2.*D2*p.rhoSN -gd3.*D3*p.rhoSN -gd4.*D4*p.rhoSN -gd5.*D5*p.rhoSN...
               + p.tau*Dd*p.rhoSN...
               + dSdt';
           dPdt= g.*P- p.l*P - p.gamma*P.^2  +dPdt'; % growth-loss intrinsic-loss zoo-J
           dD1dt= gd1.*D1- p.l*D1 - p.gammad(1)*D1.^2  +dD1dt';
           dD2dt= gd2.*D2- p.l*D2 - p.gammad(2)*D2.^2  +dD2dt';
           dD3dt= gd3.*D3- p.l*D3 - p.gammad(3)*D3.^2  +dD3dt';
           dD4dt= gd4.*D4- p.l*D4 - p.gammad(4)*D4.^2  +dD4dt';
           dD5dt= gd5.*D5- p.l*D5 - p.gammad(5)*D5.^2  +dD5dt';
           dDpdt= p.l*P + p.gamma*P.^2 -p.tau.*Dp + dDpdt';
           dDddt= p.l*(D1+D2+D3+D4+D5) + p.gammad(1)*D1.^2 + p.gammad(2)*D2.^2+ p.gammad(3)*D3.^2 + p.gammad(4)*D4.^2 + p.gammad(5)*D5.^2 ...
                -p.tau*Dd...
                + dDddt';

           dydt= [dNdt;dSdt;dPdt;dD1dt;dD2dt;dD3dt;dD4dt;dD5dt;dDpdt;dDddt];
end

toc

end

