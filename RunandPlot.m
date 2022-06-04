%Directly run
%Only can change "agg": aggregation/no aggregation.   "season" does't work.
%All simulations are seasonal forcing.
function [sim1,sim2,sim3,sim4] = RunandPlot(season,agg)

arguments
season logical =1; % true=seasonal forcing  false=constant
agg logical =1;
end

%---------------------Seasonal Run---------------------------
switch agg
    case 1
 sim1=runDiatom(5,1,1);%Si:N 1
 sim2=runDiatom(5,2,1);%Si:N 0.7
 sim3=runDiatom(5,3,1);%Si:N 0.4
 sim4=runDiatom(5,4,1);%Si:N 0.1
% sim1=sim3;
% sim2=sim3;
% sim3=sim1;
% sim4=sim1;
    case 0
sim1=runDiatom(5,1,0);%Si:N 1
sim2=runDiatom(5,2,0);%Si:N 0.7
sim3=runDiatom(5,3,0);%Si:N 0.4
sim4=runDiatom(5,4,0);%Si:N 0.1
end
%
%-------------------------------Plot------------------------------
%
switch season 
    case 0 %------------------Constant------------------
% Constant Plot
figure(1)
t1=tiledlayout(2,2);
nexttile %Non-limited each
    plot(sim1.P(end,:),-sim1.p.z,LineWidth=2)%P
    hold on 
    plot(sim1.D1(end,:),-sim1.p.z,LineWidth=2)%D
    plot(sim1.D2(end,:),-sim1.p.z,LineWidth=2)
    plot(sim1.D3(end,:),-sim1.p.z,LineWidth=2)
    plot(sim1.D4(end,:),-sim1.p.z,LineWidth=2)
    plot(sim1.D5(end,:),-sim1.p.z,LineWidth=2)
    hold off
  legend('P','D1 (Defence specialist)','D2','D3','D4','D5 (Competitor specialist)','Location','southeast','FontSize',12)
  legend boxoff
  ylabel('Depth (m)')
   switch sim1.NC  %simulate spring
        case 1
          title('Non limited')
        case 2 % simulate summer
          title('Nutrient limited')
          xlim([0 0.8])
   end
   set(gca,FontSize=20)

nexttile%Nutrient limited each
    plot(sim2.P(end,:),-sim2.p.z,LineWidth=2)%P
    hold on 
    plot(sim2.D1(end,:),-sim2.p.z,LineWidth=2)%D
    plot(sim2.D2(end,:),-sim2.p.z,LineWidth=2)
    plot(sim2.D3(end,:),-sim2.p.z,LineWidth=2)
    plot(sim2.D4(end,:),-sim2.p.z,LineWidth=2)
    plot(sim2.D5(end,:),-sim2.p.z,LineWidth=2)
    hold off
  legend('P','D1 (Defence specialist)','D2','D3','D4','D5 (Competitor specialist)','Location','southeast','FontSize',12)
  legend boxoff
   switch sim2.NC  %simulate spring
        case 1
          title('Non limited')
        case 2 % simulate summer
          title('Nutrient limited')
          xlim([0 0.8])
   end
   set(gca,FontSize=20)

nexttile %Non-limited P&Dia&De
  plot(sim1.P(end,:),-sim1.p.z,LineWidth=2)
  hold on
  plot((sim1.D1(end,:)+sim1.D2(end,:)+sim1.D3(end,:)+sim1.D4(end,:)+sim1.D5(end,:)),-sim1.p.z,LineWidth=2)
  plot(sim1.Dp(end,:),-sim1.p.z,LineWidth=1.5)
  plot(sim1.Dd(end,:),-sim1.p.z,LineWidth=1.5)
  hold off
  legend('P','D','Dp','Dd','FontSize',10,'Location','southeast')
  legend boxoff
  xlabel('Biomass (mmol N m^-^3)')
  ylabel('Depth (m)')
  set(gca,FontSize=20)
  if sim1.NC==2
  xlim([0 2])
  end

nexttile%Nutrient limited P&Dia&De
  plot(sim2.P(end,:),-sim2.p.z,LineWidth=2)
  hold on
  plot((sim2.D1(end,:)+sim2.D2(end,:)+sim2.D3(end,:)+sim2.D4(end,:)+sim2.D5(end,:)),-sim2.p.z,LineWidth=2)
  plot(sim2.Dp(end,:),-sim2.p.z,LineWidth=1.5)
  plot(sim2.Dd(end,:),-sim2.p.z,LineWidth=1.5)
  hold off
  legend('P','D','Dp','Dd','FontSize',10,'Location','southeast')
  legend boxoff
  xlabel('Biomass (mmol N m^-^3)')
  set(gca,FontSize=20)
  if sim2.NC==2
  xlim([0 2.5])
  end

figure(2)
tl1=tiledlayout(1,1);
ax1=axes(tl1);
 plot(ax1,sim1.N(end,:),-sim1.p.z,'b',LineWidth=2)
hold on
 plot(ax1,sim1.S(end,:),-sim1.p.z,'k',LineWidth=2)
 hold off
 legend('N','Si','FontSize',8,'Color','None','Box','off')
 xlabel('Concentration (mmol m^-^3)')
ylabel('Depth (m)')
set(gca,FontSize=15)
 ax2=axes(tl1);
 plot(ax2,getLight(sim1.P(end,:)',sim1.D1(end,:)',sim1.D2(end,:)',sim1.D3(end,:)',sim1.D5(end,:)',sim1.D5(end,:)',sim1.Dp(end,:)',sim1.Dd(end,:)',sim1.t(end),sim1.p,sim1.season),-sim1.p.z,'r',LineWidth=1.5)
 ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
xlabel('Light intensity (W m^{-2})')
legend('L','FontSize',8,'Location','southeast','Color','None','Box','off')
set(gca,FontSize=15)

figure(3)
tl1=tiledlayout(1,1);
ax1=axes(tl1);
 plot(ax1,sim2.N(end,:),-sim2.p.z,'b',LineWidth=2)
hold on
 plot(ax1,sim2.S(end,:),-sim2.p.z,'k',LineWidth=2)
 hold off
 legend('N','Si','FontSize',8,'Color','None','Box','off')
 xlabel('Concentration (mmol m^-^3)')
ylabel('Depth (m)')
set(gca,FontSize=15)
 ax2=axes(tl1);
 plot(ax2,getLight(sim2.P(end,:)',sim2.D1(end,:)',sim2.D2(end,:)',sim2.D3(end,:)',sim2.D4(end,:)',sim2.D5(end,:)',sim2.Dp(end,:)',sim2.Dd(end,:)',sim2.t(end),sim2.p,sim2.season),-sim2.p.z,'r',LineWidth=1.5)
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
xlabel('Light intensity (W m^{-2})')
legend('L','FontSize',8,'Location','southeast','Color','None','Box','off')
set(gca,FontSize=15)

    case 1 %--------------------Seasonal scenario-----------------------
%% ------------------Calculation----------------
%water column matrix for surface plot
all1=(sim1.P+sim1.D1+sim1.D2+sim1.D3+sim1.D4+sim1.D5);%Si:N 1
all2=(sim2.P+sim2.D1+sim2.D2+sim2.D3+sim2.D4+sim2.D5);%Si:N 0.7
all3=(sim3.P+sim3.D1+sim3.D2+sim3.D3+sim3.D4+sim3.D5);%Si:N 0.4
all4=(sim4.P+sim4.D1+sim4.D2+sim4.D3+sim4.D4+sim4.D5);%Si:N 0.1
%time calculation
spts=(length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12));%spring time span
spday=length((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)));%spring day number
suts=(length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12));%summer time span
suday=length((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)));%summer day number
%biomass calculation
% ALL
%whole last year
whall1=  sum(sim1.P((length(sim1.t)-365):length(sim1.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim1.D1((length(sim1.t)-365):length(sim1.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim1.D2((length(sim1.t)-365):length(sim1.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim1.D3((length(sim1.t)-365):length(sim1.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim1.D4((length(sim1.t)-365):length(sim1.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim1.D5((length(sim1.t)-365):length(sim1.t),1:sim1.p.xgrid*2/3),2);
whall2=  sum(sim2.P((length(sim2.t)-365):length(sim2.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim2.D1((length(sim2.t)-365):length(sim2.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim2.D2((length(sim2.t)-365):length(sim2.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim2.D3((length(sim2.t)-365):length(sim2.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim2.D4((length(sim2.t)-365):length(sim2.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim2.D5((length(sim2.t)-365):length(sim2.t),1:sim1.p.xgrid*2/3),2);
whall3=  sum(sim3.P((length(sim3.t)-365):length(sim3.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim3.D1((length(sim3.t)-365):length(sim3.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim3.D2((length(sim3.t)-365):length(sim3.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim3.D3((length(sim3.t)-365):length(sim3.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim3.D4((length(sim3.t)-365):length(sim3.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim3.D5((length(sim3.t)-365):length(sim3.t),1:sim1.p.xgrid*2/3),2);
whall4=  sum(sim4.P((length(sim4.t)-365):length(sim4.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim4.D1((length(sim4.t)-365):length(sim4.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim4.D2((length(sim4.t)-365):length(sim4.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim4.D3((length(sim4.t)-365):length(sim4.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim4.D4((length(sim4.t)-365):length(sim4.t),1:sim1.p.xgrid*2/3),2) +...
       sum(sim4.D5((length(sim4.t)-365):length(sim4.t),1:sim1.p.xgrid*2/3),2);
%spring bloom
spall1=sum(sum(sim1.P(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D1(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D2(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D3(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D4(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D5(spts,1:sim1.p.xgrid*2/3),2));
spall1=spall1/spday;
spall2=sum(sum(sim2.P(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D1(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D2(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D3(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D4(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D5(spts,1:sim1.p.xgrid*2/3),2));
spall2=spall2/spday;
spall3=sum(sum(sim3.P(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D1(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D2(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D3(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D4(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D5(spts,1:sim1.p.xgrid*2/3),2));
spall3=spall3/spday;
spall4=sum(sum(sim4.P(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D1(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D2(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D3(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D4(spts,1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D5(spts,1:sim1.p.xgrid*2/3),2));
spall4=spall4/spday;
%summer
suall1=sum(sum(sim1.P((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D1((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D2((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D3((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D4((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D5((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suall1=suall1/length((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)));
suall2=sum(sum(sim2.P((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D1((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D2((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D3((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D4((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D5((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suall2=suall2/length((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)));
suall3=sum(sum(sim3.P((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D1((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D2((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D3((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D4((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D5((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suall3=suall3/length((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)));
suall4=sum(sum(sim4.P((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D1((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D2((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D3((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D4((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D5((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suall4=suall4/length((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)));
%non-d phyto
%whole last year
whphyto1=  sum(sim1.P((length(sim1.t)-365):length(sim1.t),1:sim1.p.xgrid*2/3),2);
whphyto2=  sum(sim2.P((length(sim2.t)-365):length(sim2.t),1:sim1.p.xgrid*2/3),2);
whphyto3=  sum(sim3.P((length(sim3.t)-365):length(sim3.t),1:sim1.p.xgrid*2/3),2);
whphyto4=  sum(sim4.P((length(sim4.t)-365):length(sim4.t),1:sim1.p.xgrid*2/3),2);
%spring bloom
spphyto1=sum(sum(sim1.P((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spphyto1=spphyto1/length((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)));
spphyto1=spphyto1/spall1;%ratio
spphyto2=sum(sum(sim2.P((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spphyto2=spphyto2/length((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)));
spphyto2=spphyto2/spall2;
spphyto3=sum(sum(sim3.P((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spphyto3=spphyto3/length((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)));
spphyto3=spphyto3/spall3;
spphyto4=sum(sum(sim4.P((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spphyto4=spphyto4/length((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)));
spphyto4=spphyto4/spall4;
%summer
suphyto1=sum(sum(sim1.P((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suphyto1=suphyto1/length((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)));
suphyto1=suphyto1/suall1;
suphyto2=sum(sum(sim2.P((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suphyto2=suphyto2/length((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)));
suphyto2=suphyto2/suall2;
suphyto3=sum(sum(sim3.P((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suphyto3=suphyto3/length((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)));
suphyto3=suphyto3/suall3;
suphyto4=sum(sum(sim4.P((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suphyto4=suphyto4/length((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)));
suphyto4=suphyto4/suall4;
% alldiatom
%whole last year
whalldia1= sum(sim1.D1((length(sim1.t)-365):length(sim1.t),1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D2((length(sim1.t)-365):length(sim1.t),1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D3((length(sim1.t)-365):length(sim1.t),1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D4((length(sim1.t)-365):length(sim1.t),1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D5((length(sim1.t)-365):length(sim1.t),1:sim1.p.xgrid*2/3),2);

whalldia2= sum(sim2.D1((length(sim2.t)-365):length(sim2.t),1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D2((length(sim2.t)-365):length(sim2.t),1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D3((length(sim2.t)-365):length(sim2.t),1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D4((length(sim2.t)-365):length(sim2.t),1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D5((length(sim2.t)-365):length(sim2.t),1:sim1.p.xgrid*2/3),2);
whalldia3= sum(sim3.D1((length(sim3.t)-365):length(sim3.t),1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D2((length(sim3.t)-365):length(sim3.t),1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D3((length(sim3.t)-365):length(sim3.t),1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D4((length(sim3.t)-365):length(sim3.t),1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D5((length(sim3.t)-365):length(sim3.t),1:sim1.p.xgrid*2/3),2);
whalldia4= sum(sim4.D1((length(sim4.t)-365):length(sim4.t),1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D2((length(sim4.t)-365):length(sim4.t),1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D3((length(sim4.t)-365):length(sim4.t),1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D4((length(sim4.t)-365):length(sim4.t),1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D5((length(sim4.t)-365):length(sim4.t),1:sim1.p.xgrid*2/3),2);
%spring bloom
spalldia1=sum(sum(sim1.D1((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D2((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D3((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D4((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D5((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spalldia1=spalldia1/length((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)));
spalldia1=spalldia1/spall1;%ratio
spalldia2=sum(sum(sim2.D1((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D2((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D3((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D4((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D5((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spalldia2=spalldia2/length((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)));
spalldia2=spalldia2/spall2;
spalldia3=sum(sum(sim3.D1((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D2((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D3((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D4((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D5((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spalldia3=spalldia3/length((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)));
spalldia3=spalldia3/spall3;
spalldia4=sum(sum(sim4.D1((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D2((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D3((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D4((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D5((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spalldia4=spalldia4/length((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)));
spalldia4=spalldia4/spall4;
%summer
sualldia1=sum(sum(sim1.D1((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D2((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D3((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D4((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim1.D5((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
sualldia1=sualldia1/length((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)));
sualldia1=sualldia1/suall1;%ratio
sualldia2=sum(sum(sim2.D1((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D2((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D3((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D4((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim2.D5((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
sualldia2=sualldia2/length((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)));
sualldia2=sualldia2/suall2;
sualldia3=sum(sum(sim3.D1((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D2((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D3((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D4((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim3.D5((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
sualldia3=sualldia3/length((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)));
sualldia3=sualldia3/suall3;
sualldia4=sum(sum(sim4.D1((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D2((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D3((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D4((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2) +...
           sum(sim4.D5((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
sualldia4=sualldia4/length((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)));
sualldia4=sualldia4/suall4;
%D1-D5
%D1 defence specialist
%whole last year
whD11=sum(sim1.D1((length(sim1.t)-365):length(sim1.t),1:sim1.p.xgrid*2/3),2);
whD12=sum(sim2.D1((length(sim2.t)-365):length(sim2.t),1:sim1.p.xgrid*2/3),2);
whD13=sum(sim3.D1((length(sim3.t)-365):length(sim3.t),1:sim1.p.xgrid*2/3),2);
whD14=sum(sim4.D1((length(sim4.t)-365):length(sim4.t),1:sim1.p.xgrid*2/3),2);
%spring bloom
spD11=sum(sum(sim1.D1((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD11=spD11/length((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)));
spD11=spD11/spall1;
spD12=sum(sum(sim2.D1((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD12=spD12/length((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)));
spD12=spD12/spall2;
spD13=sum(sum(sim3.D1((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD13=spD13/length((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)));
spD13=spD13/spall3;
spD14=sum(sum(sim4.D1((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD14=spD14/length((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)));
spD14=spD14/spall4;
%summer
suD11=sum(sum(sim1.D1((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD11=suD11/length((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)));
suD11=suD11/suall1;
suD12=sum(sum(sim2.D1((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD12=suD12/length((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)));
suD12=suD12/suall2;
suD13=sum(sum(sim3.D1((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD13=suD13/length((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)));
suD13=suD13/suall3;
suD14=sum(sum(sim4.D1((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD14=suD14/length((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)));
suD14=suD14/suall4;
%D2
%whole last year
whD21=sum(sim1.D2((length(sim1.t)-365):length(sim1.t),1:sim1.p.xgrid*2/3),2);
whD22=sum(sim2.D2((length(sim2.t)-365):length(sim2.t),1:sim1.p.xgrid*2/3),2);
whD23=sum(sim3.D2((length(sim3.t)-365):length(sim3.t),1:sim1.p.xgrid*2/3),2);
whD24=sum(sim4.D2((length(sim4.t)-365):length(sim4.t),1:sim1.p.xgrid*2/3),2);
%spring bloom
spD21=sum(sum(sim1.D2((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD21=spD21/length((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)));
spD21=spD21/spall1;
spD22=sum(sum(sim2.D2((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD22=spD22/length((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)));
spD22=spD22/spall2;
spD23=sum(sum(sim3.D2((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD23=spD23/length((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)));
spD23=spD23/spall3;
spD24=sum(sum(sim4.D2((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD24=spD24/length((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)));
spD24=spD24/spall4;
%summer
suD21=sum(sum(sim1.D2((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD21=suD21/length((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)));
suD21=suD21/suall1;
suD22=sum(sum(sim2.D2((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD22=suD22/length((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)));
suD22=suD22/suall2;
suD23=sum(sum(sim3.D2((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD23=suD23/length((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)));
suD23=suD23/suall3;
suD24=sum(sum(sim4.D2((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD24=suD24/length((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)));
suD24=suD24/suall4;
%D3
%whole last year
whD31=sum(sim1.D3((length(sim1.t)-365):length(sim1.t),1:sim1.p.xgrid*2/3),2);
whD32=sum(sim2.D3((length(sim2.t)-365):length(sim2.t),1:sim1.p.xgrid*2/3),2);
whD33=sum(sim3.D3((length(sim3.t)-365):length(sim3.t),1:sim1.p.xgrid*2/3),2);
whD34=sum(sim4.D3((length(sim4.t)-365):length(sim4.t),1:sim1.p.xgrid*2/3),2);
%spring bloom
spD31=sum(sum(sim1.D3((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD31=spD31/length((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)));
spD31=spD31/spall1;
spD32=sum(sum(sim2.D3((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD32=spD32/length((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)));
spD32=spD32/spall2;
spD33=sum(sum(sim3.D3((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD33=spD33/length((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)));
spD33=spD33/spall3;
spD34=sum(sum(sim4.D3((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD34=spD34/length((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)));
spD34=spD34/spall4;
%summer
suD31=sum(sum(sim1.D3((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD31=suD31/length((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)));
suD31=suD31/suall1;
suD32=sum(sum(sim2.D3((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD32=suD32/length((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)));
suD32=suD32/suall2;
suD33=sum(sum(sim3.D3((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD33=suD33/length((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)));
suD33=suD33/suall3;
suD34=sum(sum(sim4.D3((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD34=suD34/length((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)));
suD34=suD34/suall4;
%D4
%whole last year
whD41=sum(sim1.D4((length(sim1.t)-365):length(sim1.t),1:sim1.p.xgrid*2/3),2);
whD42=sum(sim2.D4((length(sim2.t)-365):length(sim2.t),1:sim1.p.xgrid*2/3),2);
whD43=sum(sim3.D4((length(sim3.t)-365):length(sim3.t),1:sim1.p.xgrid*2/3),2);
whD44=sum(sim4.D4((length(sim4.t)-365):length(sim4.t),1:sim1.p.xgrid*2/3),2);
%spring bloom
spD41=sum(sum(sim1.D4((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD41=spD41/length((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)));
spD41=spD41/spall1;
spD42=sum(sum(sim2.D4((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD42=spD42/length((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)));
spD42=spD42/spall2;
spD43=sum(sum(sim3.D4((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD43=spD43/length((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)));
spD43=spD43/spall3;
spD44=sum(sum(sim4.D4((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD44=spD44/length((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)));
spD44=spD44/spall4;
%summer
suD41=sum(sum(sim1.D4((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD41=suD41/length((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)));
suD41=suD41/suall1;
suD42=sum(sum(sim2.D4((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD42=suD42/length((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)));
suD42=suD42/suall2;
suD43=sum(sum(sim3.D4((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD43=suD43/length((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)));
suD43=suD43/suall3;
suD44=sum(sum(sim4.D4((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD44=suD44/length((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)));
suD44=suD44/suall4;
%D5 competitor specialist
%whole last year
whD51=sum(sim1.D5((length(sim1.t)-365):length(sim1.t),1:sim1.p.xgrid*2/3),2);
whD52=sum(sim2.D5((length(sim2.t)-365):length(sim2.t),1:sim1.p.xgrid*2/3),2);
whD53=sum(sim3.D5((length(sim3.t)-365):length(sim3.t),1:sim1.p.xgrid*2/3),2);
whD54=sum(sim4.D5((length(sim4.t)-365):length(sim4.t),1:sim1.p.xgrid*2/3),2);
%spring bloom
spD51=sum(sum(sim1.D5((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD51=spD51/length((length(sim1.t)-365+floor(365/12)):(length(sim1.t)-365+floor(4*365/12)));
spD51=spD51/spall1;
spD52=sum(sum(sim2.D5((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD52=spD52/length((length(sim2.t)-365+floor(365/12)):(length(sim2.t)-365+floor(4*365/12)));
spD52=spD52/spall2;
spD53=sum(sum(sim3.D5((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD53=spD53/length((length(sim3.t)-365+floor(365/12)):(length(sim3.t)-365+floor(4*365/12)));
spD53=spD53/spall3;
spD54=sum(sum(sim4.D5((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)),1:sim1.p.xgrid*2/3),2));
spD54=spD54/length((length(sim4.t)-365+floor(365/12)):(length(sim4.t)-365+floor(4*365/12)));
spD54=spD54/spall4;
%summer
suD51=sum(sum(sim1.D5((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD51=suD51/length((length(sim1.t)-365+floor(5*365/12)):(length(sim1.t)-365+floor(7*365/12)));
suD51=suD51/suall1;
suD52=sum(sum(sim2.D5((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD52=suD52/length((length(sim2.t)-365+floor(5*365/12)):(length(sim2.t)-365+floor(7*365/12)));
suD52=suD52/suall2;
suD53=sum(sum(sim3.D5((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD53=suD53/length((length(sim3.t)-365+floor(5*365/12)):(length(sim3.t)-365+floor(7*365/12)));
suD53=suD53/suall3;
suD54=sum(sum(sim4.D5((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)),1:sim1.p.xgrid*2/3),2));
suD54=suD54/length((length(sim4.t)-365+floor(5*365/12)):(length(sim4.t)-365+floor(7*365/12)));
suD54=suD54/suall4;
%% -------------------------------------Plot-----------------------------------
%-------------------------------Surface&Profile----------------------------------
figure(1)
t1=tiledlayout(4,3);
%Si:N 1
nexttile
surface(0:365,-sim1.p.z,all1((length(sim1.t)-365):length(sim1.t),:)')
   xlim([0 365])
   shading interp
   bar=colorbar;
   bar.Label.String = 'Biomass (mmol N m^-^3)';
   xticks([0 365/12 2*365/12 3*365/12 4*365/12 5*365/12 6*365/12 7*365/12 8*365/12 9*365/12 10*365/12 11*365/12 365]);
   xticklabels([])
   ylabel({'Si:N=1';'Depth (m)'})
   set(gca,FontSize=10)
nexttile%spring profile day106
plot(sim1.P(length(sim1.t)-260,:),-sim1.p.z,LineWidth=1)% day 106
    hold on 
    plot(sim1.D1(length(sim1.t)-260,:),-sim1.p.z,LineWidth=1)%
    plot(sim1.D2(length(sim1.t)-260,:),-sim1.p.z,LineWidth=1)
    plot(sim1.D3(length(sim1.t)-260,:),-sim1.p.z,LineWidth=1)
    plot(sim1.D4(length(sim1.t)-260,:),-sim1.p.z,LineWidth=1)
    plot(sim1.D5(length(sim1.t)-260,:),-sim1.p.z,LineWidth=1)
    hold off
    title('Spring (day 106)')
    set(gca,FontSize=10)
nexttile %summer profile day181
plot(sim1.P(length(sim1.t)-185,:),-sim1.p.z,LineWidth=1)%
    hold on 
    plot(sim1.D1(length(sim1.t)-185,:),-sim1.p.z,LineWidth=1)%
    plot(sim1.D2(length(sim1.t)-185,:),-sim1.p.z,LineWidth=1)
    plot(sim1.D3(length(sim1.t)-185,:),-sim1.p.z,LineWidth=1)
    plot(sim1.D4(length(sim1.t)-185,:),-sim1.p.z,LineWidth=1)
    plot(sim1.D5(length(sim1.t)-185,:),-sim1.p.z,LineWidth=1)
    hold off
    title('Summer (day 181)')
    set(gca,FontSize=10)

%Si:N 0.7
nexttile
surface(0:365,-sim2.p.z,all2((length(sim2.t)-365):length(sim2.t),:)')
   xlim([0 365])
   shading interp
   bar=colorbar;
   bar.Label.String = 'Biomass (mmol N m^-^3)';
   xticks([0 365/12 2*365/12 3*365/12 4*365/12 5*365/12 6*365/12 7*365/12 8*365/12 9*365/12 10*365/12 11*365/12 365]);
   xticklabels([])
   ylabel({'Si:N=0.7';'Depth (m)'})
   set(gca,FontSize=10)
nexttile%spring profile day106
plot(sim2.P(length(sim2.t)-260,:),-sim2.p.z,LineWidth=1)% day 106
    hold on 
    plot(sim2.D1(length(sim2.t)-260,:),-sim2.p.z,LineWidth=1)%
    plot(sim2.D2(length(sim2.t)-260,:),-sim2.p.z,LineWidth=1)
    plot(sim2.D3(length(sim2.t)-260,:),-sim2.p.z,LineWidth=1)
    plot(sim2.D4(length(sim2.t)-260,:),-sim2.p.z,LineWidth=1)
    plot(sim2.D5(length(sim2.t)-260,:),-sim2.p.z,LineWidth=1)
    hold off
    set(gca,FontSize=10)
nexttile %summer profile day181
plot(sim2.P(length(sim2.t)-185,:),-sim2.p.z,LineWidth=1)%
    hold on 
    plot(sim2.D1(length(sim2.t)-185,:),-sim2.p.z,LineWidth=1)%
    plot(sim2.D2(length(sim2.t)-185,:),-sim2.p.z,LineWidth=1)
    plot(sim2.D3(length(sim2.t)-185,:),-sim2.p.z,LineWidth=1)
    plot(sim2.D4(length(sim2.t)-185,:),-sim2.p.z,LineWidth=1)
    plot(sim2.D5(length(sim2.t)-185,:),-sim2.p.z,LineWidth=1)
    hold off
  set(gca,FontSize=10)

%Si:N 0.4
nexttile
surface(0:365,-sim3.p.z,all3((length(sim3.t)-365):length(sim3.t),:)')
  xlim([0 365])
  shading interp
  bar=colorbar;
  bar.Label.String = 'Biomass (mmol N m^-^3)';
  xticks([0 365/12 2*365/12 3*365/12 4*365/12 5*365/12 6*365/12 7*365/12 8*365/12 9*365/12 10*365/12 11*365/12 365]);
  xticklabels([])
  ylabel({'Si:N=0.4';'Depth (m)'})
  set(gca,FontSize=10)
nexttile%spring profile day106
plot(sim3.P(length(sim3.t)-260,:),-sim3.p.z,LineWidth=1)% day 106
    hold on 
    plot(sim3.D1(length(sim3.t)-260,:),-sim3.p.z,LineWidth=1)%
    plot(sim3.D2(length(sim3.t)-260,:),-sim3.p.z,LineWidth=1)
    plot(sim3.D3(length(sim3.t)-260,:),-sim3.p.z,LineWidth=1)
    plot(sim3.D4(length(sim3.t)-260,:),-sim3.p.z,LineWidth=1)
    plot(sim3.D5(length(sim3.t)-260,:),-sim3.p.z,LineWidth=1)
    hold off
    set(gca,FontSize=10)
nexttile %summer profile day181
plot(sim3.P(length(sim3.t)-185,:),-sim3.p.z,LineWidth=1)%
    hold on 
    plot(sim3.D1(length(sim3.t)-185,:),-sim3.p.z,LineWidth=1)%
    plot(sim3.D2(length(sim3.t)-185,:),-sim3.p.z,LineWidth=1)
    plot(sim3.D3(length(sim3.t)-185,:),-sim3.p.z,LineWidth=1)
    plot(sim3.D4(length(sim3.t)-185,:),-sim3.p.z,LineWidth=1)
    plot(sim3.D5(length(sim3.t)-185,:),-sim3.p.z,LineWidth=1)
    hold off
    set(gca,FontSize=10)

%Si:N 0.1
nexttile
  surface(0:365,-sim4.p.z,all4((length(sim4.t)-365):length(sim4.t),:)')
  xlim([0 365])
  shading interp
  bar=colorbar;
  bar.Label.String = 'Biomass (mmol N m^-^3)';
  xticks([0 365/12 2*365/12 3*365/12 4*365/12 5*365/12 6*365/12 7*365/12 8*365/12 9*365/12 10*365/12 11*365/12 365]);
  xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
  xlabel('Month')
  ylabel({'Si:N=0.1';'Depth (m)'})
  set(gca,FontSize=10)
nexttile%spring profile day106
    plot(sim4.P(length(sim4.t)-260,:),-sim4.p.z,LineWidth=1)% day 106
    hold on 
    plot(sim4.D1(length(sim4.t)-260,:),-sim4.p.z,LineWidth=1)%
    plot(sim4.D2(length(sim4.t)-260,:),-sim4.p.z,LineWidth=1)
    plot(sim4.D3(length(sim4.t)-260,:),-sim4.p.z,LineWidth=1)
    plot(sim4.D4(length(sim4.t)-260,:),-sim4.p.z,LineWidth=1)
    plot(sim4.D5(length(sim4.t)-260,:),-sim4.p.z,LineWidth=1)
    hold off
    xlabel('Biomass (mmol N m^-^3)')
    set(gca,FontSize=10)
nexttile %summer profile day181
plot(sim4.P(length(sim4.t)-185,:),-sim4.p.z,LineWidth=1)%
    hold on 
    plot(sim4.D1(length(sim4.t)-185,:),-sim4.p.z,LineWidth=1)%
    plot(sim4.D2(length(sim4.t)-185,:),-sim4.p.z,LineWidth=1)
    plot(sim4.D3(length(sim4.t)-185,:),-sim4.p.z,LineWidth=1)
    plot(sim4.D4(length(sim4.t)-185,:),-sim4.p.z,LineWidth=1)
    plot(sim4.D5(length(sim4.t)-185,:),-sim4.p.z,LineWidth=1)
    hold off
  legend('P','D1','D2','D3','D4','D5','Location','southeast')
  legend boxoff
  xlabel('Biomass (mmol N m^-^3)')
  set(gca,FontSize=10)

%--------------------------------TiledPlot----------------------------
RSN=[1 0.7 0.4 0.1]; %Si:N ratio
figure(2)
t2=tiledlayout(8,3);
%All
nexttile
  plot(0:365,whall1,Linewidth=2)%,Linewidth=2
  hold on
  plot(0:365,whall2,Linewidth=2);
  plot(0:365,whall3,Linewidth=2);
  plot(0:365,whall4,Linewidth=2);
  yl=ylim;
  shadex=[365/12 5*365/12;4*365/12 7*365/12;4*365/12 7*365/12;365/12 5*365/12];
  shadey=[yl(1) yl(1);      yl(1) yl(1);      yl(2) yl(2);      yl(2) yl(2)];
  c=[211,211,211]/255;
  fill(shadex,shadey,c, 'linestyle','none','FaceAlpha','0.3')
  hold off
  xlim([0 365])
  xticks([0 365/12 2*365/12 3*365/12 4*365/12 5*365/12 6*365/12 7*365/12 8*365/12 9*365/12 10*365/12 11*365/12 365]);
  xticklabels([])
  ylabel({'Total';'biomass';' mmol N m^{-3}'})
  set(gca,FontSize=10)
nexttile
plot(RSN,[spall1 spall2 spall3 spall4],'-o',LineWidth=2,MarkerSize=4,MarkerFaceColor='#0072BD')
xticks(flip(RSN))
set(gca,'XDir','reverse',FontSize=10);
xticklabels([])
title('Spring')
nexttile
  plot(RSN,[suall1 suall2 suall3 suall4],'-o',LineWidth=2,MarkerSize=4,MarkerFaceColor='#0072BD')
  xticks(flip(RSN))
  xticklabels([])
  set(gca,'XDir','reverse',FontSize=10);
  title('Summer')

%non-dia phyto
nexttile
  plot(0:365,whphyto1,Linewidth=2)%,Linewidth=2
  hold on
  plot(0:365,whphyto2,Linewidth=2);
  plot(0:365,whphyto3,Linewidth=2);
  plot(0:365,whphyto4,Linewidth=2);
  yl=ylim;
  shadex=[365/12 5*365/12;4*365/12 7*365/12;4*365/12 7*365/12;365/12 5*365/12];
  shadey=[yl(1) yl(1);      yl(1) yl(1);      yl(2) yl(2);      yl(2) yl(2)];
  c=[211,211,211]/255;
  fill(shadex,shadey,c, 'linestyle','none','FaceAlpha','0.3')
  hold off
  xlim([0 365])
  xticks([0 365/12 2*365/12 3*365/12 4*365/12 5*365/12 6*365/12 7*365/12 8*365/12 9*365/12 10*365/12 11*365/12 365]);
  xticklabels([])
  ylabel({'Non-diatom';'phytoplankton';' mmol N m^{-3}'})
set(gca,FontSize=10)
nexttile
plot(RSN,[spphyto1 spphyto2 spphyto3 spphyto4],'-o',LineWidth=2,MarkerSize=4,MarkerFaceColor='#0072BD')
xticks(flip(RSN))
set(gca,'XDir','reverse',FontSize=10);
ylabel({'Relative';'biomass'})
xticklabels([])
  nexttile
  plot(RSN,[suphyto1 suphyto2 suphyto3 suphyto4],'-o',LineWidth=2,MarkerSize=4,MarkerFaceColor='#0072BD')
  xticks(flip(RSN))
  set(gca,'XDir','reverse',FontSize=10);
  xticklabels([])

%all diatom
nexttile
  plot(0:365,whalldia1,Linewidth=2)%,Linewidth=2
  hold on
  plot(0:365,whalldia2,Linewidth=2);
  plot(0:365,whalldia3,Linewidth=2);
  plot(0:365,whalldia4,Linewidth=2);
  yl=ylim;
  shadex=[365/12 5*365/12;4*365/12 7*365/12;4*365/12 7*365/12;365/12 5*365/12];
  shadey=[yl(1) yl(1);      yl(1) yl(1);      yl(2) yl(2);      yl(2) yl(2)];
  c=[211,211,211]/255;
  fill(shadex,shadey,c, 'linestyle','none','FaceAlpha','0.3')
  hold off
  xlim([0 365])
  xticks([0 365/12 2*365/12 3*365/12 4*365/12 5*365/12 6*365/12 7*365/12 8*365/12 9*365/12 10*365/12 11*365/12 365]);
  xticklabels([])
  ylabel({'All diatom';' mmol N m^{-3}'})
  set(gca,FontSize=10)
nexttile
plot(RSN,[spalldia1 spalldia2 spalldia3 spalldia4],'-o',LineWidth=2,MarkerSize=4,MarkerFaceColor='#0072BD')
xticks(flip(RSN))
set(gca,'XDir','reverse',FontSize=10);
ylabel({'Relative';'biomass'})
xticklabels([])
  nexttile
  plot(RSN,[sualldia1 sualldia2 sualldia3 sualldia4],'-o',LineWidth=2,MarkerSize=4,MarkerFaceColor='#0072BD')
  xticks(flip(RSN))
  set(gca,'XDir','reverse',FontSize=10);
  xticklabels([])

%D1
nexttile
  plot(0:365,whD11,Linewidth=2)%,Linewidth=2
  hold on
  plot(0:365,whD12,Linewidth=2);
  plot(0:365,whD13,Linewidth=2);
  plot(0:365,whD14,Linewidth=2);
  yl=ylim;
  shadex=[365/12 5*365/12;4*365/12 7*365/12;4*365/12 7*365/12;365/12 5*365/12];
  shadey=[yl(1) yl(1);      yl(1) yl(1);      yl(2) yl(2);      yl(2) yl(2)];
  c=[211,211,211]/255;
  fill(shadex,shadey,c, 'linestyle','none','FaceAlpha','0.3')
  hold off
  xlim([0 365])
  xticks([0 365/12 2*365/12 3*365/12 4*365/12 5*365/12 6*365/12 7*365/12 8*365/12 9*365/12 10*365/12 11*365/12 365]);
  xticklabels([])
  ylabel({'Defence specialist';'D1';' mmol N m^{-3}'})
  set(gca,FontSize=10)
nexttile
plot(RSN,[spD11 spD12 spD13 spD14],'-o',LineWidth=2,MarkerSize=4,MarkerFaceColor='#0072BD')
xticks(flip(RSN))
set(gca,'XDir','reverse',FontSize=10);
ylabel({'Relative';'biomass'})
xticklabels([])
  nexttile
  plot(RSN,[suD11 suD12 suD13 suD14],'-o',LineWidth=2,MarkerSize=4,MarkerFaceColor='#0072BD')
  xticks(flip(RSN))
  set(gca,'XDir','reverse',FontSize=10);
  xticklabels([])

%D2
nexttile
  plot(0:365,whD21,Linewidth=2)%,Linewidth=2
  hold on
  plot(0:365,whD22,Linewidth=2);
  plot(0:365,whD23,Linewidth=2);
  plot(0:365,whD24,Linewidth=2);
  yl=ylim;
  shadex=[365/12 5*365/12;4*365/12 7*365/12;4*365/12 7*365/12;365/12 5*365/12];
  shadey=[yl(1) yl(1);      yl(1) yl(1);      yl(2) yl(2);      yl(2) yl(2)];
  c=[211,211,211]/255;
  fill(shadex,shadey,c, 'linestyle','none','FaceAlpha','0.3')
  hold off
  xlim([0 365])
  xticks([0 365/12 2*365/12 3*365/12 4*365/12 5*365/12 6*365/12 7*365/12 8*365/12 9*365/12 10*365/12 11*365/12 365]);
  xticklabels([])
  ylabel({'D2';'mmol N m^{-3}'})
  set(gca,FontSize=10)
nexttile
plot(RSN,[spD21 spD22 spD23 spD24],'-o',LineWidth=2,MarkerSize=4,MarkerFaceColor='#0072BD')
xticks(flip(RSN))
set(gca,'XDir','reverse',FontSize=10);
ylabel({'Relative';'biomass'})
xticklabels([])
  nexttile
  plot(RSN,[suD21 suD22 suD23 suD24],'-o',LineWidth=2,MarkerSize=4,MarkerFaceColor='#0072BD')
  xticks(flip(RSN))
  set(gca,'XDir','reverse',FontSize=10);
  xticklabels([])

%D3
nexttile
  plot(0:365,whD31,Linewidth=2)%,Linewidth=2
  hold on
  plot(0:365,whD32,Linewidth=2);
  plot(0:365,whD33,Linewidth=2);
  plot(0:365,whD34,Linewidth=2);
  yl=ylim;
  shadex=[365/12 5*365/12;4*365/12 7*365/12;4*365/12 7*365/12;365/12 5*365/12];
  shadey=[yl(1) yl(1);      yl(1) yl(1);      yl(2) yl(2);      yl(2) yl(2)];
  c=[211,211,211]/255;
  fill(shadex,shadey,c, 'linestyle','none','FaceAlpha','0.3')
  hold off
  xlim([0 365])
  xticks([0 365/12 2*365/12 3*365/12 4*365/12 5*365/12 6*365/12 7*365/12 8*365/12 9*365/12 10*365/12 11*365/12 365]);
  xticklabels([])
  ylabel({'D3';'mmol N m^{-3}'})
  set(gca,FontSize=10)
nexttile
plot(RSN,[spD31 spD32 spD33 spD34],'-o',LineWidth=2,MarkerSize=4,MarkerFaceColor='#0072BD')
xticks(flip(RSN))
set(gca,'XDir','reverse',FontSize=10);
ylabel({'Relative';'biomass'})
xticklabels([])
  nexttile
  plot(RSN,[suD31 suD32 suD33 suD34],'-o',LineWidth=2,MarkerSize=4,MarkerFaceColor='#0072BD')
  xticks(flip(RSN))
  xticklabels([])
  set(gca,'XDir','reverse',FontSize=10);
  xticklabels([])

%D4
nexttile
  plot(0:365,whD41,Linewidth=2)%,Linewidth=2
  hold on
  plot(0:365,whD42,Linewidth=2);
  plot(0:365,whD43,Linewidth=2);
  plot(0:365,whD44,Linewidth=2);
  yl=ylim;
  shadex=[365/12 5*365/12;4*365/12 7*365/12;4*365/12 7*365/12;365/12 5*365/12];
  shadey=[yl(1) yl(1);      yl(1) yl(1);      yl(2) yl(2);      yl(2) yl(2)];
  c=[211,211,211]/255;
  fill(shadex,shadey,c, 'linestyle','none','FaceAlpha','0.3')
  hold off
  xlim([0 365])
  xticks([0 365/12 2*365/12 3*365/12 4*365/12 5*365/12 6*365/12 7*365/12 8*365/12 9*365/12 10*365/12 11*365/12 365]);
  xticklabels([])
  ylabel({'D4';'mmol N m^{-3}'})
  set(gca,FontSize=10)
nexttile
  plot(RSN,[spD41 spD42 spD43 spD44],'-o',LineWidth=2,MarkerSize=4,MarkerFaceColor='#0072BD')
  xticks(flip(RSN))
  set(gca,'XDir','reverse',FontSize=10);
  ylabel({'Relative';'biomass'})
  xticklabels([])
nexttile
  plot(RSN,[suD41 suD42 suD43 suD44],'-o',LineWidth=2,MarkerSize=4,MarkerFaceColor='#0072BD')
  xticks(flip(RSN))
  set(gca,'XDir','reverse',FontSize=10);
  xticklabels([])

%D5
nexttile
  plot(0:365,whD51,Linewidth=2)%,Linewidth=2
  hold on
  plot(0:365,whD52,Linewidth=2);
  plot(0:365,whD53,Linewidth=2);
  plot(0:365,whD54,Linewidth=2);
  yl=ylim;
  shadex=[365/12 5*365/12;4*365/12 7*365/12;4*365/12 7*365/12;365/12 5*365/12];
  shadey=[yl(1) yl(1);      yl(1) yl(1);      yl(2) yl(2);      yl(2) yl(2)];
  c=[211,211,211]/255;
  fill(shadex,shadey,c, 'linestyle','none','FaceAlpha','0.3')
  hold off
  xlim([0 365])
  legend('Si:N=1','Si:N=0.7','Si:N=0.4','Si:N=0.1')
  legend boxoff
  xticks([0 365/12 2*365/12 3*365/12 4*365/12 5*365/12 6*365/12 7*365/12 8*365/12 9*365/12 10*365/12 11*365/12 365]);
  xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
  xlabel('Month')
  ylabel({'Competitor specialist';'D5';' mmol N m^{-3}'})
  set(gca,FontSize=10)
nexttile
plot(RSN,[spD51 spD52 spD53 spD54],'-o',LineWidth=2,MarkerSize=4,MarkerFaceColor='#0072BD')
xticks(flip(RSN))
set(gca,'XDir','reverse',FontSize=10);
ylabel({'Relative';'biomass'})
xlabel('Si:N ratio')
  nexttile
  plot(RSN,[suD51 suD52 suD53 suD54],'-o',LineWidth=2,MarkerSize=4,MarkerFaceColor='#0072BD')
  xticks(flip(RSN))
  xlabel('Si:N ratio')
  set(gca,'XDir','reverse',FontSize=10);

end

end