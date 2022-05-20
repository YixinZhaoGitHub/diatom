% Figure 2 in the report.
%% light+diffusivity
tiledlayout(2,1)
%light
nexttile
t=0:365;
p.Iin=200;
plot(t,p.Iin.*(0.4*sin(2*pi*t/365-pi/2)+0.6),'k',LineWidth=3)
   xlim tight
    ylim ([0 225])
   xticks([0 365/12 2*365/12 3*365/12 4*365/12 5*365/12 6*365/12 7*365/12 8*365/12 9*365/12 10*365/12 11*365/12 365]);
   xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
%    xlabel('Month')
   ylabel('Light intensity (W m^{-2})')
   legend('Light','Location','northeast');
   legend('boxoff')
   title('Incident light intensity')
set(gca,FontSize=20)

%diffusivity
p.diffml=100;%mixed layer diff     diffusivity m2/day 
p.diffdl=4;%deep layer diff
p.mld=50;% for constant simulation
p.w=1;% themocline width
p.Zm=300; % depth
p.dz = 2; % grid cell size
p.z = (0+p.dz/2):p.dz:(p.Zm-p.dz/2); %grid vector
fash=load('mld_fasham_mat.mat');
zmld1=fash.alk3;
zmld1=zmld1(1:end-1);
zmld=[zmld1; repmat(zmld1(1:end), 365/365-1,1)]; %mixed layer depth
zmld(end+1)=zmld(1);
zmldgrid=round(zmld/p.dz);% mixed layer corresponding grid number

for i=1:365
diff=p.diffdl+((p.diffml-p.diffdl)./(1+(exp(p.z-zmld(i))./p.w)));
diffmat(:,i)=diff;
end

diffmat(:,end+1)=p.diffdl+((p.diffml-p.diffdl)./(1+(exp(p.z-zmld(end))./p.w)));

% nexttile
% plot(0:365, -zmld)
% ylim([-300 0])
% xlim([0 365])

nexttile
surface(0:365,-p.z,diffmat(:,1:366))
  xlim([0 365])
shading interp
  bar=colorbar;
  bar.Label.String = 'Diffusivity (m^2 d^{-1})';
  xticks([0 365/12 2*365/12 3*365/12 4*365/12 5*365/12 6*365/12 7*365/12 8*365/12 9*365/12 10*365/12 11*365/12 365]);
  xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
  xlabel('Month')
  ylabel('Depth (m)')
  title('Diffusivity in the water column')
  set(gca,FontSize=20)
