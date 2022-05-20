function A=getAffinity %Hansen and Visser (2019)
v=linspace(0.4,0.8,5);%vacuole fraction
r= linspace(25,5,5); % cell radius um
rb=r.*(1-v).^(1/3);%spherical equivalent radius of the cytoplasm um
cb=0.125*10^-6;%Carbon density of biomass μg C μm−3
cbN=0.125*10^-12 /12*1000*0.13/1.05;%Carbon density of biomass pg C μm−3 to mmol N to μm−3     Si:N=1.05 Si:C=0.13--- N:C=0.13/1.05
cbS=0.125*10^-12/12*1000*0.13;%Carbon density of biomass pg μm−3 mmol Si to μm−3  
dif=1.62 * 10^-4*10^12; % um2 d−1
cl=0.01*cb^(2/3); %µg C day-1 (W/m2)-1 µm3^(-2/3)) quantum yield of photosynthesis
A.nitrogen=(3*dif./(rb.^2*cbN)).*(1-v).^(-1/3).*1E-18;%m3 mmol N-1 d-1
A.silicon=(3*dif./(rb.^2*cbS)).*(1-v).^(-1/3).*1E-18; %m3 mmol Si-1 d-1
A.light=3*cl./(4*rb*cb).*(1-v).^(-2/3); % (W m-2)-1 day-1
end

