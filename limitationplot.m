%Figure 4
figure(13)%P
 nexttile
 plot(sim4.p.An.*sim4.N(length(sim4.t)-185,:)./(sim4.p.An.*sim4.N(length(sim4.t)-185,:)+sim4.p.pmax),-sim4.p.z,'--b',LineWidth=2)
  hold on
  plot(sim4.p.Al.*getLight(sim4.P(length(sim4.t)-185,:)',sim4.D1(length(sim4.t)-185,:)',sim4.D2(length(sim4.t)-185,:)',sim4.D3(length(sim4.t)-185,:)',sim4.D4(length(sim4.t)-185,:)',sim4.D5(length(sim4.t)-185,:)',sim4.Dp(length(sim4.t)-185,:)',sim4.Dd(length(sim4.t)-185,:)',sim4.t(length(sim4.t)-185),sim4.p)./(sim4.p.Al.*getLight(sim4.P(length(sim4.t)-185,:)',sim4.D1(length(sim4.t)-185,:)',sim4.D2(length(sim4.t)-185,:)',sim4.D3(length(sim4.t)-185,:)',sim4.D4(length(sim4.t)-185,:)',sim4.D5(length(sim4.t)-185,:)',sim4.Dp(length(sim4.t)-185,:)',sim4.Dd(length(sim4.t)-185,:)',sim4.t(length(sim4.t)-185),sim4.p)+sim4.p.pmax),-sim4.p.z,'--r',LineWidth=2)
  hold off
%   legend('Nitrogen','Light')
  ylabel('Depth (m)')
  title('Non-diatom phytoplankton')
  set(gca,FontSize=20)

nexttile%D5
  plot(sim4.p.And(5).*sim4.N(length(sim4.t)-185,:)./(sim4.p.And(5).*sim4.N(length(sim4.t)-185,:)+sim4.p.pmaxd(5)),-sim4.p.z,'--b',LineWidth=2)
  hold on
  plot(sim4.p.Asd(5).*sim4.S(length(sim4.t)-185,:)./(sim4.p.Asd(5).*sim4.S(length(sim4.t)-185,:)+sim4.p.pmaxd(5)),-sim4.p.z,'--k',LineWidth=2)
  plot(sim4.p.Ald(5).*getLight(sim4.P(length(sim4.t)-185,:)',sim4.D1(length(sim4.t)-185,:)',sim4.D2(length(sim4.t)-185,:)',sim4.D3(length(sim4.t)-185,:)',sim4.D4(length(sim4.t)-185,:)',sim4.D5(length(sim4.t)-185,:)',sim4.Dp(length(sim4.t)-185,:)',sim4.Dd(length(sim4.t)-185,:)',sim4.t(length(sim4.t)-185),sim4.p)./(sim4.p.Ald(5).*getLight(sim4.P(length(sim4.t)-185,:)',sim4.D1(length(sim4.t)-185,:)',sim4.D2(length(sim4.t)-185,:)',sim4.D3(length(sim4.t)-185,:)',sim4.D4(length(sim4.t)-185,:)',sim4.D5(length(sim4.t)-185,:)',sim4.Dp(length(sim4.t)-185,:)',sim4.Dd(length(sim4.t)-185,:)',sim4.t(length(sim4.t)-185),sim4.p)+sim4.p.pmaxd(5)),-sim4.p.z,'--r',LineWidth=2)
  hold off
  legend('Nitrogen','Silicate','Light',Location="southeast")
%   ylabel('Depth (m)')
  title('Competitor specialist (D5)')
 set(gca,FontSize=20)