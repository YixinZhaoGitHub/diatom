%% light function
%
function I = getLight(P,D1,D2,D3,D4,D5,Dp,Dd,t,p)

         int=(cumsum((P+D1+D2+D3+D4+D5+Dp+Dd).*p.dz))*p.k;

         I=p.Iin.*(0.4*sin(2*pi*t/365-pi/2)+0.6).*exp(-int-p.Kbg.*p.z');
end