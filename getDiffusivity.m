function [diffmat,zmld,zmldgrid] = getDiffusivity(t,p)

%load file for the mixed layer depth   Bermuda 1960     Camila NUM
fash=load('mld_fasham_mat.mat');

zmld1=fash.alk3;
zmld1=zmld1(1:end-1);
zmld=[zmld1; repmat(zmld1(1:end), t/365-1,1)]; %mixed layer depth
zmld(end+1)=zmld(1);
zmldgrid=round(zmld/p.dz);% mixed layer corresponding grid number

for i=1:t
diff=p.diffdl+((p.diffml-p.diffdl)./(1+(exp(p.z-zmld(i))./p.w)));
diffmat(:,i)=diff;
end

diffmat(:,end+1)=p.diffdl+((p.diffml-p.diffdl)./(1+(exp(p.z-zmld(end))./p.w)));
diffmat(1,:)=[];

end


