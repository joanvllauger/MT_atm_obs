% F3_correlation
% Funció de correlació entre totes les estacions que tenim
%--------------------------------------------------------------------------
function [est_corr_max,est_corr_mlag]=F3_correlation_v6(y,D_matrix,R)
[m,Ne]=size(y);

est_corr_max=zeros(Ne,Ne);
est_corr_mlag=zeros(Ne,Ne);
for ne=1:Ne
    for nei=ne:Ne
        d=D_matrix(ne,nei);
        if d<=R
            maxlag=ceil(d/10*1000/60);
        [lag,corr]=correlacion_lag(y(:,ne),y(:,nei),maxlag);
        [m,I]=max(corr);
        est_corr_max(ne,nei)=m;
        est_corr_max(nei,ne)=m;
        est_corr_mlag(ne,nei)=lag(I);
        est_corr_mlag(nei,ne)=-lag(I);
        else
        est_corr_max(ne,nei)=NaN;
        est_corr_mlag(ne,nei)=NaN;
        est_corr_max(nei,ne)=NaN;
        est_corr_mlag(nei,ne)=NaN;
        end
    end
end