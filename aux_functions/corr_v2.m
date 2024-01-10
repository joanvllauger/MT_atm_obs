function [cc]=corr_v2(x,y)
aux=not(isnan(x)|isnan(y));
cc=mean(x(aux).*y(aux))/sqrt(mean(x(aux).^2)*mean(y(aux).^2));
end
