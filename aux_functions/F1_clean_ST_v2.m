% F1_clean_ST: importa tall i ordena les series temporals de estacions
% meteorlògiques o mareografs
%--------------------------------------------------------------------------
% 1- Importar dades
% 2- Tallar per el tempts dessitjat
% 3- Ordenar
% 4- Eliminar estacions sense dades
% 5- Omplir buits
function [dades,tmax]=F1_clean_ST_v2(path_mat,arx_dades,variables,tini,tfi,t_change)
% variables ha de tenir la següent forma:
% {'estacions','lon','lat','time','Altres ST'}

dades=load([path_mat,arx_dades]);

% eliminam variables innecessaries
aux=fieldnames(dades);
aux=setdiff(aux,variables);
for nv=1:length(aux)
    dades=rmfield(dades,aux{nv});
end

%--- Reetallem el temps en l'interval dessitjat
aux=(dades.time>=tini & dades.time<=tfi);

%--- Retallem les variables que necessitem
for var=variables(4:end)
    var=char(var);
    eval(['dades.',var,'=dades.',var,'(aux,:);']);
end

%--- Eliminam les estacions que no tenen dades
eval(['aux=isnan(dades.',variables{5},');']);
aux=sum(aux,1)/length(dades.time);
aux=(aux>0.05);
for nv=1:3
    eval(['dades.',variables{nv},'(aux)=[];']);
end
for nv=5:length(variables)
    eval(['dades.',variables{nv},'(:,aux)=[];']);
end

%--- Eliminam les dades sense variança
eval(['aux=diff(dades.',variables{5},',1,1)==0;']);
aux=sum(aux,1)/length(dades.time);
aux=(aux>0.95);
for nv=1:3
    eval(['dades.',variables{nv},'(aux)=[];']);
end
for nv=5:length(variables)
    eval(['dades.',variables{nv},'(:,aux)=[];']);
end

eval(['Ne=length(dades.',variables{1},');']);

%--- Omplim els nans
for nv=5:length(variables)
    eval(['aux=dades.',variables{nv},';']);
for ne=1:Ne
    aux(:,ne)=rellenar_huecos(dades.time,aux(:,ne),120/24/60,'autorregresivo');
end
    eval(['dades.',variables{nv},'=aux;']);
end

%--- Eliminam les estacions que no tenen dades
eval(['aux=isnan(dades.',variables{5},');']);
aux=sum(aux,1);
aux=(aux>=1);
for nv=1:3
    eval(['dades.',variables{nv},'(aux)=[];']);
end
for nv=5:length(variables)
    eval(['dades.',variables{nv},'(:,aux)=[];']);
end


%--- Ordenem les estacions de mes propera a més llunyana de la estació de
%referència
[a,aux]=sort(dades.lat); 
for nv=1:3
    eval(['dades.',variables{nv},'=dades.',variables{nv},'(aux);']);
end
for nv=5:length(variables)
    eval(['dades.',variables{nv},'=dades.',variables{nv},'(:,aux);']);
end

if t_change
    aux=find(contains(dades.mareografs,'ciutadella'),1,'first');
%--- Centram la serie temporal al temps de la rissaga
[y,Amax,tmax]=F2_filt(dades.time,dades.SL(:,aux),'finestra',8*60,'f_filt_env',24*60);
tini=tmax-1; tfi=tmax+1;
tini=tmax-2; tfi=tmax+4;

%--- Reetallem el temps en l'interval dessitjat
aux=(dades.time>tini & dades.time<tfi);

%--- Retallem les variables que necessitem
for var=variables(4:end)
    var=char(var);
    eval(['dades.',var,'=dades.',var,'(aux,:);']);
end
else
    tmax=[];
end

end