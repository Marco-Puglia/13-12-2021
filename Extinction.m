function [v0]=Extinction(PSindex,ext);
load('steadystate11-T5EiEtau1.666gamma0.5.mat')

%Default (No extinction)
extindex = [];

%Choose who gets extinct

if ext == 1
    extindex = reshape(PSindex(1:length(PSindex),1),1,[]);                           %AUTOTROPHS
end
if ext == 2
    extindex = reshape(PSindex(1:length(PSindex),2:length(PSindex)-1),1,[]);         %MIXOTROPHS
end
if ext == 3
    extindex = reshape(PSindex(1:length(PSindex),length(PSindex)),1,[]);             %HETEROPHS
end
if ext == 4
    extindex = reshape(PSindex(1:length(PSindex),1:length(PSindex)-1),1,[]);         %AUT-MIX
end
if ext == 5
    extindex = reshape(PSindex(1:length(PSindex),2:length(PSindex)),1,[]);           %HETER-MIX
end
if ext == 6
    extindex = reshape([PSindex(1:length(PSindex),1) PSindex(1:length(PSindex),length(PSindex))],1,[]);   %AUT-HETER
end

if ext == 7
    extindex = [];                                                      %NONE
end

if ext == 8
    glue = cat(2,[PSindex(1:length(PSindex),1:5) PSindex(1:length(PSindex),9:12)]);              %AUT-MIX-HETER
    extindex = reshape([glue PSindex(1:length(PSindex),17:length(PSindex))],1,[]);               %AUT-MIX-HETER
end


v0(1)=sum(v0(extindex));
v0(extindex+1) = 0;

clear glue

end