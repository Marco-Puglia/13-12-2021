A = 1:441;
A = reshape(A,[21,21]);
A = A';

%CENTER              iseed=[12,33,54,75,96,117,138,159];
%ALL                 iseed=[1:441];
%SIDES-HORIZONTAL    iseed=[1:147 315:441];
%SIDES-VERTICAL      iseed=reshape([A(1:21,1:7);A(1:21,15:21)],1,[]);
%MIDDLE-VERTICAL     iseed=reshape(A(1:21,8:14),1,[]);

deathindex=reshape(A(1:21,15:21),1,[]);

    if k == 10
       v0(1)=sum(v0(deathindex));
       v0(deathindex)=0;
    end    
    if k == 50
        deathindex=reshape(A(11,1:21),1,[]);
        deathindex=deathindex+1;
       v0(1)=sum(v0(deathindex));
       v0(deathindex)=0;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%to use delete v0=[N_0;P_0]; in MCM.m

load("steadystate.mat")


A = 1:441;
A = reshape(A,[21,21]);
A = A';
deathindex=reshape(A(1:21,8:15),1,[]);
deathindex=deathindex+1;

    if k == 10
       v0(1)=sum(v0(deathindex));
       v0(deathindex)=0;
    end