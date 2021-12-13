clear

%%
load('matlab.mat')

a=reshape(twotwenty,[2,20]);
b=reshape(twentywenty,[2,20]);

c=cat(2,a(1,:),a(2,:));
d=cat(2,b(1,:),b(2,:));

plot(c)
hold on
plot(d)
legend('Grid 2x20','Grid 20x20')
text(1,-0.01,'NANOPHYTOPLANKTON            MICROPHYTOPLANKTON           MESOPHYTOPLANKTON')
text(21,-0.01,'NANOZOOPLANKTON            MICROZOOPLANKTON            MESOZOOPLANKTON')
text(20,-0.01,'|')
xlabel('Compartments') 
ylabel('Biomass quantity') 




