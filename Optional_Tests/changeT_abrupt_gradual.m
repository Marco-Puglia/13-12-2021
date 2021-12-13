    %%%%%%%%%%%%%%%%%%%%%%%%%%%%CHOSE BETWEEN GRADUAL AND ABRUPT CHANGE OF TEMPERATURE%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%abrupt 20-40
    if k<10
        eco_pars.T=20;
    else eco_pars.T=40;
    end
    
    %%%gradual 20-30 in 100 years
    CC=20:0.1:30;
    if k<100
     = CC(k);
    else eco_pars.T = 30;
    end        
    
    %%%gradual 30-20 in 100 years
    CC=20:0.1:30;
    complete=cat(2,CC,flip(CC));
    if k<201
         eco_pars.T = complete(k);
    else eco_pars.T = 20;
    end
    
    %%%abrupt 20-50
    if k<10
         eco_pars.T=20;
    else eco_pars.T=50;
    end
    
    %%%abrupt 20-30
    if k<10
         eco_pars.T=20;
    else eco_pars.T=30;
    end
    
    %%%abrupt 20-10
    if k<10
         eco_pars.T=20;
    else eco_pars.T=10;
    end
    
    %%%abrupt 20-15
    if k<10
         eco_pars.T=20;
    else eco_pars.T=15;
    end
    
    %%%abrupt 20-18
    if k<10
         eco_pars.T=20;
    else eco_pars.T=18;
    end
    
    %%%abrupt 20-21
    if k<10
         eco_pars.T=20;
    else eco_pars.T=21;
    end   
    
    %%%abrupt 20-23
    if k<10
         eco_pars.T=20;
    else eco_pars.T=23;
    end      
    
    %%%abrupt 20-22
    if k<10
        eco_pars.T=20;
    else eco_pars.T=22;
    end        
    
    %%%abrupt 20-20.5
    if k<10
         eco_pars.T=20;
    else eco_pars.T=20.5;
    end           
    
    %%%abrupt 20-20.8
    if k<10
         eco_pars.T=20;
    else eco_pars.T=20.8;
    end           
    
    %%%abrupt 20-20.9
    if k<10
         eco_pars.T=20;
    else eco_pars.T=20.9;
    end              
    
    %%%abrupt 20-30
    if k<10
         eco_pars.T=20;
    else eco_pars.T=30;
    end   
    
    %ext: Kill 1=AUTOTROPHS, 2=MIXOTROPHS, 3=HETEROPHS, 4=AUT-MIX, 5=HETER-MIX, 6=AUT-HETER, 7=NONE, 8=AUT-MIX-HETER
    %gradual 20-30 in 100 years
    CC=20:0.1:30;
    if k<100
        T = CC(k);
    else T = 30;
    end        
    if k==100
        ext = 1;
        [v0] = Extinction(PSindex,ext);
    end
            
    if k==10
        ext = 6;
        [v0] = Extinction(PSindex,ext);
    end