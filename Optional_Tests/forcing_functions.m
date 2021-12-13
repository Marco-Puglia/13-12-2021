function [ forcing_fcns ] = forcing_functions ( )
    forcing_fcns.stepfunction   = @stepfunction;
    forcing_fcns.sinusoidal   = @sinusoidal;
    forcing_fcns.squarewave   = @squarewave;
end

%%
function T = stepfunction(t)
    T(t>10) =25;
    T(t<=10)=14;
end

%%
function T = sinusoidal(t)
    T = 20 + 15 * sind(t/5);
end

%%
function T = squarewave(t)
    T = 15 + 5 * square(t.*2.*pi);
    
    T(t<10) =15;
end




