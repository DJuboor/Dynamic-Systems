%% Laplace Workshop

%% Part A, Definition of the Laplace Transform.

% Given a function f(t) in the time domain, its one-sided Laplace Transform 
% defined by the following integral:
% 
% L{f(t)} = int(exp(-s*t) * f, 0,inf);
% 
% The function exp(-s*t) is called the kernel or nucleus of the transform,
% There are many other useful integral transforms including the Fourier
% transform, Fourier sine and cosine transforms, Hartley transform, Mellin
% transform, Weierstrass transform, Hankel transform, Abel transform amd the
% Hilbert transform; all defined using different kernels.
% 
% 
% In this section, we find a few Laplace transforms by directly evaluating
% the above integral.
% 
% EXAMPLE: Show that the Laplace transform of the constant function 
% f(t) = 1    is    L{1} = 1 / s
% Directly evaluate the above integral. Tip: Use inf to denote infinity in
% the limits of integration.

        syms s t; assume(real(s)>0)
        f = 1;
        L = int( exp(-s*t) * f, 0, inf);

% After evaluating this code, you should see that L = 1 / s
% Completely agreeing with a Table of Laplace Transforms

% Lets make this a function for later use.
        Temp_Laplace = @(f)[int( exp(-s*t) * f, 0, inf)];
% This will work as a generic function as long as we remember to define
% t and s as a symbolic and put the assumptions on s

% Lets try it out on a few more examples:
        syms s t; assume(real(s)>0)
        A = Temp_Laplace(t);
        B = Temp_Laplace(t^2);
        C = Temp_Laplace(exp(2*t));

%% Assumptions
% You may have noticed that the answer to C looks odd and does not match your
% Table of Laplace Transforms. That is because we have not told MATLAB enough
% information about the variable s. Above, we only said it was a "symbol",
% and that reals(s) was positive. 
% 
% Isaac Asimov once said: "Your assumptions are your windows on the world.
% Scrub them off every once in a while, or the light won't come in." 
% 
% Lets try "scrubbing them off" and add the assumption (after declaring 
% s and t as symbols):
%         >> assume(s>2)
% And rerun to see what we get.
        syms s t; assume(s>2)
        C = Temp_Laplace(exp(2*t));

%% More Assumptions
% Lets clear our assumptions and try some other examples where
% our answers will require suitable assumptions to reproduce the form
% given in a Table of Laplace Transforms.
        assume(s,'clear')

        syms s t n; assume(n>-t); assume(s>0)
        A = Temp_Laplace(t^n);
        % In our case-- for positive integer aruments-- gamma(n+1) = n!
        % Therefore our answer matches the Table of Laplace Transfoms!  
        
        assume(s,'clear')
        syms s t a; assume(t<1); assume(a<1); assume(s>0)
        B = Temp_Laplace(cos(a*t));
        % In this case, two extra assumptions were needed for the restricitons on
        % the Cosine function.
        
        assume(s,'clear')
        syms s t a; assume(a>0); assume(t>0); assume(s > a)
        C = Temp_Laplace(exp(a*t));
        % In this case, two extra assumptions were needed for the restricitons on
        % the Exponential function.

%% Laplace Transform in MATLAB
% The good news, is that MATLAB has a built-in command named laplace(), so
% you won't need to find these transforms using the defining integral, which
% is best used for demonstrating the fundamental properties of the transform.
% But be sure, you can always use a Laplace Transform Table at any time.
%
% Lets try some general functions with laplace():
clear, clc

        syms t;
        A = laplace(sinh(2*t));
        B = laplace(sqrt(t));
        
        %How about a unit-step fuction? We'll use heaviside(t) to denote
        %this.
        C = laplace(heaviside(t-10));
 
% Be sure to check your solutions with the Transform Table.
% Easy, huh?

%% Inverse Transform
% Just as important as the Laplace transform, is its inverse transform. In
% MATLAB, this is found using the ilaplace() command.
% 
% Lets find the inverse Laplace transform for some example functions that are
% already defined in the s-domain.
        
        syms s;
        A = ilaplace((s+4)/(s^2+4));
        B = ilaplace(1/sqrt(s));
        
        % Because ilaplace is undefined for input arguments of the type
        % "double", we'll have to use (1+0*s) 
        C = ilaplace(1+0*s);

% Hint: You can print an "ugly" mathematical answer (f) in a nicer format
% using:    
% >> pretty(f)

%% Partial Fraction Expansions
% Partial fraction expansions are absolutely necessary so you can find
% inverse Laplace transforms using the standard Laplace tables. Fortunately,
% MATLAB now has the build-in command partfrac(). Let's see how that works
% now.
%
% %WARNING: partfrac is a fairly recent command. Users with an older version
% of MATLAB may beed to use the following command instead:
%    >>    feval(symengine, 'partfrac', F)
%   
% Let's do some examples:
        syms s;
        A = partfrac(1/(s^2+s));
        B = partfrac((10*(s+3))/((s-2)*(s+1)*(s-5)));
        C = partfrac((2*s+4)/(s*(s^2+1)));
        
%% Solving a Differential Equation using the Laplace Transform
%
% Given a DE: y'' + y = sin(2*t) and the ICs of y(0)=2, y'(0)=1
% We can use the method of Laplace Transforms to find the solution to the
% linear differential equation.
%
% First, let's enter the differential equasion as usual, and find the exact
% solution using dsolve(). We will use this later to compare to the answer
% found using the Laplace method to confirm the same answer is given.

        syms y(t)
        Dy = diff(y,t) ; D2y = diff(y,t,t);
        DE = D2y + y == sin(2*t);

        Soln = dsolve(DE, y(0)==2, Dy(0) == 1);
        % Let's simplify this:
        Smpl_ds_Soln = simplify(Soln); 

        
% % If you're the type who likes to visualize, uncomment this block and
% % run!
%   
%     % Vectorized Exact solutions found from above:
%     Y = @(t)cos(t) + (8.*sin(t))./3 - (2.*cos(t).*sin(t))./3;
%     Dy = @(t)(8.*cos(t))./3 - (2.*cos(2.*t))./3 - sin(t);
%     
%     % Let's plot our solution!
%     figure
%     subplot(2,2,[1,3])
%     t = 0: 0.01 : 17;
%     
%     hold on
%     plot(t, Y(t), 'b', 'LineWidth',2)
%     plot(t, Dy(t), 'r', 'LineWidth',2)
%     plot(0, 2, 'blueo', 'MarkerFaceColor', 'Yellow') %IC
%     plot(0, 1, 'redo', 'MarkerFaceColor', 'Yellow') %IC
%     axis([0 16 -4 4])
%     title('Exact Solution and Derivative')
%     legend('Y(t)', 'Dy(t)', 'Location', 'Best')
%     grid on
%     box on
%     
%     subplot(1,2,2)
%     hold on
%     plot(Y(t), Dy(t), 'k', 'LineWidth', 2)
%     title('Phase Plot')
%     xlabel('Y(t)')
%     ylabel('Dy/Dt')
%     plot(1, 2, 'blacko', 'MarkerFaceColor', 'Cyan') %IC
%     grid on
%     box on
% %  

%% Laplace Transform Method for Solution
%
% Given the exact same DE and initial conditions, let's use our newfound
% skills in the Laplace Method to do the exact same thing!


        syms s t Y 
        
        % In this example, we will use Y(s) to denote the transform of the
        % unknown function y(t).

        % Next, we'll find the Laplace transform of y'(t): Y1 = s Y - y(0)
        % This is necessary, even though this term does not appear in the
        % LHS of the differential equation.
        
        % y(0)=2
        Y1 = s*Y - 2;
        
        % y'(0)=1
        Y2 = s*Y1 - 1; 
        
        % Find the Laplace transform F of the forcing term f(t) = sin(2*t)
        F = laplace( sin(2*t) );
        
        % Next, we'll combine all of the terms into the transform of the entire equation,
        % which we will name LTofDE for Laplace Transform of DE.
        LTofDE = Y2 + Y == F;

        % Lastly, we just solve as Normal
        Soln = solve(LTofDE, Y); 
        LT_Soln = ilaplace(Soln); LT_Soln = simplify(A);
        

%% Let's Compare:
disp 'The simplified general solution using the dsolve method'
pretty(Smpl_ds_Soln)

disp 'The simplified feneral solution using the Laplace Transform method'
pretty(LT_Soln)

