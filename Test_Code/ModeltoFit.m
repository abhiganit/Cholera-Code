clear;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % General model specifications
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Specify XU the model being fit
                    % XU(1)- beta_0
                    % XU(2) - population density
                    % XU(3) - number of health facilities 
                    % XU(4) - Past incidence
                    % XU(5) - Product of incidence and attacks
                    % XU(6) - Product of incidence and conflict
                    % XU(7) - Product of incidence and rainfall
                    % XU(8) - Rainfall only        
                    % XU(9) - Incidence in other govnorates
                    % XU(10) - Attacks only
                    % XU(11) - Rebel control
            XU=zeros(1,11);
            XU([1 2 4 7 8 9])=1;
            % Specify the lag
            % tau(1) - population density incidence
            % tau(2) - health zone incidence
            % tau(3) - Past incidence
            % tau(4) - Product of incidence and attacks
            % tau(5) - Product of incidence and conflict
            % tau(6) - Product of incidence and rainfall
            % tau(7) - Perciptiation only
            % tau(8) - Incidence in other govneroates
            % tau(9)- Attack only
            % tau(10)- Rebel control

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % The formation of the environmental function
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %Specify the attack function to be used
             % AF=0 attack only has effect before; AF=1 Attack has effect only after; AF=2; Attack has effect before and after

            %Specify the conflict function to be used
             % CF=0 linear effect; CF=1 Hill function with n=1; CF=2; Full hill function; 

            % Specify the rainfall function to be used in rainfall*incidence
             % RF=0 Increased incidence for low-rainfall; RF=1 increased incidence for high rainfall; RF=2 increased incidence for high and low rain fall


            % Specify the rainfall function to be used
            % RF=0 Increased incidence for low-rainfall; RF=1 increased incidence for high rainfall; RF=2 increased incidence for high and low rain fall

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            % Specify plots generated
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Specify the govenorates want to examine (1-22)
            %1- 'Abyan' 
            %2- Aden' 
            %3- Al Bayda' 
            %4- Al Dhale''e' 
            %5- Al Hudaydah' 
            %6- Al Jawf' 
            %7- Al Maharah' 
            %8- Al Mahwit' 
            %9- Amanat Al Asimah' 
            %10- Amran' 
            %11- Dhamar' 
            %12- Hadramaut' 
            %13- Hajjah' 
            %14- Ibb' 
            %15- Lahj' 
            %16- Marib' 
            %17- Raymah' 
            %18- Sa''ada' 
            %19- Sana''a' 
            %20- Shabwah' 
            %21- Socotra' 
            %22- Taizz'
            G=[]; 

            % Specify if want plot of fit
            PF=0; %PF=0 does not produce plot

            % Specify if want plot of function
            PE=1; %PE=0 does not produce plot

            % Specify if want plot of projection and fit
            PP=1; %PP=0 does not produce plot

            % Specigy if table and fitting information displayed 
            DT=1; % DT=0; does not display this information

            
            ProFittingGA(XU,G,PE,PP,DT)