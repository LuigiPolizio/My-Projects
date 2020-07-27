%--------------------------------------------------------------------
% Dinamica e Simulazione di Volo - Coiro, De Marco A.A. 2010-2011
%--------------------------------------------------------------------
% ver. 1.1, 03-May-2011
%
% Updates: 
%          ver. 1.0, initial implementation
%          ver. 1.1, added e_T
%
classdef DSVAircraft % < handle
    % DSVAircraft class for aircraft data management
    %   Detailed explanation goes here
    
    properties
        %------------------------------------------------------------------
        % IDs and misc.
        %------------------------------------------------------------------
        Name = 'DSVAircraft - <Put a name here>';
        g = 9.81;
        err = 0;
        %------------------------------------------------------------------
        % Geometry
        %------------------------------------------------------------------
        S;
        b;
        mac; 
        %------------------------------------------------------------------
        % Mass, inertia, etc
        %------------------------------------------------------------------
        mass;
        W;
        k_y;
        mu_x; 
        Xcg_adim;
        Xn_adim;
        %------------------------------------------------------------------
        % Aerodynamics
        %------------------------------------------------------------------
        CD_0;
        K;
        m;
        CL_alpha;
        CL_delta_e;
        CL_delta_s;
        CL_alpha_dot;
        CL_q;
        Cm_0;
        Cm_alpha;
        Cm_alpha_dot;
        Cm_delta_s;
        Cm_delta_e;
        Cm_delta_e_dot;
        Cm_q;
        %------------------------------------------------------------------
        % Elevator
        %------------------------------------------------------------------
        S_e;
        x_C_e;
        Lambda_e;
        mac_e;
        mass_e;
        W_e;
        ec_adim;
        k_e;
        I_e;
        I_ey;
        Ch_e_0;
        Ch_e_alpha;
        Ch_e_delta_s;
        Ch_e_delta_e;
        Ch_e_delta_e_dot;
        Ch_e_q;
        Ch_e_alpha_dot;
        eps_0;
        deps_dalpha;
        Rs_e;
        Rg_e;
        delta_e_max;
        delta_e_min;
        %------------------------------------------------------------------
        % Propulsion
        %------------------------------------------------------------------
        T;
        Cm_T_0;
        Cm_T_alpha;
        mu_T;
        e_T;
        %------------------------------------------------------------------
        % Limitations
        %------------------------------------------------------------------
        CL_max;
        CL_min;
        n_max;
        n_min;
        Fe_max;
        Fe_min;
    end
    
    methods
        %% Constructor, populate class properties reading from file
        function obj = DSVAircraft(dataFileName)
            f_id = fopen(dataFileName,'r');
            % check if file opening was ok
            if (f_id==-1)
                obj.err = -1; % opening failed
                disp(['DSVAircraft :: initFromFile __ Could NOT open file ', dataFileName, ' ...'])
            else
                disp(['DSVAircraft :: initFromFile __ Opening file ', dataFileName, ' ... OK.'])
                %% File open, OK
                for i=1:5 % read 5 dummy rows
                    temp = fgetl(f_id);
                end
                %% Geometric data
                obj.S = fscanf(f_id,'%f ');
                temp = fgetl(f_id);% read the rest of the line as dummy text
                obj.b = fscanf(f_id,'%f '); temp=fgetl(f_id);
                obj.mac = fscanf(f_id,'%f '); temp=fgetl(f_id);
                for i=1:2
                    temp = fgetl(f_id);
                end
                %% Mass data
                obj.mass = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.W = obj.mass*obj.g;
                obj.k_y = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.mu_x = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.Xcg_adim = fscanf(f_id,'%f '); temp = fgetl(f_id);
                for i=1:3
                    temp = fgetl(f_id);
                end
                %% Aerodynamics
                % Neutral point of the aircraft, non-dimensional
                obj.Xn_adim = fscanf(f_id,'%f '); temp = fgetl(f_id);
                for i=1:1
                    temp = fgetl(f_id);
                end
                % Polare
                obj.CD_0 = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.K = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.m = fscanf(f_id,'%f '); temp = fgetl(f_id);
                for i=1:1
                    temp = fgetl(f_id);
                end
                % Aerodynamic gradients
                obj.CL_alpha = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.CL_delta_e = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.CL_delta_s = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.CL_alpha_dot = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.CL_q = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.Cm_0 = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.Cm_delta_e = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.Cm_delta_s = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.Cm_alpha_dot = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.Cm_q = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.Cm_delta_e_dot = fscanf(f_id,'%f '); temp = fgetl(f_id);
                % calculate stability derivative
                %obj.Cm_alpha = -obj.CL_alpha*(obj.Xn_adim - obj.Xcg_adim);
                obj.Cm_alpha = fscanf(f_id,'%f '); temp = fgetl(f_id);
                for i=1:3
                    temp = fgetl(f_id);
                end
                %% Elevator data
                % Geometry
                obj.S_e = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.Lambda_e = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.x_C_e = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.mac_e = fscanf(f_id,'%f '); temp = fgetl(f_id);
                for i=1:1
                    temp = fgetl(f_id);
                end
                % Mass, inertia, etc
                obj.mass_e = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.W_e = obj.mass_e*obj.g;
                obj.ec_adim = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.k_e = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.I_e = (obj.mass_e)*(obj.k_e^2);
                obj.I_ey = obj.mass_e*obj.ec_adim*obj.mac_e*obj.x_C_e - obj.I_e*sin(obj.Lambda_e);
                for i=1:1
                    temp = fgetl(f_id);
                end
                % Aerodynamics
                obj.Ch_e_0 = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.Ch_e_alpha = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.Ch_e_delta_s = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.Ch_e_delta_e = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.Ch_e_delta_e_dot = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.Ch_e_q = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.Ch_e_alpha_dot = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.eps_0 = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.deps_dalpha = fscanf(f_id,'%f '); temp = fgetl(f_id);
                for i=1:4
                    temp = fgetl(f_id);
                end
                %% Command linkage characteristics
                % Irreversible Type
                obj.Rs_e = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.Rg_e = fscanf(f_id,'%f '); temp = fgetl(f_id);
                for i=1:1
                    temp = fgetl(f_id);
                end
                % Limitasions: angular excursion range
                obj.delta_e_max = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.delta_e_min = fscanf(f_id,'%f '); temp = fgetl(f_id);
                for i=1:2
                    temp = fgetl(f_id);
                end
                %% Propulsion daqta
                obj.T = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.T = obj.T*obj.g;
                obj.Cm_T_0 = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.Cm_T_alpha = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.mu_T = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.e_T = fscanf(f_id,'%f '); temp = fgetl(f_id);
                for i=1:2
                    temp = fgetl(f_id);
                end
                %% Limitations (aerodynamic and structural)
                obj.CL_max = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.CL_min = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.n_max = fscanf(f_id,'%f ');temp = fgetl(f_id);
                obj.n_min = fscanf(f_id,'%f '); temp = fgetl(f_id);
                for i=1:2
                    temp = fgetl(f_id);
                end
                %% Limitations on piloting stick force
                obj.Fe_max = fscanf(f_id,'%f '); temp = fgetl(f_id);
                obj.Fe_min = fscanf(f_id,'%f '); temp = fgetl(f_id);
                %% finally, set the error tag
                obj.err = 0;
            end % if f_id not -1            
        end % constructor, init. from file
    end % methods
    
end

