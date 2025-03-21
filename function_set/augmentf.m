classdef augmentf < handle
    properties
        D        % 
        R        % data after imputed
        M       % the measurement data, sampled every second
    end
    properties (Access = private)
        driving_  % Similar driver data
    end
    methods
        function obj=augmentf(D)
            obj.D=D; % 
            obj.initial; % initialization
            obj.vagment; % Our imputation model
            obj.linear; % Linear interpolation model
            obj.spline; % Spline interpolation model
            obj.speed_prcess; % BBI and markov-based interpolation model
            obj.SOCprocess; % Battery data interpolation model
            
        end
        function linear(obj)
            data = obj.D.data_trip;
            x1 = cumsum(data(:,1)); y1 = data(:,6);
            xq = 0:obj.D.ds:x1(end);
            v_linear = interp1(x1, y1, xq);
            v_linear(v_linear<0)=0;
            obj.R.v_linear = v_linear;
        end

        function initial(obj)
            Mm = obj.D.data_trip1;
            obj.M.v = Mm(:,6)';
            obj.M.t = 0:1:length(obj.M.v)-1;
            obj.driving_.drlim = obj.D.drlim;
            obj.driving_.fitf = obj.D.F2;
            obj.driving_.soc1e = obj.D.soc1_toenegy;
        end

        function spline(obj)
            data = obj.D.data_trip;
            x1 = cumsum(data(:,1)); y1 = data(:,6);
            xq = 0:obj.D.ds:x1(end);
            v_spline = spline(x1, y1, xq);
            v_spline(v_spline<0)=0;
            obj.R.v_spline = v_spline;
        end

        function vagment(obj)
            data = obj.D.data_trip;
            id = obj.D.id; ds = obj.D.ds; 
            index13 = data(:, 13) > 0; index16 = data(:, 16) < 0; a_index = find(index13 & index16); %Accelerator pedal inconsistent with positive acceleration
            index14 = data(:, 13) < 0; index16 = data(:, 16) > 0; b_index = find(index14 & index16); %Brake pedal inconsistent with negative acceleration
            ss = [0];
            vt = [];
            chindex = [];
            a_base = [];
            aa_max = max(data(:,16)); ab_max = min(data(:,16));
            v_max = max(data(:,6));
            sv = [diff(data(:,17)); 0]; %Ture distance
            pedal = data(:,13);
            %% Activate the calibration model
            for i = 1:length(data)-1
                sm = sv(i);
                sd = sum(ss);
                k = data(i+1,1);
                st = sm - sd;
                vc = data(i,6); v_0 = vc; % Current vehicle speed
                vcni = data(i+1,6); v_end = vcni; % Next moment speed
                vrr = [vc]; vcnir = [vcni];
                vr = [vc];
                a1 = data(i,16);
                aat = data(i+1,16);
                if heaviside(pedal(i+1)*aat)==1
                    a2 = aat;
                else
                    a2 = -aat;
                end
                if i>1
                    a0 = a_copy;
                else
                    a0 = 0;
                end
                %% Special situations
                if a1<0 && a0<=0 && a1>0.8*a0 && a2>=0 && data(i+1,6)<1.5 % Model7
                    beta = 0.2; arf = 1.1;
                    W2 = arf.^(k-1:-1:0);
                    a1 = 4/k*a0.*W2;
                    for j = 1:floor(data(i+1,1)/ds)
                        vv = vc + a1(j)*ds;
                        vrr = [vrr vv];
                        vc = vrr(end);
                    end
                    vcnir = data(i+1,6).*ones(1,length(vrr));
                    [m indexv] = min(abs(vrr-vcnir));
                    vr = [vrr(1:indexv-1) vcnir(indexv:end)];
                    vr(vr<0)=0;
                    chindex = [chindex;i];
                elseif a1>=0 && a0<=0 && data(i,6)<1 % Model8
                    a1 = -aa_max.*ones(1,10);
                    vc = data(i+1,6);
                    vrr = vc;
                    for j = 1:floor(data(i+1,1)/ds)
                        vv = vc + a1(j)*ds;
                        vrr = [vrr vv];
                        vc = vrr(end);
                    end
                    vcnir = data(i,6).*ones(1,length(vrr));
                    [m indexv] = min(abs(vrr-vcnir));
                    vr = flip([vrr(1:indexv) vcnir(indexv+1:end)]);
                    vr(vr<0)=0;
                    chindex = [chindex;i];
                    %% Dual-mode Scenario
                elseif heaviside(-a0*a1)==1 && heaviside(-a0*a2)==1 && (any(find(a_index==i)) || any(find(b_index==i))) % Model1
                    myConstraints = @(x) constraintFunction2(x, obj.driving_.drlim);
                    objectiveFcn = @(x) optv2(sv(i), ds, x(1), x(2), a0, vc, vcni, vrr, vcnir, aa_max, k);
                    initial_guess = [0.95, 0.95];  %
                    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');
                    [x_optimized, fval] = fmincon(objectiveFcn, initial_guess, [], [], [], [], [], [], myConstraints, options);
                    disp('Optimized x1:');
                    disp(x_optimized);
                    disp('Function value:');
                    disp(fval);
                    W1 = x_optimized(1) - (1.01.^(1:1:k)-1);
                    W2 = x_optimized(2) - (1.01.^(1:1:k)-1);
                    a1 = a0.*W1;
                    for j = 1:floor(k/ds)
                        vv = vc + a1(j)*ds;
                        vrr = [vrr vv];
                        vc = vrr(end);
                    end
                    if a0<0
                        a2 = a0.*W2;
                    else
                        a2 = aa_max.*W2;
                    end
                    for j = 1:floor(k/ds)
                        vv = vcni + a2(j)*ds;
                        vcnir = [vcnir vv];
                        vcni = vcnir(end);
                    end
                    vcnir = flip(vcnir);
                    [m indexv] = min(abs(vrr-vcnir));
                    vr = [vrr(1:indexv) vcnir(indexv+1:end)];
                    vr(vr<0)=0;
                    chindex = [chindex;i];
                elseif heaviside(a0*a1)==1 && heaviside(-a0*a2)==1 && ~(any(find(a_index==i)) || any(find(b_index==i)))  %Model3
                    myConstraints = @(x) constraintFunction2(x, obj.driving_.drlim);
                    objectiveFcn = @(x) optv2(sv(i), ds, x(1), x(2), a0, vc, vcni, vrr, vcnir, aa_max, k);
                    initial_guess = [0.95, 0.95];  %
                    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');
                    [x_optimized, fval] = fmincon(objectiveFcn, initial_guess, [], [], [], [], [], [], myConstraints, options);
                    disp('Optimized x1:');
                    disp(x_optimized);
                    disp('Function value:');
                    disp(fval);
                    W1 = x_optimized(1) - (1.01.^(1:1:k)-1);
                    W2 = x_optimized(2) - (1.01.^(1:1:k)-1);
                    a1 = a0.*W1;
                    for j = 1:floor(k/ds)
                        vv = vc + a1(j)*ds;
                        vrr = [vrr vv];
                        vc = vrr(end);
                    end
                    if a0<0
                        a2 = a0.*W2;
                    else
                        a2 = aa_max.*W2;
                    end
                    for j = 1:floor(k/ds)
                        vv = vcni + a2(j)*ds;
                        vcnir = [vcnir vv];
                        vcni = vcnir(end);
                    end
                    vcnir = flip(vcnir);
                    [m indexv] = min(abs(vrr-vcnir));
                    vr = [vrr(1:indexv) vcnir(indexv+1:end)];
                    vr(vr<0)=0;
                    chindex = [chindex;i];
                elseif heaviside(a0*a1)==1 && heaviside(a0*a2)==1 && (any(find(a_index==i)) || any(find(b_index==i))) %Model5
                    myConstraints = @(x) constraintFunction2(x, obj.driving_.drlim);
                    objectiveFcn = @(x) optv3(sv(i), ds, x(1), x(2), a0, vc, vcni, vrr, vcnir, aa_max, k, data(i,16), a2);
                    initial_guess = [0.95, 0.95];  %
                    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');
                    [x_optimized, fval] = fmincon(objectiveFcn, initial_guess, [], [], [], [], [], [], myConstraints, options);
                    disp('Optimized x1:');
                    disp(x_optimized);
                    disp('Function value:');
                    disp(fval);
                    W1 = x_optimized(1) - (1.01.^(1:1:k)-1);
                    W2 = x_optimized(2) - (1.01.^(1:1:k)-1);
                    a1 = -a0.*W1;
                    for j = 1:floor(data(i+1,1)/ds)
                        vv = vc + a1(j)*ds;
                        vrr = [vrr vv];
                        vc = vrr(end);
                    end
                    if -a0<0
                        a2 = min([-a0,-data(i,16),-a2]).*W2;
                    else
                        a2 = aa_max.*W2;
                    end
                    for j = 1:floor(data(i+1,1)/ds)
                        vv = vcni + a2(j)*ds;
                        vcnir = [vcnir vv];
                        vcni = vcnir(end);
                    end
                    vcnir = flip(vcnir);
                    [m indexv] = min(abs(vrr-vcnir));
                    vr = [vrr(1:indexv) vcnir(indexv+1:end)];
                    vr(vr<0)=0;
                    chindex = [chindex;i];
                %% Tri-mode Scenario
                elseif heaviside(-a0*a1)==1 && heaviside(a0*a2)==1 && (any(find(a_index==i)) || any(find(b_index==i))) %Model2
                    myConstraints = @(x) constraintFunction(x, vc, vcni, obj.driving_.drlim, a0, a2);
                    objectiveFcn = @(x) optv(sv(i), ds, x(1), x(2), a0, a2, vc, vcni);
                    initial_guess = [1, 3];
                    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');
                    [x_optimized, fval] = fmincon(objectiveFcn, initial_guess, [], [], [], [], [], [], myConstraints, options);
                    disp('Optimized x3 and x4:');
                    disp(x_optimized);
                    disp('Function value:');
                    disp(fval);
                    n1 = [0; vc];
                    n2 = [10; vcni];
                    nc1 = [x_optimized(1); vc+x_optimized(1)*a0];
                    nc2 = [x_optimized(1)+x_optimized(2); vcni-a2*(10-x_optimized(1)-x_optimized(2))];
                    A = [n1(1)^3 n1(1)^2 n1(1)^1 1; n2(1)^3 n2(1)^2 n2(1)^1 1;...;
                        nc1(1)^3 nc1(1)^2 nc1(1)^1 1; nc2(1)^3 nc2(1)^2 nc2(1)^1 1];
                    Y =  [n1(2); n2(2); nc1(2); nc2(2)];
                    cofs = inv(A)*Y;
                    ac = @(xx)cofs(1)*xx.^3 + cofs(2)*xx.^2 + cofs(3)*xx.^1 + cofs(4);
                    x1 = 0:1:10;
                    vr = ac(x1);
                    vr(vr<0)=0;
                    chindex = [chindex;i];
                elseif heaviside(a0*a1)==1 && heaviside(a0*a2)==1 && ~(any(find(a_index==i)) || any(find(b_index==i))) %Model4
                    myConstraints = @(x) constraintFunction(x, vc, vcni, obj.driving_.drlim, a0, a2);
                    objectiveFcn = @(x) optv(sv(i), ds, x(1), x(2), a0, a2, vc, vcni);
                    initial_guess = [1, 3];
                    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');
                    [x_optimized, fval] = fmincon(objectiveFcn, initial_guess, [], [], [], [], [], [], myConstraints, options);
                    disp('Optimized x3 and x4:');
                    disp(x_optimized);
                    disp('Function value:');
                    disp(fval);
                    n1 = [0; vc];
                    n2 = [10; vcni];
                    nc1 = [x_optimized(1); vc+x_optimized(1)*a0];
                    nc2 = [x_optimized(1)+x_optimized(2); vcni-a2*(10-x_optimized(1)-x_optimized(2))];
                    A = [n1(1)^3 n1(1)^2 n1(1)^1 1; n2(1)^3 n2(1)^2 n2(1)^1 1;...;
                        nc1(1)^3 nc1(1)^2 nc1(1)^1 1; nc2(1)^3 nc2(1)^2 nc2(1)^1 1];
                    Y =  [n1(2); n2(2); nc1(2); nc2(2)];
                    cofs = inv(A)*Y;
                    ac = @(xx)cofs(1)*xx.^3 + cofs(2)*xx.^2 + cofs(3)*xx.^1 + cofs(4);
                    x1 = 0:1:10;
                    vr = ac(x1);
                    vr(vr<0)=0;
                    chindex = [chindex;i];
                elseif heaviside(a0*a1)==1 && heaviside(-a0*a2)==1 && (any(find(a_index==i)) || any(find(b_index==i))) %Model6
                    myConstraints = @(x) constraintFunction(x, vc, vcni, obj.driving_.drlim, a0, a2);
                    objectiveFcn = @(x) optv(sv(i), ds, x(1), x(2), -a0, a2, vc, vcni);
                    initial_guess = [1, 3];
                    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');
                    [x_optimized, fval] = fmincon(objectiveFcn, initial_guess, [], [], [], [], [], [], myConstraints, options);
                    disp('Optimized x3 and x4:');
                    disp(x_optimized);
                    disp('Function value:');
                    disp(fval);
                    n1 = [0; vc];
                    n2 = [10; vcni];
                    nc1 = [x_optimized(1); vc+x_optimized(1)*-a0];
                    nc2 = [x_optimized(1)+x_optimized(2); vcni-a2*(10-x_optimized(1)-x_optimized(2))];
                    A = [n1(1)^3 n1(1)^2 n1(1)^1 1; n2(1)^3 n2(1)^2 n2(1)^1 1;...;
                        nc1(1)^3 nc1(1)^2 nc1(1)^1 1; nc2(1)^3 nc2(1)^2 nc2(1)^1 1];
                    Y =  [n1(2); n2(2); nc1(2); nc2(2)];
                    cofs = inv(A)*Y;
                    ac = @(xx)cofs(1)*xx.^3 + cofs(2)*xx.^2 + cofs(3)*xx.^1 + cofs(4);
                    x1 = 0:1:10;
                    vr = ac(x1);
                    vr(vr<0)=0;
                    chindex = [chindex;i];
                else
                    for j = 1:floor(data(i+1,1)/ds)
                        vv = vc + a1*ds;
                        vr = [vr vv];
                        vc = vr(end);
                    end
                    vr(vr<0)=0;
                end
                vr(1) = v_0; vr(end) = v_end;
                if any(diff(vr)>3)
                    vr = linspace(v_0, v_end, 11);
                end
                srt = 0:ds:k;
                sr = trapz(srt, vr);
                ss = [ss;sr];
                if i>1
                    vt = [vt vr(2:end)];
                else
                    vt = [vt vr];
                end
                a_copy = (vt(end)-vt(end-1))./ds;
                a_base = [a_base;a_copy];
            end
            obj.R.v_prop = vt; 
        end

        function speed_prcess(obj)
            TPM = creattpm(obj.M.v); % Generate TPM matrix using high frequency data 
            interl = 1:10:length(obj.M.v);
            t_low = obj.M.t(interl);
            v_low = obj.M.v(interl);
            v_highref = obj.M.v;
            a_limit = diff(obj.M.v);
            s_limit = [0 (obj.M.v(1:end-1) + obj.M.v(2:end))/2];
            t_low_res = t_low;
            v_low_res = v_low;
            % Call BBI function for interpolation
            [t_high_res, v_high_res] = bbi_imputation(t_low_res, v_low_res, a_limit, v_highref);
            % Call markov function for interpolation
            [t_high_resm, v_high_resm] = markov_imputation(t_low_res, v_low_res, a_limit, v_highref, TPM, obj.driving_.drlim);
            obj.R.v_bbi = v_high_res;
            obj.R.v_markov = v_high_resm;
        end
        function SOCprocess(obj)
            I_ture1 = obj.D.data_trip1(:,9)';
            U_ture1 = obj.D.data_trip1(:,8)';
            SOC_10 = obj.D.data_trip(:,10)';
            a_prop = [diff(obj.R.v_prop) 0];
            I_prop_fit = obj.driving_.fitf(obj.R.v_prop', a_prop');
            [I_1_unique, ~, idx] = unique(I_ture1');
            U_1_agg = accumarray(idx, U_ture1', [], @mean);
            U_prop_fit = interp1(I_1_unique', U_1_agg', I_prop_fit, 'linear', 'extrap');
            SOC_fit = zeros(length(obj.R.v_prop),1);
            x = cumsum(obj.D.data_trip(:,1))+1;
            for i = 1:length(x) - 1
                SOC_fit(x(i)) = obj.D.data_trip(i,10);
            end
            SOC_fit(end) = obj.D.data_trip(end,10);
            soc_deg_end = [10*find(diff(SOC_10)) length(obj.R.v_prop)];
            soc_deg_start = [1, 10*find(diff(SOC_10))+1];
            soc_deg_index = [soc_deg_start; soc_deg_end];
            soc1_toenegy = obj.driving_.soc1e;  % Capacity (As) corresponding to each 1% SOC
            % Processing of complete unit soc battery data
            soc = [];
            ad_parameters = [];
            for i = 2:length(soc_deg_start)-1
                istart = soc_deg_start(i);
                iend = soc_deg_end(i);
                Iad = I_prop_fit(istart:iend);
                tad = 0:1:length(Iad)-1;
                k_initial = 1;
                options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');
                lb = -2;
                ub = 2;
                [k_optimal, fval] = fmincon(@(k) objectiveFunction(k, tad, Iad, soc1_toenegy), k_initial, [], [], [], [], lb, ub, @(k) nonlincon(k, Iad), options);
                fprintf('Optimal k: %f\n', k_optimal);
                fprintf('Minimum soc1_toenergy: %f\n', fval);
                ad_parameters = [ad_parameters k_optimal];
                I_prop_fit(istart:iend) = k_optimal.*Iad;
                cs_current = cumsum(Iad);
                try
                    soc_seg = interp1([cs_current(1),cs_current(end)],[SOC_fit(istart),SOC_fit(istart)-1],cs_current, 'linear', 'extrap');
                catch ME
                    disp(cs_current)
                    disp(i)
                end
                soc = [soc; soc_seg];
            end
            % Processing battery data for the first part of the trip
            current = I_prop_fit(soc_deg_start(1):soc_deg_end(1));
            I_nagetive = current(current<0);
            ttt = 0:1:length(current)-1;
            energy_ofb1soc = cumtrapz(ttt,current);
            if energy_ofb1soc(end)>soc1_toenegy
                error = soc1_toenegy - energy_ofb1soc(end);
                ad_paramter = error/trapz(current(current<0));
                current(current<0) = I_nagetive.*(1+ad_paramter);
                I_prop_fit(soc_deg_start(1):soc_deg_end(1)) = current;
                cs_current = cumsum(current);
                soc_1b = interp1([cs_current(1),cs_current(end)],[SOC_fit(soc_deg_start(1)),SOC_fit(soc_deg_end(1)+1)],cs_current, 'linear', 'extrap');
            else
                percent1_soc = energy_ofb1soc./soc1_toenegy;
                soc_1b = flip(SOC_fit(soc_deg_end(1)+1) + percent1_soc);
            end

            % Processing battery data at the end of the trip
            current = I_prop_fit(soc_deg_start(end):soc_deg_end(end));
            I_nagetive = current(current<0);
            ttt = 0:1:length(current)-1;
            energy_ofa1soc = cumtrapz(ttt,current);
            if energy_ofa1soc(end)>soc1_toenegy
                error = soc1_toenegy - energy_ofa1soc(end);
                ad_paramter = error/trapz(current(current<0));
                current(current<0) = I_nagetive.*(1+ad_paramter);
                I_prop_fit(soc_deg_start(end):soc_deg_end(end)) = current;
                cs_current = cumsum(current);
                soc_1a = interp1([cs_current(1),cs_current(end)],[SOC_fit(soc_deg_start(end)),SOC_fit(soc_deg_end(end))-1],cs_current, 'linear', 'extrap');
            else
                percent1_soc = energy_ofa1soc./soc1_toenegy;
                soc_1a = SOC_fit(soc_deg_end(end)) - percent1_soc;
            end

            soc = [soc_1b; soc; soc_1a];
            soc_k = SOC_fit;
            %% baseline
            op_I = [];
            for i = 2:length(soc_deg_start)-1
                n = soc_deg_index(2,i) - soc_deg_index(1,i) + 1;
                initial_values = zeros(n, 1);
                objective = @(x) soc1_toenegy - sum(x);
                Aeq = ones(1, n);
                beq = soc1_toenegy;
                lb = ones(n, 1)*-100;
                ub = ones(n, 1)*150;
                options = optimoptions('fmincon', 'Display', 'iter'); 
                optimized_I = fmincon(objective, initial_values, [], [], Aeq, beq, lb, ub, [], options);
                op_I = [op_I; optimized_I];
            end
            I_end = I_prop_fit(soc_deg_start(end):soc_deg_end(end));
            I_0 = I_prop_fit(soc_deg_start(1):soc_deg_end(1));
            I_based = [I_0; op_I; I_end];
            U_based = mean(U_prop_fit)*ones(length(I_based),1);
            obj.R.I = I_prop_fit; obj.R.U = U_prop_fit; obj.R.SOC = soc;
            obj.R.Ib = I_based; obj.R.Ub = U_based;
            obj.M.I = I_ture1; obj.M.U = U_ture1; obj.M.SOC = soc_k;
        end
    end
end


function saveWorkspaceVars(filename)
    vars = evalin('base', 'who');  
    dataR = struct();  
    for i = 1:length(vars)
        varName = vars{i};
        dataR.(varName) = evalin('base', varName); 
    end
    save(filename, '-struct', 'dataR');  
end

function [t_high_res, v_high_res] = bbi_imputation(t_low_res, v_low_res, a_limit, v_highref)
q = max(a_limit); 
sigma = 0.3; 
alpha = 0.5; 
epsilon = 0.15; 
t_high_res = 0:1:max(t_low_res); 
v_high_res = zeros(size(t_high_res)); 
for i = 1:length(t_low_res)-1
    t_start = t_low_res(i);
    t_end = t_low_res(i+1);
    v_start = v_low_res(i);
    v_end = v_low_res(i+1);
    if v_start == v_end && v_start == 0
        t_indices = t_start+1:t_end;
        v_high_res(t_indices) = zeros(1,10);
    else
        delta_t = t_end - t_start;
        delta_d  = sum(v_highref(t_start+1: t_end));
        v_interval = zeros(1, delta_t);
        v_interval(1) = v_start;
        % Brownian Bridge Interpolation Process
        for tau = 2:delta_t
            % Calculating convergence drift
            lambda = (v_end - v_interval(tau-1)) / (delta_t - tau + 2);
            % Calculate distance drift
            delta_d_estimated = sum(v_interval(1:tau-1)) + (v_end + v_interval(tau-1)) / 2 * (delta_t - tau + 1);
            eta = alpha * (delta_d - delta_d_estimated) / (delta_d + delta_d_estimated);
            % Combination drift
            mu = lambda + eta;
            % Add gaussian noise
            r = sigma * rand;
            % Calculating acceleration
            a = mu + r; 
            a = max(min(a, q), -q);
            % Update speed
            v_interval(tau) = v_interval(tau-1) + a;
            % Applying a speed ​​constraint
            v_interval(tau) = max([min([v_interval(tau), v_start + q * tau, v_end + q * (delta_t - tau)]), ...
                min([v_interval(tau), v_start - q * tau, v_end - q * (delta_t - tau)]), 0]);
        end
        t_indices = t_start+1:t_end;
        v_high_res(t_indices) = v_interval;
    end
end
end

function [t_high_res, v_high_res] = markov_imputation(t_low_res, v_low_res, a_limit, v_highref, TPM, drlim)
acceleration_sequence = [diff(v_highref) 0];
speed_bins = floor(min(v_highref)):1:ceil(max(v_highref));
acceleration_bins = floor(min(acceleration_sequence)):0.1:ceil(max(acceleration_sequence));
q = max(a_limit); 
t_high_res = 0:1:max(t_low_res); 
v_high_res = zeros(size(t_high_res)); 
for i = 1:length(t_low_res)-1
    t_start = t_low_res(i);
    t_end = t_low_res(i+1);
    v_start = v_low_res(i);
    v_end = v_low_res(i+1);
    if v_start == v_end && v_start == 0
        t_indices = t_start+1:t_end;
        v_high_res(t_indices) = zeros(1,10);
    else
        delta_t = t_end - t_start;
        a_asist = (v_end-v_start)./delta_t;
        v_interval = zeros(1, delta_t);
        v_interval(1) = v_start;

        for tau = 2:delta_t
            speed_idx = find(v_interval(tau-1) >= speed_bins(1:end-1) & v_interval(tau-1) < speed_bins(2:end));
            TPM_row = TPM(speed_idx, :);
            cumulative_prob = cumsum(TPM_row);
            rand_num = rand();
            selected_acceleration_idx = find(cumulative_prob >= rand_num, 1);
            a = acceleration_bins(selected_acceleration_idx) + a_asist;
            a = max(min(a, q), -q);
            v_interval(tau) = v_interval(tau-1) + a;
            v_interval(tau) = max([min([v_interval(tau), v_start + q * tau, v_end + q * (delta_t - tau)]), ...
                min([v_interval(tau), v_start - q * tau, v_end - q * (delta_t - tau)]), 0]);
        end
        t_indices = t_start+1:t_end;
        v_high_res(t_indices) = v_interval;
    end

    
end
end

function obvalue = objectiveFunction(k, tad, Iad, soc1_toenegy)
    obvalue = abs(soc1_toenegy - trapz(tad, k .* Iad));
end

function [c, ceq] = nonlincon(k, Iad)
    product = k * Iad;
    c = [product - 200;  
        -100 - product]; 
    ceq = [];
end

function [c, ceq] = constraintFunction(x, vc, vcni, drlim, a0, a2)
    nc1 = [x(1); vc+x(1)*a0];
    nc2 = [x(1)+x(2); vcni-a2*(10-x(1)-x(2))];
    k_c = (nc2(2) - nc1(2))./ (nc2(1) - nc1(1));
    c(1) = x(1) - 9;
    c(2) = x(2) - 10 + x(1);
    c(3) = -x(1) + 2;
    c(4) = -x(2) + 5;
    c(5) = k_c - drlim.au;
    c(6) = -k_c +  drlim.ad;
    ceq = [];
end

function [c, ceq] = constraintFunction2(x, drlim)
    c(1) = x(1) - drlim.bataau;
    c(2) = -x(1) + drlim.bataad;
    c(3) = x(2) - drlim.batadu;
    c(4) = -x(2) + drlim.batadd;
    ceq = [];
end
