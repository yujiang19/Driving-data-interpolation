classdef valdationpolt < handle
    properties
        data 
        fea
        ferr
    end
    
    methods
        function obj = valdationpolt(data)
            obj.data = data;
            obj.fea = table();
            obj.ferr = table();
        end
        function calfea(obj)
            RowNames = {'Maximum driving speed', 'Average speed', 'Average driving speed', 'Speed ​​standard deviation', 'Average acceleration',...
            'Acceleration standard deviation', 'Average positive acceleration', 'Maximum positive acceleration', 'Maximum negative acceleration',...
            'Average negative acceleration', 'Distance','Acceleration ratio','Deceleration ratio','Idle speed ratio'};
            x_labels = {'v_{max}', 'v_{m}', 'v_{mt}', '\sigma_{v}', 'a_{m}', ...
            '\sigma_{a}', 'a_{am}', 'a_{amax}', 'a_{dmax}', 'a_{dm}', 'd_{t}', ...
            'P_{a}', 'P_{d}', 'P_{i}'};
            features_linear = cal_features(obj.data.M.t, obj.data.R.v_linear);
            features_prop = cal_features(obj.data.M.t, obj.data.R.v_prop);
            features_1s = cal_features(obj.data.M.t, obj.data.M.v);
            features_spline = cal_features(obj.data.M.t, obj.data.R.v_spline);
            features_v3 = cal_features(obj.data.M.t, obj.data.R.v_bbi);
            features_markov = cal_features(obj.data.M.t, obj.data.R.v_markov);
            linear_err = err_features(features_1s, features_linear);
            spline_err = err_features(features_1s, features_spline);
            v3_err = err_features(features_1s, features_v3);
            markov_err = err_features(features_1s, features_markov);
            prop_err = err_features(features_1s, features_prop);
            combined = [features_1s, features_linear, features_spline, features_v3, features_markov, features_prop];
            combined_err = [linear_err, spline_err, v3_err, markov_err, prop_err]; obj.ferr = array2table(combined_err, 'RowNames', RowNames);
            obj.fea = array2table(combined, 'RowNames', RowNames);
            obj.fea.Properties.VariableNames = {'features_1s', 'features_linear', 'features_spline', 'features_bbi', 'features_markov', 'features_prop'};
            obj.ferr.Properties.VariableNames = {'features_linear', 'features_spline', 'features_bbi', 'features_markov', 'features_prop'};
            figure(98);
            hold on
            plot(1:length(x_labels), linear_err', 'b-d', 'DisplayName', 'Baseline approach1', 'LineWidth', 2.5, 'MarkerSize', 10);
            plot(1:length(x_labels), spline_err', 'g-d', 'DisplayName', 'Baseline approach2', 'LineWidth', 2.5, 'MarkerSize', 10);
            plot(1:length(x_labels), v3_err', 'r-^', 'DisplayName', 'Baseline approach3', 'LineWidth', 2.5, 'MarkerSize', 10);
            plot(1:length(x_labels), markov_err', 'Color', [0.6, 0.3, 0.1], 'Marker', '^', 'DisplayName', 'Baseline approach4', 'LineWidth', 2.5, 'MarkerSize', 10);
            plot(1:length(x_labels), prop_err', 'o-', 'Color', [1, 0.5, 0], 'DisplayName', 'Proposed approach', 'LineWidth', 2.5, 'MarkerSize', 10);
            grid on;
            title('Comparison of the error between the reconstructed speed and the original speed in characteristic parameters using different methods.')
            xlabel('Parameter');
            ylabel('Error [%]');
            legend('Baseline method1', 'Baseline method2', 'Baseline method3', 'Baseline method4', 'Proposed method');
            % [linear_method spline_method bbi_method markov_method proposed_method]
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 24);
        end
        function vplot(obj)
            t_1 = obj.data.M.t;
            v_ture1 = obj.data.M.v;
            v_linear = obj.data.R.v_linear;
            v_v3 = obj.data.R.v_bbi;
            v_markov = obj.data.R.v_markov;
            vv_spline = obj.data.R.v_spline;
            vt = obj.data.R.v_prop;
            colors = {[0.871, 0.196, 0.298], [0.957, 0.537, 0.373], [0.973, 0.882, 0.435],... 
                [0.588, 0.337, 0.635], [0.584, 0.812, 0.573], [0.212, 0.604, 0.8]};
            figure(1);
            subplot(3,1,1);
            hold on;
            plot(t_1, v_ture1, 'LineWidth', 2, 'Color', colors{1});  % 
            plot(t_1, v_linear, 'LineWidth', 2, 'Color', colors{2});  % 
            plot(t_1, vv_spline, 'LineWidth', 2, 'Color', colors{3});  % 
            plot(t_1, v_v3, 'LineWidth', 2, 'Color', colors{5});  % 
            plot(t_1, v_markov, 'LineWidth', 2, 'Color', colors{6});  % 
            plot(t_1, vt, 'LineWidth', 2, 'Color', colors{4});  %
            xlabel('Time (s)');
            ylabel('Velocity (m/s)');
            legend('Original data','Linear method','Spline method','Gassi method','Markov method','Proposed method');
            title('Speed reconstruction results under different methods.')
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
            zoom_x1 = 400:500;  zoom_x2 = 200:290; zoom_x3 = 1014:1156;
            subplot(3,3,4); 
            hold on;
            plot(t_1(zoom_x1), v_ture1(zoom_x1), 'LineWidth', 2, 'Color', colors{1}); %
            plot(t_1(zoom_x1), v_linear(zoom_x1), 'LineWidth', 2, 'Color', colors{2}); %
            plot(t_1(zoom_x1), vv_spline(zoom_x1), 'LineWidth', 2, 'Color', colors{3}); %
            plot(t_1(zoom_x1), vt(zoom_x1), 'LineWidth', 2, 'Color', colors{4}); %
            ylabel('Velocity (m/s)');
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
            subplot(3,3,5); 
            hold on;
            plot(t_1(zoom_x2), v_ture1(zoom_x2), 'LineWidth', 2, 'Color', colors{1}); % 
            plot(t_1(zoom_x2), v_linear(zoom_x2), 'LineWidth', 2, 'Color', colors{2}); % 
            plot(t_1(zoom_x2), vv_spline(zoom_x2), 'LineWidth', 2, 'Color', colors{3}); % 
            plot(t_1(zoom_x2), vt(zoom_x2), 'LineWidth', 2, 'Color', colors{4}); % 
            ylabel('Velocity (m/s)');
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
            subplot(3,3,6); 
            hold on;
            plot(t_1(zoom_x3), v_ture1(zoom_x3), 'LineWidth', 2, 'Color', colors{1}); % 
            plot(t_1(zoom_x3), v_linear(zoom_x3), 'LineWidth', 2, 'Color', colors{2}); % 
            plot(t_1(zoom_x3), vv_spline(zoom_x3), 'LineWidth', 2, 'Color', colors{3}); % 
            plot(t_1(zoom_x3), vt(zoom_x3), 'LineWidth', 2, 'Color', colors{4}); % 
            ylabel('Velocity (m/s)');
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
            subplot(3,3,7);
            hold on;
            plot(t_1(zoom_x1), v_ture1(zoom_x1), 'LineWidth', 2, 'Color', colors{1}); % 
            plot(t_1(zoom_x1), v_v3(zoom_x1), 'LineWidth', 2, 'Color', colors{5}); % 
            plot(t_1(zoom_x1), v_markov(zoom_x1), 'LineWidth', 2, 'Color', colors{6}); % 
            plot(t_1(zoom_x1), vt(zoom_x1), 'LineWidth', 2, 'Color', colors{4}); % 
            xlabel('Time (s)');
            ylabel('Velocity (m/s)');
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
            subplot(3,3,8);
            hold on;
            plot(t_1(zoom_x2), v_ture1(zoom_x2), 'LineWidth', 2, 'Color', colors{1}); % 
            plot(t_1(zoom_x2), v_v3(zoom_x2), 'LineWidth', 2, 'Color', colors{5}); % 
            plot(t_1(zoom_x2), v_markov(zoom_x2), 'LineWidth', 2, 'Color', colors{6}); % 
            plot(t_1(zoom_x2), vt(zoom_x2), 'LineWidth', 2, 'Color', colors{4}); % 
            xlabel('Time (s)');
            ylabel('Velocity (m/s)');
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
            subplot(3,3,9);
            hold on;
            plot(t_1(zoom_x3), v_ture1(zoom_x3), 'LineWidth', 2, 'Color', colors{1}); % 
            plot(t_1(zoom_x3), v_v3(zoom_x3), 'LineWidth', 2, 'Color', colors{5}); % 
            plot(t_1(zoom_x3), v_markov(zoom_x3), 'LineWidth', 2, 'Color', colors{6}); % 
            plot(t_1(zoom_x3), vt(zoom_x3), 'LineWidth', 2, 'Color', colors{4}); % 
            xlabel('Time (s)');
            ylabel('Velocity (m/s)');
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 15);
        end
        function vaplot(obj)
            v_ture1 = obj.data.M.v;
            vt = obj.data.R.v_prop;
            v_linear = obj.data.R.v_linear;
            am = [0 diff(v_ture1)]; Xm = [v_ture1', am'];
            at = [0 diff(vt)]; Xt = [vt', at'];
            al = [0 diff(v_linear)]; Xl = [v_linear', al'];
            max_v = max(v_ture1);
            bar_width = 0.15; 
            speed_intervals = 0:2:ceil(max_v); 
            acceleration_intervals = -3:0.2:3; 
            figure(3);
            num_intervals = length(speed_intervals) - 1; 
            num_rows = ceil(num_intervals / 3); 
            num_cols = min(3, num_intervals); 

            for i = 1:num_intervals-1
                lower_speed = speed_intervals(i);
                upper_speed = speed_intervals(i + 1);
                Xm_in = find(Xm(:,1) >= lower_speed & Xm(:,1) < upper_speed);
                Xl_in = find(Xl(:,1) >= lower_speed & Xl(:,1) < upper_speed);
                Xt_in = find(Xt(:,1) >= lower_speed & Xt(:,1) < upper_speed);
                [f1, x1] = histcounts(Xm(Xm_in,2), acceleration_intervals, 'Normalization', 'probability');
                [f2, x2] = histcounts(Xl(Xl_in,2), acceleration_intervals, 'Normalization', 'probability');
                [f3, x3] = histcounts(Xt(Xt_in,2), acceleration_intervals, 'Normalization', 'probability');
                x_centers = (x1(1:end-1) + x1(2:end)) / 2; 
                subplot(num_rows, num_cols, i); 
                hold on;
                bar(x_centers, f1, 'FaceColor', 'b', 'FaceAlpha', 1, 'DisplayName', 'origin data');
                bar(x_centers, f2, 'FaceColor', 'r', 'FaceAlpha', 1, 'DisplayName', 'linear-based');
                bar(x_centers, f3, 'FaceColor', 'g', 'FaceAlpha', 1, 'DisplayName', 'ours')
                pd1 = fitdist(Xm(Xm_in, 2), 'Normal'); 
                pd2 = fitdist(Xl(Xl_in, 2), 'Normal');
                pd3 = fitdist(Xt(Xt_in, 2), 'Normal'); 
                x_fit = linspace(min(acceleration_intervals), max(acceleration_intervals), 100);
                y_fit1 = pdf(pd1, x_fit);
                y_fit2 = pdf(pd2, x_fit);
                y_fit3 = pdf(pd3, x_fit);
                plot(x_fit, y_fit1, 'b-', 'LineWidth', 2, 'DisplayName', 'Fit: Origin Data');
                plot(x_fit, y_fit2, 'r-', 'LineWidth', 2, 'DisplayName', 'Fit: Linear-based');
                plot(x_fit, y_fit3, 'g-', 'LineWidth', 2, 'DisplayName', 'Fit: Ours');
                title(sprintf('%d-%d m/s', lower_speed, upper_speed));
                legend('show');
                hold off;
                set(gca, 'FontName', 'Palatino Linotype', 'FontSize', 14);
            end
                xlabel('Acceleration (m/s^2)');
                ylabel('Fraction/f(x)');
            [x_1, f_1] = ksdensity(Xm(:,2));  
            [x_2, f_2] = ksdensity(Xl(:,2)); 
            [x_3, f_3] = ksdensity(Xt(:,2));  
            f_max = max([f_1(2)-f_1(1); f_2(2)-f_2(1); f_3(2)-f_3(1)]);
            entropy1 =  calentropy(x_1, f_1, f_max);
            entropy2 =  calentropy(x_2, f_2, f_max);
            entropy3 =  calentropy(x_3, f_3, f_max);
            disp('The entropy of driving data：');
            disp(['The entropy of 1 second raw data is：', num2str(entropy1)]);
            disp(['The entropy of baseline：', num2str(entropy2)]);
            disp(['The entropy of proposed method：', num2str(entropy3)]);
            figure(4);
            hold on;
            fill([f_1 fliplr(zeros(size(f_1)))], [x_1 fliplr(x_1)], [0, 0, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Original Data');
            fill([f_2 fliplr(zeros(size(f_2)))], [x_2 fliplr(x_2)], [1, 0, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Linear-based');
            fill([f_3 fliplr(zeros(size(f_3)))], [x_3 fliplr(x_3)], [0, 1, 0], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Ours');
            % 
            plot(f_1, x_1, 'Color', [0, 0, 1], 'LineWidth', 2);
            plot(f_2, x_2, 'Color', [1, 0, 0], 'LineWidth', 2);
            plot(f_3, x_3, 'Color', [0, 1, 0], 'LineWidth', 2);
            text(15, 0.3, '(number of original samples: 500)', 'Color', 'k', 'FontSize', 10, 'HorizontalAlignment', 'center');
            xlabel('Acceleration (m/s^2)');
            ylabel('Probability Density');
            title('The kernel density function curves of acceleration across the entire speed range.');
            legend('show', 'Location', 'northeast');
            set(gca, 'Color', [1 1 1]); 
            hold off;
            set(gca, 'FontName', 'Palatino Linotype', 'FontSize', 24);
            title('The kernel density function curves of acceleration across the entire speed range.')
        end
        function socplot(obj)
            %% SOC
            soc = obj.data.R.SOC;
            soc_k = obj.data.M.SOC;
            v_ture1 = obj.data.M.v;    
            n = length(soc_k);
            for i = 2:n
                if soc_k(i) == 0
                    soc_k(i) = soc_k(i-1); 
                end
            end
            figure(5);
            subplot(2,1,1);
            yyaxis left;
            plot(v_ture1, 'LineWidth', 1.5);
            ylabel('Velocity (m/s)');
            yyaxis right;
            title('SOC resampling results.')
            plot(soc, 'LineWidth', 1.5);
            hold on;
            plot(soc_k, 'LineWidth', 1.5);
            legend('Original Velocity', 'Proposed method', 'Original SOC');
            xlabel('Time (s)');
            ylabel('SOC');
            set(gca, 'FontName', 'Palatino Linotype', 'FontSize', 24);
            zoom_x1 = 332:538; 
            zoom_x2 = 740:915;
            subplot(2,2,3);
            yyaxis left;
            plot(v_ture1(zoom_x1), 'LineWidth', 1.5);
            ylabel('Velocity (m/s)');
            yyaxis right;
            plot(soc(zoom_x1), 'LineWidth', 1.5);
            legend('Original Velocity', 'Proposed method');
            xlabel('Time (s)');
            ylabel('SOC');
            set(gca, 'FontName', 'Palatino Linotype', 'FontSize', 16);
            subplot(2,2,4);
            yyaxis left;
            plot(v_ture1(zoom_x2), 'LineWidth', 1.5);
            ylabel('Velocity (m/s)');
            yyaxis right;
            plot(soc(zoom_x2), 'LineWidth', 1.5);
            legend('Original Velocity', 'Proposed method');
            xlabel('Time (s)');
            ylabel('SOC');
            set(gca, 'FontName', 'Palatino Linotype', 'FontSize', 16);
            
            %% current
            figure(6);
            subplot(1,2,1)
            I_ture1 = obj.data.M.I;
            I_prop_fit = obj.data.R.I;
            plot(I_ture1, 'b-', 'LineWidth', 1.5); %
            hold on;
            plot(I_prop_fit, 'r-', 'LineWidth', 1.5); %
            xlabel('Time (s)');
            ylabel('Current (A)');
            legend('Original current', 'Proposed current');
            % title('Resampling results of current and voltage under different methods.')
            set(gca, 'FontName', 'Palatino Linotype', 'FontSize', 24);
            ax1 = gca;
            grid on;
            bar_width = 1;
            bins = -100:20:200;
            [Iture_hist, ~] = histcounts(I_ture1, bins); %
            [Iprop, ~] = histcounts(I_prop_fit, bins); %
            Iture_ratio = Iture_hist / length(I_ture1); %
            Iprop_ratio = Iprop / length(I_prop_fit); %
            x_centers = (bins(1:end-1) + bins(2:end)) / 2;
            ax2=axes('Position',get(ax1,'Position'),...
                'XAxisLocation','top',...
                'YAxisLocation','right',...
                'Color','none',...
                'XColor','k','YColor','k');
            hold on;
            barh(x_centers, Iture_ratio, 'FaceColor', 'b', 'FaceAlpha', 0.2, 'BarWidth', bar_width, 'DisplayName', 'Original Data','Parent',ax2);
            barh(x_centers, Iprop_ratio, 'FaceColor', 'r', 'FaceAlpha', 0.2, 'BarWidth', 0.8, 'DisplayName', 'Proposed Based','Parent',ax2);
            grid on;
            legend('Original current', 'Proposed current');
            set(gca, 'FontName', 'Palatino Linotype', 'FontSize', 24); 
            xlabel('Proportion');
            %% Voltage
            subplot(1,2,2)
            U_ture1 = obj.data.M.U;
            U_prop_fit = obj.data.R.U;
            plot(U_ture1, 'b-', 'LineWidth', 1.5); 
            hold on;
            plot(U_prop_fit, 'r-', 'LineWidth', 1.5); 
            xlabel('Time (s)');
            ylabel('Voltage (V)');
            ax1 = gca;
            set(ax1, 'FontName', 'Palatino Linotype', 'FontSize', 24);
            grid on;
            bar_width = 1;
            bins = 275:5:315;
            [Uture_hist, ~] = histcounts(U_ture1, bins); %
            [Uprop, ~] = histcounts(U_prop_fit, bins); %
            Uture_ratio = Uture_hist / length(U_ture1); %
            Uprop_ratio = Uprop / length(U_prop_fit); %
            x_centers = (bins(1:end-1) + bins(2:end)) / 2;
            ax2=axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
            hold on;
            barh(x_centers, Uture_ratio, 'FaceColor', 'b', 'FaceAlpha', 0.2, 'BarWidth', bar_width, 'DisplayName', 'Original Data','Parent',ax2);
            barh(x_centers, Uprop_ratio, 'FaceColor', 'r', 'FaceAlpha', 0.2, 'BarWidth', 0.8, 'DisplayName', 'Proposed Based','Parent',ax2);
            grid on;
            legend('Original voltage', 'Proposed voltage');
            set(gca, 'FontName', 'Palatino Linotype', 'FontSize', 24); 
            xlabel('Proportion');
        end
        function energyplot(obj)
            energy_total = [];
            for ik = 1:length(obj.data.D)
                pdata = obj.data;
                if ~isempty(pdata)
                    odindex = splitt(pdata);
                else
                    continue;
                end
                energy_seg = [];
                for kk = 1:size(odindex,2)
                    start_index = odindex(1,kk); end_index = odindex(2,kk);
                    Iture_seg = pdata.M.I(start_index:end_index);
                    Uture_seg = pdata.M.U(start_index:end_index);
                    Ipre_seg = pdata.R.I(start_index:end_index);
                    Upre_seg = pdata.R.U(start_index:end_index);
                    Ibs_seg = pdata.R.Ib(start_index:end_index);
                    Ubs_seg = pdata.R.Ub(start_index:end_index);
                    t = 0:1:length(Upre_seg)-1;
                    energy_ture = trapz(t, Iture_seg.*Uture_seg);
                    energy_pre = trapz(t, Ipre_seg.*Upre_seg);
                    energy_bs = trapz(t, Ibs_seg.*Ubs_seg);
                    energy_seg = [energy_seg [energy_ture;energy_pre;energy_bs]];
                end
                energy_total = [energy_total energy_seg];
            end
            energy_totalkwh = energy_total/3600/1000;
            samples = 1:1:length(energy_totalkwh);
            y = energy_totalkwh(1, :); y_predl = energy_totalkwh(3, :); y_predp = energy_totalkwh(2, :);
            error_linear = 100.*abs(y - y_predl)./abs(y);
            error_prop = 100.*abs(y - y_predp)./abs(y);
            num_seg = max(1,floor(length(error_prop)/15)); index_seg = 1:num_seg:length(error_prop);
            prop_sef = [];
            for ij = 1:length(index_seg)-1
                seg = error_prop(index_seg(ij):index_seg(ij+1) -1)';
                prop_sef = [prop_sef seg];
            end
            % Error scatter plot
            figure(23);
            scatter(y, y_predp, 'o', 'DisplayName', 'Proposed method');
            hold on;
            scatter(y, y_predl, 'x', 'DisplayName', 'Baseline');
            xlabel('Measurements (kWh)', 'FontName', 'Times New Roman', 'FontSize', 16);
            ylabel('Predictions (kWh)', 'FontName', 'Times New Roman', 'FontSize', 16);
            x_limits = xlim;
            y_limits = ylim;
            line([x_limits(1), x_limits(2)], [x_limits(1), x_limits(2)], 'Color', 'b', 'LineStyle', '-', 'DisplayName', 'y = x');
            legend('Proposed method', 'Baseline');
            set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
            title('Comparison between simulated energy consumption and measurements.')
            grid on;
            rmsel = sqrt(mean((y - y_predl).^2)); rmsep = sqrt(mean((y - y_predp).^2));
            mapel = mean(abs((y - y_predl) ./ abs(y))) * 100; mapep = mean(abs((y - y_predp) ./ abs(y))) * 100;
            fprintf('Error of energy consumption: \n');
            fprintf('RMSE of baseline: %.2f\n', rmsel);
            fprintf('MAPE of baseline: %.2f\n',mapel);
            fprintf('RMSE of proposed method: %.2f\n', rmsep);
            fprintf('MAPE of proposed method: %.2f\n',mapep);
            % Error box plot
            figure(25);
            boxplot(prop_sef, 'Colors', 'b', 'OutlierSize', 4);
            h = findobj(gca, 'Tag', 'Box');
            for j = 1:length(h)
                patch(get(h(j), 'XData'), get(h(j), 'YData'), 'w', 'FaceAlpha', 0, 'EdgeColor', 'b');
            end
            hMedian = findobj(gca, 'Tag', 'Median');
            for j = 1:length(hMedian)
                hMedian(j).Color = 'r';
            end
            set(findobj(gca, 'type', 'line', 'Tag', 'Box'), 'LineWidth', 3);
            set(findobj(gca, 'type', 'line', 'Tag', 'Median'), 'LineWidth', 3);
            set(findobj(gca, 'type', 'line', 'Tag', 'Whisker'), 'LineWidth', 3);
            set(findobj(gca, 'type', 'line', 'Tag', 'Upper Adjacent Value'), 'LineWidth', 3);
            set(findobj(gca, 'type', 'line', 'Tag', 'Lower Adjacent Value'), 'LineWidth', 3);
            xlabel('Samples');
            ylabel('Error (%)');
            set(gca, 'FontSize', 24);
            title('Error box plot of battery calibration model results.')
            grid on;
            hold off;
        end
    end
end

function features = cal_features(t, v)
    dt = diff(t);
    dv = diff(v);
    acc = [dv ./ dt 0];
    mean_speed = mean(v);
    v_std = std(v);
    mean_vm = mean(v(v>0.1));
    max_vm = max(v(v>0.1));
    mean_accel = mean(acc(acc > 0.01));
    mean_a = mean(acc);
    std_a = std(acc);
    mean_decel = mean(acc(acc < -0.01));
    max_speed = max(v);
    max_accel = max(acc);
    max_decel = min(acc);
    dis = cumtrapz(t,v); dis = dis(end);
    Ti = length(find(v < 0.1 & acc > -0.01 & acc < 0.01))/length(v);
    Td = length(find(v >= 0.1 & acc <= -0.01))/length(v);
    Ta = 1 - Ti - Td;
    features = [max_vm; mean_speed; mean_vm; v_std; mean_a; std_a; mean_accel; max_accel; max_decel; ...
        mean_decel; dis; Ta; Td; Ti];

end
function errs = err_features(meas, pre)
    errs = 100.*(pre-meas)./meas;
end