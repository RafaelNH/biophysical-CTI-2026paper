close all
clear all
clc

% add paths
fs = filesep;

denoised = true;
complex = true;
if denoised
    den = '_den';
else
    den = '';
end
if complex
    c = '_complex';
else
    c = '';
end

% load data
for whichd = 1
    dinfo = fun_data_dirs(whichd);
    
    if dinfo(1).folder>9
        load([dinfo(1).savename(1:end-6), 'MASKS.mat'], 'mask_final')
    else
        load([dinfo(1).savename(1:end-5), 'MASKS.mat'], 'mask_final')
    end
    
    for di = 1:length(dinfo)
        file_name = dinfo(di).savename;
        pt=[file_name, c, den, '_align'];
        load(pt)
        
        bt = gtab.bval1 + gtab.bval2;
        
        lowb = min(bt);
        
        % create some variable for QA assessment
        siz = size(data);
        x = 1:siz(4); % experiements indexes
        xb0 = x(bt==lowb); % low b-valud indexs
        xdi = x(~(bt==lowb));
        b0s = data(:, :, :, bt==lowb);
        dis = data(:, :, :, ~(bt==lowb));
        sizb0 = size(b0s);
        sizdi = size(dis);
        
        % First QA - observe if mouse bed moved during acquisitions by plotting
        % first and last lowb acquisitions and their difference
        first_b0 = squeeze(data(end:-1:1, :, :, xb0(1)));
        last_b0 = squeeze(data(end:-1:1, :, :, xb0(end)));
        diff_b0 = last_b0 - first_b0;
        
        figure('color', [1 1 1])
        subplot(3, 1, 1)
        imagesc(first_b0(:, :))
        subplot(3, 1, 2)
        imagesc(last_b0(:, :))
        subplot(3, 1, 3)
        imagesc(diff_b0(:, :))
        
        % Second QA - plot mean b0 values for each slice
        mean_int = zeros(sizb0(3:4));
        mean_di = zeros(sizdi(3:4));
        mean_all = zeros(siz(3:4));
        
        for si = 1:sizb0(3)
            mask_si = mask_final(:, :, si);
            for vi = 1:siz(4)
                dii = squeeze(data(:, :, si, vi));
                mean_all(si, vi) = mean(dii(mask_si(:)==1));
            end
            mean_int(si, :) = mean_all(si, bt==lowb);
            mean_di(si, :)  = mean_all(si, ~(bt==lowb));
        end
        
        figure('color', [1 1 1])
        subplot(3, 1, 1)
        plot(xb0, mean_int')
        hold on
        xline(159, 'red')
        xline(317, 'red')
        xline(475, 'red')
        xlim([0, siz(4)])
        legend('s=1', 's=2', 's=3', 's=4', 's=5', 's=6')
        title('b0s')
        
        subplot(3, 1, 2)
        plot(xdi, mean_di')
        xline(159, 'red')
        xline(317, 'red')
        xline(475, 'red')
        xlim([0, siz(4)])
        title('DWIs')
        
        subplot(3, 1, 3)
        plot(mean_int')
        xlim([0, sizb0(4)])
        
        % Code to inspect issues
        %     db0s = b0s;
        %     meanb0 = mean(b0s, 4);
        %     for vi = 1:sizb0(4)
        %         db0s(:, :, :, vi) = b0s(:, :, :, vi) - meanb0;
        %     end
        %
        
        %     figure('color', [1 1 1])
        %     subplot(2, 1, 2)
        %     b0_sel = db0s(:, :, 4, 1:9:92);
        %     imagesc(b0_sel(:, :))
        %
        %     subplot(2, 1, 1)
        %     b0_sel = b0s(:, :, 4, 1:9:92);
        %     imagesc(b0_sel(:, :))
        %
        % figure('color', [1 1 1])
        % subplot(2, 1, 2)
        % b0_sel = db0s(:, :, 1:6, 11);
        % imagesc(b0_sel(:, :))
        % %
        % subplot(2, 1, 1)
        % b0_sel = b0s(:, :, 1:6, 11);
        % imagesc(b0_sel(:, :))
        
        % Detection of outliers b0s
        b0s_out = false(sizb0(3:4));
        
        b0_trend = mean(mean_int, 1);
        mean_slice = mean(mean_int, 2);
        mean_diff = mean_slice - mean(mean_slice);
        
        % corrected profiles
        mean_int_c = mean_int;
        
        figure('color', [1 1 1])
        for si = 1:siz(3)
            pred =  mean_diff(si) + b0_trend;
            
            subplot(siz(3), 1, si)
            plot(xb0, mean_int(si, :))
            hold on
            plot(xb0, pred)
            subplot(siz(3), 1, si)
            
            res = mean_int(si, :) - pred;
            stddif = std(res);
            %plot(xb0, abs(res/stddif))
            
            b0s_out(si, :) = abs(res/stddif) > 4;
            mean_int_c(si, b0s_out(si, :)) = pred(b0s_out(si, :));
        end
        
        figure('color', [1 1 1])
        for si = 1:siz(3)
            plot(xb0, mean_int(si, :))
            hold on
            plot(xb0(b0s_out(si, :)), mean_int(si, b0s_out(si, :)), 'o')
            disp(xb0(b0s_out(si, :)))
        end
        xline(159, 'red')
        xline(317, 'red')
        xline(475, 'red')
        xlim([0, siz(4)])
        title('b0s')
        
        % Detection of outliers diff
        data_out = false(sizdi(3:4));
        
        data_trend = mean(mean_di, 1);
        slice_mean = mean(mean_di, 2);
        slice_diff = slice_mean - mean(slice_mean);
        
        figure('color', [1 1 1])
        for si = 1:siz(3)
            subplot(siz(3), 1, si)
            plot(xdi, mean_di(si, :))
            hold on
            plot(xdi, slice_diff(si) + data_trend)
            subplot(siz(3), 1, si)
            res = mean_di(si, :) - slice_diff(si) - data_trend;
            stddif = std(res);
            plot(xdi, abs(res/stddif))
            
            data_out(si, :) = abs(res/stddif) > 3;
        end
        
        figure('color', [1 1 1])
        for si = 1:siz(3)
            plot(xdi, mean_di(si, :))
            hold on
            plot(xdi(data_out(si, :)), mean_di(si, data_out(si, :)), 'o')
            disp(xdi(data_out(si, :)))
        end
        xline(159, 'red')
        xline(317, 'red')
        xline(475, 'red')
        xlim([0, siz(4)])
        title('all')
        
        %n=4;
        % Signal drift estimation
        %figure('color', [1 1 1])
        % for si = 1:siz(3)
        %     xb0_sel = xb0(~b0s_out(si, :));
        %     b0_sel = mean_int(si, ~b0s_out(si, :));
        %     p = polyfit(xb0_sel, b0_sel, n);
        %     s0_pred = polyval(p, x);
        %
        %     plot(xb0, mean_int(si, :))
        %     hold on
        %     plot(xb0(b0s_out(si, :)), mean_int(si, b0s_out(si, :)), 'o')
        %     plot(x, s0_pred)
        % end
        % xline(159, 'red')
        % xline(317, 'red')
        % xline(475, 'red')
        % xlim([0, siz(4)])
        % title('b0s')
        s0_raw = 0;
        nv = sum(mask_final(:));
        for si=1:siz(3)
            mask_si = mask_final(:, :, si);
            w = sum(mask_si(:)==1) / nv;
            s0_raw = s0_raw + mean_int_c(si, :) * w;
            disp(w)
        end
        
        figure('color', [1 1 1])
        plot(xb0, s0_raw, 'color', 'red', 'linewidth', 5)
        hold on
        
        %n=2;
        %p = polyfit(xb0, s0_raw, n);
        movs0 = movmean(s0_raw, 11);
        %s0_pred = polyval(p, x);
        s0_pred = interp1(xb0, movs0, x, 'pchip');
        plot(x, s0_pred)
        
        s0_final = s0_pred;
        %s0_final = spline(xb0, s0_raw, x);
        plot(x, s0_final, 'color', 'black', 'linewidth', 5)
        %plot(x, s0_pred, 'color', 'blue', 'linewidth', 5)
        %ylim([17 19])
        
        data_c = zeros(size(data));
        for vi=1:siz(4)
            data_sel = data(:, :, :, vi);
            data_sel(mask_final(:)==1) = data_sel(mask_final(:)==1)/s0_final(vi);
            data_sel(mask_final(:)==0) = 0;
            data_c(:, :, :, vi) = data_sel;
        end
        
        data = data_c;
        
        % Last QA - repeat plot mean b0 values for each slice
        mean_int = zeros(sizb0(3:4));
        mean_di = zeros(sizdi(3:4));
        mean_all = zeros(siz(3:4));
        
        for si = 1:sizb0(3)
            mask_si = mask_final(:, :, si);
            for vi = 1:siz(4)
                dii = squeeze(data(:, :, si, vi));
                mean_all(si, vi) = mean(dii(mask_si(:)==1));
            end
            mean_int(si, :) = mean_all(si, bt==lowb);
            mean_di(si, :)  = mean_all(si, ~(bt==lowb));
        end
        
        figure('color', [1 1 1])
        subplot(3, 1, 1)
        plot(xb0, mean_int')
        hold on
        xline(159, 'red')
        xline(317, 'red')
        xline(475, 'red')
        xlim([0, siz(4)])
        legend('s=1', 's=2', 's=3', 's=4', 's=5', 's=6')
        title('b0s')
        
        subplot(3, 1, 2)
        plot(xdi, mean_di')
        xline(159, 'red')
        xline(317, 'red')
        xline(475, 'red')
        xlim([0, siz(4)])
        title('DWIs')
        
        save([file_name, c, den, '_align_drift'], 'data', 'gtab')
    end
end