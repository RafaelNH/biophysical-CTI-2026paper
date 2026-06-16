close all
clear all
clc

fs = filesep;
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'MC_utils/'])
addpath(['..', fs, 'DWIu_toolbox_v1.11', fs, 'DKI_toolbox/'])

f = 0.7;
for sub = 1:2
    ai = 1;
    for seed = 123
        for As = 0:0.1:0.6
            as = As * 100;
            extraD = 150;
            intraD = 250;
            
            substrate = {'intra', 'extra'};
            
            % load data
            pth = ['MC_results_from_step001', fs];
            
            if sub == 1
                pos0 = load([pth,'pos_',substrate{sub},'_as',num2str(as),'_seed',num2str(seed),'_intraD',num2str(intraD),'_0.txt']);
                pos1 = load([pth,'pos_',substrate{sub},'_as',num2str(as),'_seed',num2str(seed),'_intraD',num2str(intraD),'_1.txt']);
            else
                pos0 = load([pth,'pos_',substrate{sub},'_as',num2str(as),'_seed',num2str(seed),'_f',num2str(f*100),'_extraD',num2str(extraD),'_0.txt']);
                pos1 = load([pth,'pos_',substrate{sub},'_as',num2str(as),'_seed',num2str(seed),'_f',num2str(f*100),'_extraD',num2str(extraD),'_1.txt']);
            end
            
            % convert R from m to um
            pos1 = pos1 * 1000000;
            pos0 = pos0 * 1000000;
            
            figure('color', [1 1 1])
            subplot(2, 2, 1)
            plot3(pos0(:, 1), pos0(:, 2), pos0(:, 3), '.')
            view(0, 90)
            
            R1 = pos1 - pos0;
            
            time = 10; %ms
            
            d11 = fun_reconst_d_from_pos(R1, time, 1, 1);
            d22 = fun_reconst_d_from_pos(R1, time, 2, 2);
            d33 = fun_reconst_d_from_pos(R1, time, 3, 3);
            d12 = fun_reconst_d_from_pos(R1, time, 1, 2);
            d13 = fun_reconst_d_from_pos(R1, time, 1, 3);
            d23 = fun_reconst_d_from_pos(R1, time, 2, 3);
            dt = [d11, d22, d33, d12, d13, d23];
            
            md = (d11 + d22 + d33)/3;
            
            [MD, AD, RD, FA, V1, V2, V3, L2, L3] = dti_metrics(dt, 1);
            AD_all(sub, ai) = AD;
            RD_all(sub, ai) = RD;
            
%             subplot(2, 2, 2)
%             plot3(pos0(:, 1), pos0(:, 2), pos0(:, 3), '.')
%             view(0, 90)
%             
             subplot(2, 2, 3)
             plot3(pos1(:, 1), pos1(:, 2), pos1(:, 3), '.')
             view(0, 90)
%             
%             subplot(2, 2, 4)
%             [y, x] = hist(R1(:, 1), -30:0.5:30);
%             plot(x, y)
%             hold on
%             [y, x] = hist(R1(:, 2), -30:0.5:30);
%             plot(x, y)
%             [y, x] = hist(R1(:, 3), -30:0.5:30);
%             plot(x, y)
            
            %plot3(pos1(:, 1)*1000000, pos1(:, 2)*1000000, pos1(:, 3)*1000000, '.')
            k1111 = fun_reconst_k_from_pos(R1, 1, 1, 1, 1);
            k2222 = fun_reconst_k_from_pos(R1, 2, 2, 2, 2);
            k3333 = fun_reconst_k_from_pos(R1, 3, 3, 3, 3);
            k1112 = fun_reconst_k_from_pos(R1, 1, 1, 1, 2);
            k1113 = fun_reconst_k_from_pos(R1, 1, 1, 1, 3);
            k1222 = fun_reconst_k_from_pos(R1, 1, 2, 2, 2);
            k2223 = fun_reconst_k_from_pos(R1, 2, 2, 2, 3);
            k1333 = fun_reconst_k_from_pos(R1, 1, 3, 3, 3);
            k2333 = fun_reconst_k_from_pos(R1, 2, 3, 3, 3);
            k1122 = fun_reconst_k_from_pos(R1, 1, 1, 2, 2);
            k1133 = fun_reconst_k_from_pos(R1, 1, 1, 3, 3);
            k2233 = fun_reconst_k_from_pos(R1, 2, 2, 3, 3);
            k1123 = fun_reconst_k_from_pos(R1, 1, 1, 2, 3);
            k1223 = fun_reconst_k_from_pos(R1, 1, 2, 2, 3);
            k1233 = fun_reconst_k_from_pos(R1, 1, 2, 3, 3);
            
            k = [k1111, k2222, k3333, k1112, k1113, ...
                k1222, k2223, k1333, k2333, k1122,...
                k1133, k2233, k1123, k1223, k1233];
            
            ss2 = mean(dot(R1', R1'))^2;
            % sqrt(ss2)/(2*3*time) is the mean diffusivity
            wt = 9*k/ss2;
            
            [MK, AK, RK] = dki_tensor_metrics(dt, wt, 1);
            
            AK_all(sub, ai) = AK;
            RK_all(sub, ai) = RK;
            MK_all(sub, ai) = MK;
            ai = ai + 1;
            
            %figure('color', [1 1 1])
            % 3D geometries
            np=100;
            theta=linspace(0,2*pi,np+1);
            phi=linspace(0,2*pi,np+1);
            [theta, phi]=meshgrid(theta, phi);
            
            scal = (length(parula):-1:1)/length(parula);
            mycmap = [scal', zeros(length(parula), 2); 0, 0, 0; parula];
            
            
            subplot(1, 4, 1+2*(sub-1))
            rd = DirectionalDiff_2D(dt, theta, phi);
            [xx, yy, zz]=sph2cart(theta, phi, rd);
            surf(xx, yy, zz, rd, 'EdgeColor', 'none')
            title('Dapp'), xlabel('x'), ylabel('y'), zlabel('z')
            view(20,20)
            caxis([-0.6 0.6])
            xlim([-2 2])
            ylim([-2 2])
            zlim([-2 2])
            axis square
            grid off
            colormap(mycmap)
            colorbar
            
            subplot(1, 4, 2+2*(sub-1))
            r = DirectionalKurt_2D(dt, wt, theta, phi);
            [xx, yy, zz]=sph2cart(theta, phi, r);
            plot3([0 0], [0, 0], [-1 1], 'linewidth', 2, 'color', 'black')
            hold on
            surf(xx, yy, zz, r, 'EdgeColor', 'none')
            title('KT'), xlabel('x'), ylabel('y'), zlabel('z')
            view(-40,40)
            caxis([-0.6 0.6])
            xlim([-0.8 0.8])
            ylim([-0.8 0.8])
            zlim([-0.8 0.8])
            axis square
            grid off
            colormap(mycmap)
            
            
            if sub == 1
                dt_in = dt;
                wt_in = wt;
                save([pth,'tensors_as',num2str(as),'_intraD',num2str(intraD),'_seed',num2str(seed)], 'wt_in', 'dt_in');
            else
                dt_ex = dt;
                wt_ex = wt;
                save([pth,'tensors_as',num2str(as),'_f',num2str(f*100),'_extraD',num2str(extraD),'_seed',num2str(seed)], 'wt_ex', 'dt_ex');
            end
        end
        close all
    end
end
close all

figure
subplot(2, 2, 1)
plot(0:0.1:0.6, AD_all(1, 1:7),'r')
hold on
plot(0:0.1:0.6, RD_all(1, 1:7),'b')
legend('AD', 'RD')
title('intra')

subplot(2, 2, 2)
plot(0:0.1:0.6, AK_all(1, 1:7),'r')
hold on
plot(0:0.1:0.6, RK_all(1, 1:7),'b')
legend('AK', 'RK')
title('intra')

subplot(2, 2, 3)
plot(0:0.1:0.6, AD_all(2, 1:7),'r')
hold on
plot(0:0.1:0.6, RD_all(2, 1:7),'b')
legend('AD', 'RD')
title('extra')

subplot(2, 2, 4)
plot(0:0.1:0.6, AK_all(2, 1:7),'r')
hold on
plot(0:0.1:0.6, RK_all(2, 1:7),'b')
legend('AK', 'RK')
title('extra')
ylim([-0.05 0.2])
