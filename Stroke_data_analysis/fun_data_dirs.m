function [dinfo] = fun_data_dirs(whichd)

if whichd == 1
    folders = [25, 28:33];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke1\20200626_131058_RA_stroke_exvivo1_RA_stroke_exvivo_trial2_1_2\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke1\data', num2str(folder)];
    end
elseif whichd == 2
    folders = 36:42;
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke2\20200628_114309_RA_stroke_exvivo2_RA_stroke_exvivo_trial2_1_2\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke2\data', num2str(folder)];
    end
elseif whichd == 3
    folders = [16:21,23];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke3\20200630_101648_RA_stroke_exvivo3_RA_stroke_exvivo_trial2_1_2\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke3\data', num2str(folder)];
    end
elseif whichd == 4
    folders = 57:63;
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke4\20200702_092644_RA_stroke_exvivo4_RA_stroke_exvivo_trial2_1_2\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke4\data', num2str(folder)];
    end
elseif whichd == 5
    folders = 45:51;
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke5\20200704_103824_RA_stroke_exvivo5_RA_control_exvivo_trial2_1_2\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke5\data', num2str(folder)];
    end
elseif whichd == 6
    folders = [32:38,46]; %39 has 180 experiements
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke6\20210620_125510_RA_stroke1_n_RA_stroke1_n_1_1\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke6\data', num2str(folder)];
    end
elseif whichd == 7
    folders = [125, 128:133];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke7\20210624_105052_RA_stroke2v2_n_RA_stroke2v2_n_1_1\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke7\data', num2str(folder)];
    end
elseif whichd == 8
    folders = [41, 45:50]; %51-52 has 180 experiements
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke8\20210622_124936_RA_stroke3_n_RA_stroke3_n_1_1\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke8\data', num2str(folder)];
    end
elseif whichd == 9
    folders = [27:33];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke9\20210706_143727_RA_stroke4_n_RA_stroke4_n_1_1\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke9\data', num2str(folder)];
    end
elseif whichd == 10
    folders = [76,95:98, 100, 103]; %101, 102 has 180 experiements (99 corrupted replaced by 103)
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke10\20210611_192535_RA_stroke1_batch2_test_RA_stroke1_batch2_te_1_1\20210611_192535_RA_stroke1_batch2_test_RA_stroke1_batch2_te_1_1\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_exvivo\Stroke10\data', num2str(folder)];
    end
elseif whichd == 1001
    folders = 34:40; % extra 44
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_exvivo\Control1\20200627_125324_RA_control_exvivo1_RA_control_exvivo_trial2_1_2\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_exvivo\Control1\data', num2str(folder)];
    end
elseif whichd == 1002
    folders = [26:32];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_exvivo\Control2\20200629_102210_RA_control_exvivo2_RA_control_exvivo_trial2_1_2\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_exvivo\Control2\data', num2str(folder)];
    end
elseif whichd == 1003
    folders = 19:25;
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_exvivo\Control3\20200701_115302_RA_control_exvivo3_RA_control_exvivo_trial2_1_2\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_exvivo\Control3\data', num2str(folder)];
    end
elseif whichd == 1004
    folders = [46:50, 53]; % 52 was removed because it seems corrupted by ghosts
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_exvivo\Control4\20200703_101252_RA_control_exvivo4_RA_control_exvivo_trial2_1_2\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_exvivo\Control4\data', num2str(folder)];
    end
elseif whichd == 1005
    folders = [52, 37:38, 40, 49:51]; % 42 DtiEpi
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_exvivo\Control5\20200707_094826_RA_control_exvivo5_RA_control_exvivo_trial2_1_2\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_exvivo\Control5\data', num2str(folder)];
    end
elseif whichd == 21
    folders = [26];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_invivo_longitudinal\Stroke1\3d\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_invivo_longitudinal\Stroke1\3d\data', num2str(folder)];
    end
elseif whichd == 22
    folders = [11];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_invivo_longitudinal\Stroke2\3d\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_invivo_longitudinal\Stroke2\3d\data', num2str(folder)]';
    end
elseif whichd == 23
    folders = [16];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_invivo_longitudinal\Stroke3\3d\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_invivo_longitudinal\Stroke3\3d\data', num2str(folder)]';
    end
elseif whichd == 24
    folders = [13];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_invivo_longitudinal\Stroke4\3d\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_invivo_longitudinal\Stroke4\3d\data', num2str(folder)]';
    end
elseif whichd == 25
    folders = [17];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_invivo_longitudinal\Stroke5\3d\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_invivo_longitudinal\Stroke5\3d\data', num2str(folder)]';
    end
elseif whichd == 26
    folders = [22];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Stroke_invivo_longitudinal\Stroke6\3d\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Stroke_invivo_longitudinal\Stroke6\3d\data', num2str(folder)]';
    end
end

