function [dinfo] = fun_data_dirs(whichd)

if whichd == 1
    folders = [113, 112];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Glio\ExVivo_cryo\Pilot01_20220608_CT2A_MB1\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Glio\ExVivo_cryo\Pilot01_20220608_CT2A_MB1\data', num2str(folder)];
    end
elseif whichd == 2
    folders = [14, 13];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Glio\ExVivo_cryo\Pilot02_20220615_GL261_MB1\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Glio\ExVivo_cryo\Pilot02_20220615_GL261_MB1\data', num2str(folder)];
    end
elseif whichd == 3
    folders = [11, 12];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Glio\ExVivo_cryo\Pilot03_20220617_GL261_MB2\20220617_114624_RNH_CRP_CTI_Glio_GL261_MB2_1_2\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Glio\ExVivo_cryo\Pilot03_20220617_GL261_MB2\20220617_114624_RNH_CRP_CTI_Glio_GL261_MB2_1_2\data', num2str(folder)];
    end
elseif whichd == 4
    folders = [9, 10];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Glio\ExVivo_cryo\Pilot04_20220619_CT2A_MB2\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Glio\ExVivo_cryo\Pilot04_20220619_CT2A_MB2\data', num2str(folder)];
    end
elseif whichd == 5
    folders = [14, 15];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Glio\ExVivo_cryo\Pilot05_20220621_GL261_MB3\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Glio\ExVivo_cryo\Pilot05_20220621_GL261_MB3\data', num2str(folder)];
    end
elseif whichd == 6
    folders = [14, 15];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Glio\ExVivo_cryo\Pilot06_20220623_CT2A_MB3\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Glio\ExVivo_cryo\Pilot06_20220623_CT2A_MB3\data', num2str(folder)];
    end
elseif whichd == 7
    folder_name = 'C:\Users\rafae\Data\Glio\ExVivo_cryo\Pilot07_20220707_CT2A_MB1';
    folders = [17, 18];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Glio\ExVivo_cryo\Pilot07_20220707_CT2A_MB1\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Glio\ExVivo_cryo\Pilot07_20220707_CT2A_MB1\data', num2str(folder)];
    end
elseif whichd == 8
    folders = [11, 10];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).name = ['C:\Users\rafae\Data\Glio\ExVivo_cryo\20220711_131834_RNH_CRP_CTI_Glio_MB_control_1_9\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Glio\ExVivo_cryo\20220711_131834_RNH_CRP_CTI_Glio_MB_control_1_9\data', num2str(folder)];
    end
end

