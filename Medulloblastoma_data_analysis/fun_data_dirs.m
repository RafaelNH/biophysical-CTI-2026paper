function [dinfo] = fun_data_dirs(whichd)

if whichd==0.001 % one ave
    folders = 12;
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).folder = folder;
        dinfo(fi).name = ['C:\Users\rafae\Data\Cerebrum\Tumour1\Pilot_acquisitions\room_temp\12\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Cerebrum\Tumour1\Pilot_acquisitions\room_temp\12\data', num2str(folder)];
    end
elseif whichd == 0.01 % 3 ave
    folders = 18;
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).folder = folder;
        dinfo(fi).name = ['C:\Users\rafae\Data\Cerebrum\Tumour1\Pilot_acquisitions\room_temp\18\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Cerebrum\Tumour1\Pilot_acquisitions\room_temp\18\data', num2str(folder)];
    end
elseif whichd == 0.02 % 3 ave
    folders = 19;
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).folder = folder;
        dinfo(fi).name = ['C:\Users\rafae\Data\Cerebrum\Tumour1\Pilot_acquisitions\room_temp\18\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Cerebrum\Tumour1\Pilot_acquisitions\room_temp\18\data', num2str(folder)];
    end
elseif whichd == 0.1 % 3 ave 2 rep
    folders = [18, 19];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).folder = folder;
        dinfo(fi).name = ['C:\Users\rafae\Data\Cerebrum\Tumour1\Pilot_acquisitions\room_temp\18\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Cerebrum\Tumour1\Pilot_acquisitions\room_temp\18\data', num2str(folder)];
    end
elseif whichd == 1 % tumour 1
    folders = 25;
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).folder = folder;
        dinfo(fi).name = ['C:\Users\rafae\Data\Cerebrum\Tumour1\25\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Cerebrum\Tumour1\25\data', num2str(folder)];
    end
elseif whichd == 2 % Control 1
    folders = 38;
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).folder = folder;
        dinfo(fi).name = ['C:\Users\rafae\Data\Cerebrum\Control1\38\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Cerebrum\Control1\38\data', num2str(folder)];
    end
elseif whichd == 3 % tumour 2
    folders = 47;
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).folder = folder;
        dinfo(fi).name = ['C:\Users\rafae\Data\Cerebrum\Tumour2\47\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Cerebrum\Tumour2\47\data', num2str(folder)];
    end
elseif whichd == 4 % tumour 3
    folders = 77;
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).folder = folder;
        dinfo(fi).name = ['C:\Users\rafae\Data\Cerebrum\Tumour3\77\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Cerebrum\Tumour3\77\data', num2str(folder)];
    end
elseif whichd == 5 % control 2
    folders = 89;
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).folder = folder;
        dinfo(fi).name = ['C:\Users\rafae\Data\Cerebrum\Control2\89\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Cerebrum\Control2\89\data', num2str(folder)];
    end
elseif whichd == 6 % control 3
    folders = 103;
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).folder = folder;
        dinfo(fi).name = ['C:\Users\rafae\Data\Cerebrum\Control3\103\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Cerebrum\Control3\103\data', num2str(folder)];
    end
elseif whichd == 7 % tumour 4
    folders = 113;
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).folder = folder;
        dinfo(fi).name = ['C:\Users\rafae\Data\Cerebrum\Tumour4\113\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Cerebrum\Tumour4\113\data', num2str(folder)];
    end
elseif whichd == 7.1 % tumour 4 (180 test)
    folders = 114;
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).folder = folder;
        dinfo(fi).name = ['C:\Users\rafae\Data\Cerebrum\Tumour4\113\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Cerebrum\Tumour4\113\data', num2str(folder)];
    end
elseif whichd == 7.2 % tumour 4 (90 and 180 test)
    folders = [113, 114];
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).folder = folder;
        dinfo(fi).name = ['C:\Users\rafae\Data\Cerebrum\Tumour4\113\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Cerebrum\Tumour4\113\data', num2str(folder)];
    end
elseif whichd == 8 % Control 4 
    folders = 126; % 126 test, 127 retest, 128 data with 180 
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).folder = folder;
        dinfo(fi).name = ['C:\Users\rafae\Data\Cerebrum\Control4\126\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Cerebrum\Control4\126\data', num2str(folder)];
    end
elseif whichd == 9
    % Control 5
    folders = 137;  
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).folder = folder;
        dinfo(fi).name = ['C:\Users\rafae\Data\Cerebrum\Control5\137\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Cerebrum\Control5\137\data', num2str(folder)];
    end
elseif whichd == 10
    % Tumour 5
    folders = 149;  % 149 data 90, 146 data 180
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).folder = folder;
        dinfo(fi).name = ['C:\Users\rafae\Data\Cerebrum\Tumour5\149\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Cerebrum\Tumour5\149\data', num2str(folder)];
    end
elseif whichd == 11
    % Control 6
    folders = 160;  % 126 test, 127 retest, 128 data with 180 
    nf = length(folders);
    for fi=1:nf
        folder = folders(fi);
        dinfo(fi).folder = folder;
        dinfo(fi).name = ['C:\Users\rafae\Data\Cerebrum\Control6\160\', num2str(folder)];
        dinfo(fi).savename = ['C:\Users\rafae\Data\Cerebrum\Control6\160\data', num2str(folder)];
    end
end
