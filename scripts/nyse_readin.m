clc % clears the command window
close all %removes all figures
clear
fileName='ahp.csv';
fullfileName=fullfile('../data/nyseold', fileName);
files=dir('../data/nyseold/*.csv');
tt = readtable(fullfileName);
tt.Var1=string(tt.Var1);
tt.Var1=strcat("19", tt.Var1);
tt.Var1=str2double(tt.Var1);
tt.Var1=datetime(tt.Var1, 'ConvertFrom','yyyymmdd');
tt=renamevars(tt, "Var2" ,fileName(1:end-4));
tt=table2timetable(tt);
tt=retime(tt, 'daily', 'fillwithconstant', 'Constant', 1);
for i=1:length(files)
    if files(i).name=="ahp.csv"
        continue
    elseif files(i).name=="NyseTicker.csv"
        continue
    elseif files(i).name=="NyseTickerMerged.csv"
        continue
    else
        fileName=fullfile('../data/nyseold', files(i).name);
        newtt = readtable(fileName);
        newtt.Var1=string(newtt.Var1);
        newtt.Var1=strcat("19", newtt.Var1);
        newtt.Var1=str2double(newtt.Var1);
        newtt.Var1=datetime(newtt.Var1, 'ConvertFrom','yyyymmdd');
        newtt=renamevars(newtt, "Var2" ,files(i).name(1:end-4));
        newtt=table2timetable(newtt);
        if files(i).name=="tex.csv"
            mask=newtt.Var1=="21-Jan-0000";
            newtt=newtt(~mask,:);
        end
        newtt=retime(newtt, 'daily', 'fillwithconstant', 'Constant', 1);
        tt=synchronize(tt, newtt);
    end
end
nyseTimeTable=tt;
save("nyse.mat", "nyseTimeTable");