clc % clears the command window
close all %removes all figures
clear
fileName='ahp.csv';
fullfileName=fullfile('../data/nysemerged', fileName);
files=dir('../data/nysemerged/*.csv');
tt = readtable(fullfileName);
tt.Var1=string(tt.Var1);
tt.Var1= cellfun(@(x)sprintf('%06s',x),tt.Var1,'uni',false);
tt.Var1=convertCharsToStrings(tt.Var1);
tt.Var1(1:9443)=strcat("19", tt.Var1(1:9443));
tt.Var1(9444:end)=strcat("20", tt.Var1(9444:end));
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
    elseif files(i).name=="Bond.csv"
        continue
    elseif files(i).name=="Cash.csv"
        continue
    else
        fileName=fullfile('../data/nysemerged', files(i).name);
        newtt = readtable(fileName);
        newtt.Var1=string(newtt.Var1);
        newtt.Var1= cellfun(@(x)sprintf('%06s',x),newtt.Var1,'uni',false);
        newtt.Var1=convertCharsToStrings(newtt.Var1);
        newtt.Var1(1:9443)=strcat("19", newtt.Var1(1:9443));
        newtt.Var1(9444:end)=strcat("20", newtt.Var1(9444:end));
        newtt.Var1=str2double(newtt.Var1);
        newtt.Var1=datetime(newtt.Var1, 'ConvertFrom','yyyymmdd');
        newtt=renamevars(newtt, "Var2" ,files(i).name(1:end-4));
        newtt=table2timetable(newtt);
        newtt=retime(newtt, 'daily', 'fillwithconstant', 'Constant', 1);
        tt=synchronize(tt, newtt);
    end
end
nyseMergedTimeTable=tt;
save("nysemerged.mat", "nyseMergedTimeTable");
