function saveCellArrayToExcel(array,filename,sheetNumber)

T = cell2table(array(2:end,:),'variablenames',array(1,:));
writetable(T,filename,'sheet',sheetNumber);

end