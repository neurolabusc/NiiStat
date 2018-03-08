function roiNumVec = roiNames2roiNumbers(cellNameArray)
%Takes an array of cells that specify text read in from the ROI file
%and extracts the numbers at the beginning (before the first | character)
%and places them in a 1xX vector of integers (note uses str2double as that
%works on cell array whereas str2int does not automatically work in that
%way.

    roiNumVec = [];
    
    for i = 1:length(cellNameArray)
       roiNumVec = cat(2,roiNumVec,round(str2double(strtok(cellNameArray(i),'|')))); 
        
    end
    
    %roiNumVec

end

