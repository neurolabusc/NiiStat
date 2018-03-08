function outty = GetROINamesFromFile(fname,indx)

          fileID = fopen(fname,'r');
          C = textscan(fileID,'%s %s %s %s %s %s','Delimiter','|');
          
          fclose(fileID);  
          
          for i = 1:length(C{1})
              C{4}{i} = strcat(   C{1}{i} ,'|',C{2}{i},'|',C{3}{i} );
          end
          
          outty = C{4};
                 
end

% 
%  fileID = fopen('jhu.txt');
%  while fgetl(fileID) ~= -1
%     ln = fgetl(fileID);
%     fprintf(ln);
%     %fprintf('\n');
%     
%  end
%     