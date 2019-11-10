function sout = nii_mergestruct(new,old)
%merges structures new and old
% if same fields exist in both: values from new overwrite old
%Chris Rorden - inspired by http://stackoverflow.com/questions/38645/what-are-some-efficient-ways-to-combine-two-structures-in-matlab
sout = new;
for f = fieldnames(old)'
    cf = char(f);
    if (not(isfield(sout,cf)))
        %this is unique
        sout = setfield(sout,cf,getfield(old,cf));
    else
        %both old and new have this field
        for f2 = fieldnames(old.(cf))'
            cf2 = char(f2);
            if (not(isfield(sout.(cf),cf2)))
                sout.(cf) = setfield(sout.(cf),cf2,getfield(old.(cf),cf2));
            end
        end
    end
end

% sout = new;
% for f = fieldnames(old)'
%     cf = char(f);
%     if (not(isfield(sout,cf)))
%         sout = setfield(sout,cf,getfield(old,cf));
%     end
% end
%end nii_mergestruct()