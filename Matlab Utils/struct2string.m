function str = struct2string(strct)

fnames = fieldnames(strct);
vals = struct2cell(strct);

% Sort case-insensitive.
[~, I] = sort(lower(fnames));
fnames = fnames(I);
vals = vals(I);



% Convert sub-structs to string first.
for(i = 1:length(vals))
    if(isstruct(vals{i}))
        vals{i} = struct2string(vals{i});
    end
    vals{i} = [num2str(vals{i}, '%.2g')];
end

dat = [fnames, vals].';
str = [dat{:}];