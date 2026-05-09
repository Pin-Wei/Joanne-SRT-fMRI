addpath('/home/aclexp/Software/conn');
S = load('/home/aclexp/pinwei/Joanne_SRT_fMRI/conn_out/no_PM_260421/results/bookmarks/2026_04_23_094032219.bookmark.mat');
disp('--- class ---'); disp(class(S.conn_args));
disp('--- size ---'); disp(size(S.conn_args));
if iscell(S.conn_args)
    for k = 1:numel(S.conn_args)
        fprintf('[%d] class=%s size=%s\n', k, class(S.conn_args{k}), mat2str(size(S.conn_args{k})));
        if ischar(S.conn_args{k}), disp(['    value: ' S.conn_args{k}]); end
        if isnumeric(S.conn_args{k}) && numel(S.conn_args{k})<40, disp(S.conn_args{k}); end
        if isstruct(S.conn_args{k})
            disp('    fields:'); disp(fieldnames(S.conn_args{k}));
            try, disp(S.conn_args{k}); end
        end
        if iscell(S.conn_args{k})
            fprintf('    cell %d elems\n', numel(S.conn_args{k}));
            for j = 1:min(numel(S.conn_args{k}),15)
                v = S.conn_args{k}{j};
                if ischar(v), fprintf('      {%d} char: %s\n', j, v);
                elseif isnumeric(v), fprintf('      {%d} num: %s\n', j, mat2str(v));
                else, fprintf('      {%d} class=%s\n', j, class(v)); end
            end
        end
    end
end
