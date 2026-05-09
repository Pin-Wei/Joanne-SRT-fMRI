function [n_all, n_sig] = write_structured_stats(data, csv_path)
% Write structured CSV from data.list2txt / data.list2 / data.list2visible.
% Strips CONN's \NNN\<name>\ markers and parses each row into typed columns.
%
% Columns:
%   type, cluster_id, roi1, roi1_xyz, roi2, roi2_xyz,
%   stat_name, stat_value, df, p_unc, p_FDR, p_FWE, sig
%
% sig=1 when the row index appears in data.list2visible (CONN's idx list of
% "visible" / supra-threshold items at the current preset).
%
% Cluster header formats vary by preset (CONN sprintf padded with %-24s/%-20s,
% so fields are always separated by 2+ spaces — split on \s{2,}):
%   FNC : "Cluster K/N    F(2,18) = 8.41    0.002632    0.036849"     (2 p-vals)
%   SPC : "Cluster K/N    Mass = 12.34      0.001       0.020   0.04" (2 or 3 p-vals)
%   TFCE: "Cluster K/N    TFCE = 5.92       0.000013    0.013183 0.025" (3 p-vals)
% Connection rows always: "Connection roi1 (xyz1)-roi2 (xyz2)  STAT(df) = val  p_unc  p_FDR"

    n_all = 0; n_sig = 0;
    fcsv = fopen(csv_path, 'w');
    fprintf(fcsv, ['type,cluster_id,roi1,roi1_xyz,roi2,roi2_xyz,' ...
                   'stat_name,stat_value,df,p_unc,p_FDR,p_FWE,sig\n']);

    if ~isfield(data,'list2txt') || isempty(data.list2txt)
        fclose(fcsv);
        return;
    end

    raw      = data.list2txt;
    list2    = []; if isfield(data,'list2'),        list2    = data.list2;        end
    list2vis = []; if isfield(data,'list2visible'), list2vis = data.list2visible; end

    % CONN's own marker-strip regex (see conn_displayroi.m line 1741)
    cleaned = regexprep(raw, '\\\\(\d*)\\\\(.*?)\\\\', '$2');

    cluster_id = 0;
    % stat header inside cluster row — MATLAB regexp doesn't emit a token for
    % a capturing group nested inside an unmatched optional (?:...)?, so we
    % split into two patterns instead of one with optional df.
    stathdr_with_df_re = '^(\w+)\(([^)]+)\)\s*=\s*(\S+)\s*$';
    stathdr_no_df_re   = '^(\w+)\s*=\s*(\S+)\s*$';
    conn_re    = ['^Connection\s+(.+?)\s*' ...
                  '\((-?\d+,-?\d+,-?\d+)\)\s*-\s*(.+?)\s*' ...
                  '\((-?\d+,-?\d+,-?\d+)\)\s+' ...
                  '(\w+)\(([^)]+)\)\s*=\s*(\S+)\s+' ...
                  '(\S+)\s+(\S+)$'];

    for i = 1:numel(cleaned)
        s = strtrim(cleaned{i});
        if isempty(s), continue; end

        is_cluster_row = false;
        if size(list2,1) >= i && list2(i,1)==0 && list2(i,2)==0
            is_cluster_row = true;
        elseif startsWith(s, 'Cluster ') || startsWith(s, 'Network ') || startsWith(s, 'ROI ')
            is_cluster_row = true;
        end

        is_sig = 0;
        if ~isempty(list2vis) && any(list2vis == i), is_sig = 1; end
        n_all = n_all + 1;
        if is_sig, n_sig = n_sig + 1; end

        if is_cluster_row
            parts = regexp(s, '\s{2,}', 'split');
            cid_tok = regexp(parts{1}, '^Cluster\s+(\d+)/', 'tokens', 'once');
            if ~isempty(cid_tok), cluster_id = str2double(cid_tok{1}); end

            stat_name = ''; df = ''; stat_value = '';
            if numel(parts) >= 2
                p2 = strtrim(parts{2});
                st = regexp(p2, stathdr_with_df_re, 'tokens', 'once');
                if ~isempty(st)
                    stat_name = st{1}; df = st{2}; stat_value = st{3};
                else
                    st = regexp(p2, stathdr_no_df_re, 'tokens', 'once');
                    if ~isempty(st)
                        stat_name = st{1}; stat_value = st{2};
                    end
                end
            end
            p_unc=''; p_FDR=''; p_FWE='';
            if numel(parts) >= 3, p_unc = strtrim(parts{3}); end
            if numel(parts) >= 4, p_FDR = strtrim(parts{4}); end
            if numel(parts) >= 5, p_FWE = strtrim(parts{5}); end

            fprintf(fcsv, 'cluster,%d,,,,,%s,%s,"%s",%s,%s,%s,%d\n', ...
                cluster_id, stat_name, stat_value, df, p_unc, p_FDR, p_FWE, is_sig);
        else
            tok = regexp(s, conn_re, 'tokens', 'once');
            if ~isempty(tok)
                roi1 = strtrim(tok{1}); roi2 = strtrim(tok{3});
                fprintf(fcsv, 'connection,%d,"%s","%s","%s","%s",%s,%s,"%s",%s,%s,,%d\n', ...
                    cluster_id, roi1, tok{2}, roi2, tok{4}, ...
                    tok{5}, tok{7}, tok{6}, tok{8}, tok{9}, is_sig);
            else
                fprintf(fcsv, 'connection,%d,,,,,,,,,,,%d\n', cluster_id, is_sig);
            end
        end
    end
    fclose(fcsv);
end
