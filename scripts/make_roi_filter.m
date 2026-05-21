function roiFilter = make_roi_filter(roiNames, ROI_PREFIX)
    if matches(ROI_PREFIX, 'joanne')
        roiSetFS = { ...
            'precentral', ...               % SL
            'Caudate', ...                  % SL
            'Putamen', ...                  % SL
            'Cerebellum', ...               % SL
            'inferiorparietal', ...         % IF
            'superiorparietal', ...         % IF
            'insula', ...                   % IF
            'caudalanteriorcingulate', ...  % IF
            'rostralanteriorcingulate', ... % IF
            'Amygdala', ...                 % IF
            'Thalamus' ...                  % IF
        };
        roiSetBA = { ...
            'BA.ctx-rh-BA6',  'BA.ctx-lh-BA6', ... % supplementary motor area (SMA); SL system
            'BA.ctx-rh-BA4',  'BA.ctx-lh-BA4', ... % primary motor cortex (M1); SL system
            'BA.ctx-rh-BA9',  'BA.ctx-lh-BA9', ... % dorsolateral prefrontal cortex (dlPFC); IF system
            'BA.ctx-rh-BA46', 'BA.ctx-lh-BA46' ... % also dlPFC
        };
        roiFilter = cellfun(@(s) ...
            startsWith(s, 'FS') && any(contains(s, roiSetFS)) || ...
            any(matches(s, roiSetBA)), ...
            roiNames ...
        );
    else
        roiFilter = cellfun(@(s) startsWith(s, ROI_PREFIX), roiNames);
    end
end