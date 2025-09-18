
time0=now;
forceflag = false;
%ADNI_load_data;
disp((now-time0)*24*3600)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Perform cross-validated prediction
% Generate cvind keeping observations for given subject together
% Include all timepoints i.e., DXvec==1 for MCI, or DXvec==0 for CN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is now included in the loop
% k = 10;
% cvind_subj = crossvalind('Kfold',length(RIDlist),k); % key is that these are unique RIDs
% cvind = NaN(size(RIDvec_adnimerge));

% Ensures observations for given subject together
% for ki = 1:k
%   cvind(ismember(RIDvec_adnimerge,RIDlist(cvind_subj==ki))) = ki;
% end

% Progressing to DEM vs Not Demented (CN/MCI)
% length(unique(RIDvec_adnimerge(~isnan(AGEvec_DEM)))) 446
% size(tmp1) = 3190,1. Time to dementia diagnosis for all visits min(DXvec(ivec))<2 & mi = min(find(DXvec(ivec)==2));
% length(unique(RIDvec_adnimerge(~isnan(AGEmaxvec_notDEM)))) 2329
% isnan 1956

y = NaN(size(RIDvec_adnimerge));
% Age at Dem- Age at each visit;
tmp1 = AGEvec_DEM-AGEvec_adnimerge; y(isfinite(tmp1)) = tmp1(isfinite(tmp1));       % These are all non-positive, Time to dementia diagnosis
tmp2 = AGEmaxvec_notDEM-AGEvec_adnimerge; y(isfinite(tmp2)) = tmp2(isfinite(tmp2)); % These are all non-negative, Time to follow-up without dementia

sexvec = SEXvec_adnimerge;
agevec = AGEvec_adnimerge;

T = y; % Time from observation to conversion / progression

convvec = isfinite(AGEvec_DEM); % conversion vector
censorvec = ~convvec; % censored vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PCA
modellist = {
  {[agevec] 'age'} % 1
  {[MMSCORE] 'MMSE'} % 2
  {[measmat_nq] 'IHS'} % 3
  %{[AMYLOID_STATUS] 'AMYLOID_STATUS'} % 4
  %{[AMYLOID_CENTI] 'AMYLOID_CENTI'} % 5
  %{[APOE4sumvec APOE2sumvec] 'APOE'} % 6
  %{[AMYLOID_SUVR] 'AMYLOID_SUVR'} %
  %{[CSF_PTAU181_ABETA42] 'CSF_PTAU181_ABETA42'} %
  %{[MOCA] 'MOCA'} %
  %{[RAVLT] 'RAVLT'} %
  %{[RAVLT_DEL] 'RAVLT_DEL'} %
  %{[APOE4sumvec APOE2sumvec PHSvec_EADB_noAPOE] 'Genetics'} %
  %{[CLOCKSCOR] 'CLOCKSCOR'} %
  %{[TRAASCOR] 'TRAASCOR'} %
  %{[TRABSCOR] 'TRABSCOR'} %
  %{[BNTTOTAL] 'BNTTOTAL'} %
  %{[CDRSB] 'CDRSOB'} %  1
  %{[FAQTOTAL] 'FAQTOTAL'} %  1
  %{[CLOCK_ALL] 'CLOCK_ALL'}
  %{[COPYSCOR] 'COPYSCOR'}
  %{[COPY_ALL] 'COPY_ALL'}
  %{[BNT_ALL] 'BNT_ALL'}
  %{[APOE4sumvec APOE2sumvec] 'APOE'}
  %{[PHSvec_EADB_noAPOE] 'PHS_EADB_noAPOE'}
  %{[PHSvec_Desikan_noAPOE] 'PHS_Desikan_noAPOE'}
  %{[pT217_AB42_F] 'plasma_pT217_AB42_F'}
  %{[PTAU_ABETA42] 'CSF_PTAU_ABETA42'}
};


% constrain validation to observations available across all models
validation_defvec = true(size(RIDvec_adnimerge));
for modeltype = 1:length(modellist)
  measmat = modellist{modeltype}{1};
  validation_defvec = validation_defvec & isfinite(sum(measmat, 2));
end

% Try computing scores for each domain first, then try combinations
%  age   ihs   genet. mmse clock tra  trb    BN    RAVLT RAVLT_DEL CDR  FAQ
incmat=...
[   1   0   0
    1   1   0
    1   0   1
    1   1   1
]


%incmat=incmat(end, :);

% generate all possible model subset combinations
nmodels = length(modellist);
%incmat = [];
%for msum = 1:length(modellist)
%  incmat = cat(1,incmat,fliplr(dec2bin(sum(nchoosek(2.^(0:nmodels-1),msum),2)) - '0'));
%end
%incmat = eye(nmodels);
assert(size(incmat, 2) == nmodels)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% actual cox PROP
% creates a row vector named agevals that contains values starting from 0 to 15 with an increment of 0.1.
agevals = 0:0.1:15;
%prctile_list = [1 5 20 50 80 95 99];
prctile_list = [5 20 50 80 95];
critpct_list = 0;

for pcti = 2:length(prctile_list)+1
  critpct_list(pcti) = 2*prctile_list(pcti-1)-critpct_list(pcti-1);
end

sexlist = {[0] [1] [0 1]}; sexstr = {'Females' 'Males' 'Both Sexes'};
%sexlist = { [0 1]}; sexstr = {'Both Sexes'};


% only models of interest
%incmat = incmat(1:length(modellist),:);


rng("default"); % Re-sets random number generator, to ensure consistency

nrep = 1; k = 10;
predvecs = repmat({NaN(size(agevec))},[nmodels length(sexlist)]);
predvecs_combo = repmat({NaN(size(agevec))},[size(incmat,1) length(sexlist) nrep]);
oddsmat = NaN(size(incmat,1),length(sexlist),length(prctile_list),nrep); % Change to prctile_list from
%fh_roc=sfigure(300); clf;
fh_pca=sfigure(299); clf;
% we create CVind before removing NaN


fprintf(1, 'AUCvec_now and AUCvec_5y gives AUC of model predictin now (at the point of observtion) and in 5 years into the future, for three scenarios: DEM vs non-DEM; DEM vs MCI, DEM vs CN.\n')

for repi = 1:nrep
  cvind_subj = crossvalind('Kfold',length(RIDlist),k);
  cvind = NaN(size(RIDvec_adnimerge));
  for ki = 1:k
    cvind(ismember(RIDvec_adnimerge,RIDlist(cvind_subj==ki))) = ki;
  end
  for ki = 1:k
    for sexi = 1:length(sexlist)
      for modeltype = 1:length(modellist)
        % disp(modellist{modeltype}{2})

        % Perform singular value decomposition (SVD) and extracts principal components.
        measmat = modellist{modeltype}{1};
        X = select_first_PCs(measmat, cvind~=ki, 40, fh_pca);
        if length(sexlist{sexi})>1 % Include Sex as covariate if both sexes included in fit
          X = cat(2,sexvec,X);
        end
        defvec = isfinite(sum(X,2));
        includevec_fit = defvec&T>0&DXvec<2&ismember(sexvec,sexlist{sexi}); % Include all subjects / timpoints without DEM diagnosis
        % includevec = defvec&T>0&DXvec==1&ismember(sexvec,sexlist{sexi}); % Restrict to subjects / timpoints with MCI diagnosis
        ivec_disc = find(includevec_fit&cvind~=ki); % ensures that not in the same CV
        ivec_repl = find(includevec_fit&cvind==ki);
        [b,logl,H,stats] = coxphfit(X(ivec_disc,:),T(ivec_disc),'censoring',censorvec(ivec_disc));
        %if size(X, 2) > 41
        %    b = cvglmnet_R(X(ivec_disc,:),T(ivec_disc),censorvec(ivec_disc), 10, 0.1);
        %    fprintf(1, '%i/%i non-zero betas\n', sum(b~=0), length(b));
        %end
        predvecs{modeltype,sexi} = X*b; % apply weights to entire cohort
      end
      for combi = 1:size(incmat,1)
        predmat = [];
        for j = find(incmat(combi,:))
          predmat = cat(2,predmat,predvecs{j,sexi});
        end
        X = predmat;
        if length(sexlist{sexi})>1 % Include Sex as covariate if both sexes included in fit
          X = cat(2,sexvec,X);
        end
        defvec = isfinite(sum(X,2));
        includevec_fit = defvec&T>0&DXvec<2&ismember(sexvec,sexlist{sexi}); % Include all subjects / timpoints without DEM diagnosis
%       includevec = defvec&T>0&DXvec==1&ismember(sexvec,sexlist{sexi}); % Restrict to subjects / timpoints with MCI diagnosis
        ivec_disc = find(includevec_fit&cvind~=ki); ivec_repl = find(includevec_fit&cvind==ki);
        [b,logl,H,stats] = coxphfit(X(ivec_disc,:),T(ivec_disc),'censoring',censorvec(ivec_disc));
        xvec = H(:,1); fvec = H(:,2);
        predvecs_combo{combi,sexi,repi}(cvind==ki) = X(cvind==ki,:)*b;
      end
    end
  end
  for combi = 1:size(incmat,1)
    jlist = find(incmat(combi,:));
    modelname = modellist{jlist(1)}{2};
    for ji = 2:length(jlist)
      modelname = sprintf('%s + %s',modelname,modellist{jlist(ji)}{2});
      % disp(modelname)
    end
%    sfigure(fh_roc); subplot(length(sexlist),size(incmat,1),(sexi-1)*size(incmat,1)+combi); AUC = plot_ROC(colvec(predvec(defvec&T<=5)),colvec(predvec(defvec&tmp2>=5))); h=title(sprintf('%s in %s (AUC = %0.3f)',modelname,sexstr{sexi},AUC),'Interpreter','none'); axis image; % Not sure if this is correct
    ColorOrder = get(gca,'ColorOrder'); legends = {};
    fh_survival = sfigure(200+combi); clf;
    for sexi = 1:length(sexlist)
      predvec = predvecs_combo{combi,sexi,repi};
      predvec(~ismember(sexvec,sexlist{sexi})) = nan;
      AUCvec = ADNI_MHS_report_AUC(predvec, DXvec, AGEvec_adnimerge, AGEvec_DEM, AGEmaxvec_notDEM, AGEvec_MCI, AGEmaxvec_CN);

      for cnt = 1:length(prctile_list)
        predvec = predvecs_combo{combi,sexi,repi};
        defvec = isfinite(predvec);
        [sv si] = min(abs(agevals-5));
        % includevec = defvec&T>0&DXvec<2&ismember(sexvec,sexlist{sexi}); % Include all subjects / timpoints without DEM diagnosis
        if 0
          includevec = defvec&T>0&DXvec==1&ismember(sexvec,sexlist{sexi}); % Restrict to subjects / timpoints with MCI diagnosis
        else
          includevec = defvec&T>0&DXvec==1&ismember(sexvec,sexlist{sexi}); % Restrict to subjects who would qualify for Lecanemab at timepoint in question -- Iris to check with Jimbo / Doug on this?
%          includevec = defvec&T>0&DXvec==1&ismember(sexvec,sexlist{sexi})&AMYLOID_STATUS==0; % Restrict to subjects who would *not* qualify for Lecanemab at timepoint in question -- Iris to check with Jimbo / Doug on this?
        end

        % Uncomment to enable consistent set of subjects used for
        % validation across all models; this doesn't seem to work well if plasma & CSF
        % are included in modellist (too few observations)
        %includevec = includevec & validation_defvec;

        prctvalvec = prctile(predvec(includevec),prctile_list);
        critvalvec = prctile(predvec(includevec),critpct_list);
        legends{cnt} = sprintf('%4s Risk %%ile',num2ordinal(prctile_list(cnt)));
        ivec_tmp = includevec & predvec>=critvalvec(cnt) & predvec<=critvalvec(cnt+1);
        % disp(cnt)
        [f1,x1]=ecdf(T(ivec_tmp),'censoring',censorvec(ivec_tmp),'function','survivor'); [sv si] = min(abs(agevals-5));
        if length(f1)>0
          try
            [dummy IA] = unique(x1); x1 = x1(IA); f1 = f1(IA);
            s = interp1(x1,f1,5,'linear','extrap'); oddsmat(combi,sexi,cnt,repi) = (1-s)/s;; % Should look into pooling x1 and f1 across reps?
            sfigure(fh_survival); subplot(length(sexlist),1,sexi); plot(x1,f1,'LineWidth',3,'color',ColorOrder(cnt,:)); xlim([min(agevals) max(agevals)]); hold on; ylim([0 1]);
          catch ME
            PrintErrorStruct(ME);
          end
        end
      end
      ageCutoff = 5;
      % DEM vs non-DEM; DEM vs MCI, DEM vs CN
      includevec = defvec&T>0&DXvec<2&ismember(sexvec,sexlist{sexi}); % Restrict to subjects / timpoints with non-AD diagnosis
      [cIndex, odds_ratio, hazard_ratio] = phs_validate(predvec(includevec), T(includevec), 1-censorvec(includevec), ageCutoff);

      includevec = defvec&T>0&DXvec==1&ismember(sexvec,sexlist{sexi}); % Restrict to subjects / timpoints with MCI diagnosis
      [cIndex2, odds_ratio2, hazard_ratio2] = phs_validate(predvec(includevec), T(includevec), 1-censorvec(includevec), ageCutoff);

      includevec = defvec&T>0&DXvec==0&ismember(sexvec,sexlist{sexi}); % Restrict to subjects / timpoints representing healthy controls
      [cIndex3, odds_ratio3, hazard_ratio3] = phs_validate(predvec(includevec), T(includevec), 1-censorvec(includevec), ageCutoff);

      includevec = defvec&T>0&DXvec==1&ismember(sexvec,sexlist{sexi})&AMYLOID_STATUS==1; % Restrict to subjects who would qualify for Lecanemab at timepoint in question -- Iris to check with Jimbo / Doug on this?
      [cIndex4, odds_ratio4, hazard_ratio4] = phs_validate(predvec(includevec), T(includevec), 1-censorvec(includevec), ageCutoff);

      OR_80_20 = mean(oddsmat(combi,sexi,find(prctile_list==80),1:repi),4)./mean(oddsmat(combi,sexi,find(prctile_list==20),1:repi),4);
      h=xlabel('Time (yrs)'); set(h,'FontSize',18,'FontWeight','bold'); h=ylabel('Probability of Dementia free'); set(h,'FontSize',16,'FontWeight','bold');
      h=title(sprintf('%s in %s (OR80/20 = %0.1f)',modelname,sexstr{sexi},OR_80_20),'Interpreter','none'); set(h,'FontSize',18,'FontWeight','bold');
      h=legend(legends,'Location','SW'); set(h,'FontSize',10,'FontWeight','bold'); drawnow;
      fprintf(1, 'combi=%d/%d sexlist=[%s] repi=%d/%d cIndex=[%s], OR80/20=[%s], %i events, %i censored, AUCvec_now=[%s], AUCvec_5y=[%s] - %s\n', ...
                 combi,size(incmat,1),num2str(sexlist{sexi}, '%i '),repi,nrep,num2str([cIndex cIndex2 cIndex3 cIndex4], '%.4f '), ...
                 num2str([odds_ratio odds_ratio2 odds_ratio3], '%.2f '), ...
                 sum(~censorvec(includevec)), sum(censorvec(includevec)), num2str(AUCvec(1:3)', '%.3f '), num2str(AUCvec(4:6)', '%.3f '), modelname);
    end
  end
end

% ToDos
%   Operationalize Lecanemab candidate definition (exclude homozygous E4)
%   Check performance of MHS in subjects with neuropath data -- with vs. without AD pathology ; also visualize as "Glass Brain"
%   Save results & plots -- Iris

