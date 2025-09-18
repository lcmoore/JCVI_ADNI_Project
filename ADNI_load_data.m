if ~exist('forceflag','var') || isempty(forceflag)
  forceflag = false;
end

% add matlab paths -- should not be needed
addpath(genpath('/home/dale/matlab'));
addpath("/space/broce-syn01/1/data/GWAS/xinwang/scripts") % AMD: What's needed from here? -- Iris: move functions needed to this directory

% datasets 2025-01-31
tbl_demo = readtable('/space/ceph/1/CMIG/ADNI/LAKE/Demog/PTDEMOG.csv');
tbl_dx = readtable('/space/ceph/1/CMIG/ADNI/LAKE/Diagnosis/DXSUM.csv');
tbl_mmse = readtable('/space/ceph/1/CMIG/ADNI/LAKE/Cognition/MMSE.csv');
tbl_moca = readtable('/space/ceph/1/CMIG/ADNI/LAKE/Cognition/MOCA.csv');
tbl_faq = readtable('/space/ceph/1/CMIG/ADNI/LAKE/Cognition/FAQ.csv');
tbl_cdr = readtable('/space/ceph/1/CMIG/ADNI/LAKE/Cognition/CDR.csv');
tbl_neurobat = readtable('/space/ceph/1/CMIG/ADNI/LAKE/Cognition/NEUROBAT.csv');
tbl_tau_csf = readtable('/space/ceph/1/CMIG/ADNI/LAKE/Biomarker/UPENNBIOMK_ROCHE_ELECSYS.csv');
tbl_tau_plasma = readtable('/space/ceph/1/CMIG/ADNI/LAKE/Biomarker/ADNI4_PLASMA_BIOMARKER_DATA.csv');
tbl_apoeres = readtable('/space/ceph/1/CMIG/ADNI/DATA/Genetic/APOE.csv');
tbl_phs_eadb = readtable('/space/ceph/1/CMIG/ADNI/PHS/ADNI_EADB_PHS_20250130.csv');
tbl_phs_desikan = readtable('/space/ceph/1/CMIG/ADNI/PHS/ADNI_DESIKAN_PHS.csv');
tbl_tau = readtable('/space/ceph/1/CMIG/ADNI/LAKE/PET/UCBERKELEY_TAU_6MM.csv');
tbl_taupvc = readtable('/space/ceph/1/CMIG/ADNI/LAKE/PET/UCBERKELEY_TAUPVC_6MM.csv');
tbl_amyloid = readtable('/space/ceph/1/CMIG/ADNI/LAKE/PET/UCBERKELEY_AMY_6MM.csv'); % ib 06/25
tbl_amysuvr = readtable('/space/ceph/1/CMIG/ADNI/LAKE/PET/AMYLOID_DER.csv');
fname_nqlq = '/space/ceph/1/CMIG/ADNI/NQLQ/measures_nqlq.mat';

% parsed NQ and retain more variable (pICV, z-scores, etc)
% good to check if NQ z-scores outperform lh+rh
fname_nqlq_of = '/space/ceph/1/CMIG/ADNI/NQLQ/measures_nqlq_of.mat';

% Rename all Exam Dates to 'EXAMDATE' add PTAU/ABETA42 ratio
tbl_mmse.Properties.VariableNames{'VISDATE'} = 'EXAMDATE';
tbl_moca.Properties.VariableNames{'VISDATE'} = 'EXAMDATE';
tbl_cdr.Properties.VariableNames{'VISDATE'} = 'EXAMDATE';
tbl_faq.Properties.VariableNames{'VISDATE'} = 'EXAMDATE';
tbl_neurobat.Properties.VariableNames{'VISDATE'} = 'EXAMDATE';
tbl_amysuvr.Properties.VariableNames{'SCANDATE'} = 'EXAMDATE'; % added this line
tbl_tau.Properties.VariableNames{'SCANDATE'} = 'EXAMDATE'; % Xin added this line
tbl_taupvc.Properties.VariableNames{'SCANDATE'} = 'EXAMDATE'; % Xin added this line
tbl_tau_csf(:,{'COMMENT','update_stamp'}) = [] % Xin added, to remove unrelated columns
tbl_tau_csf.CSF_PTAU181_ABETA42 = tbl_tau_csf.PTAU./tbl_tau_csf.ABETA42; % Xin changed "CSF_PTAU_ABETA42" to "CSF_PTAU181_ABETA42"
tbl_amyloid.Properties.VariableNames{'SCANDATE'} = 'EXAMDATE'; % Xin added this line

% Xin added these lines below to differentia plasma and csf biomarkers
tbl_tau_csf.Properties.VariableNames{'ABETA40'} = 'CSF_ABETA40';
tbl_tau_csf.Properties.VariableNames{'ABETA42'} = 'CSF_ABETA42';
tbl_tau_csf.Properties.VariableNames{'TAU'} = 'CSF_Total_Tau';
tbl_tau_csf.Properties.VariableNames{'PTAU'} = 'CSF_pTau181';
tbl_tau_plasma.Properties.VariableNames{'pT217_F'} = 'Plasma_pT217';
tbl_tau_plasma.Properties.VariableNames{'AB42_F'} = 'Plasma_AB42';
tbl_tau_plasma.Properties.VariableNames{'AB40_F'} = 'Plasma_AB40';
tbl_tau_plasma.Properties.VariableNames{'AB42_AB40_F'} = 'Plasma_AB42_AB40_F';
tbl_tau_plasma.Properties.VariableNames{'pT217_AB42_F'} = 'Plasma_PTAU217_ABETA42';
% Xin added these lines below to differentia tau pet data without PVC and PVC
tbl_taupvc.Properties.VariableNames{'META_TEMPORAL_SUVR'} = 'PVC_META_TEMPORAL_SUVR'; 
tbl_taupvc.Properties.VariableNames{'CTX_ENTORHINAL_SUVR'} = 'PVC_CTX_ENTORHINAL_SUVR'; 

% Xin added these lines to re-intenstity normalise the suvr data
suvr_jvec = endsWith(tbl_tau.Properties.VariableNames, '_SUVR'); 
tbl_tau{:, suvr_jvec} = tbl_tau{:, suvr_jvec} ./tbl_tau{:,'ERODED_SUBCORTICALWM_SUVR'};
pvcsuvr_jvec = endsWith(tbl_taupvc.Properties.VariableNames, '_SUVR'); 
tbl_taupvc{:, pvcsuvr_jvec} = tbl_taupvc{:, pvcsuvr_jvec} ./tbl_taupvc{:,'CEREBRAL_WHITE_MATTER_SUVR'};


% Remove duplicates from Demographics
[~, uniqueIdx] = unique(tbl_demo.RID, 'stable');  % 'stable' keeps the first occurrence
tbl_demo = tbl_demo(uniqueIdx, :); % Create a table with only the unique RIDs
tbl_demo.DOB = datetime(tbl_demo.PTDOB, 'InputFormat', 'MM/yyyy', 'Format', 'yyyy-MM-dd');

% Add Age and Sex to tbl_dx ; Male = 1 , Female = 2

RIDlist = unique(tbl_dx.RID);
tbl_dx.AGE = NaN(size(tbl_dx.RID));
tbl_dx.SEX = NaN(size(tbl_dx.RID));

for RIDi = 1:length(RIDlist)
 RID = RIDlist(RIDi);
 ivec_dx = find(tbl_dx.RID==RID);
 ivec_demo = find(tbl_demo.RID==RID);
 
 SEXvec = tbl_demo.PTGENDER(ivec_demo);
 DOBvec = tbl_demo.DOB(ivec_demo);

 if length(SEXvec) > 0
   tbl_dx.SEX(ivec_dx) = SEXvec;
 end
 
 if length(DOBvec) > 0
 % Fixed: avoid rounding this to fulf years
   tbl_dx.AGE(ivec_dx) = years(tbl_dx.EXAMDATE(ivec_dx) - DOBvec);
 end
end


tbl_adnimerge = tbl_dx; % replace adni_merge with tbl_dx
%clear tbl_dx
% clear RIDlist
                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Some changes to AMD CODE
% change: Female coding from 2 to 0 to work with code below
% change: DX is coded 1 = CN, 2 = MCI, 3 = Dementia ; DXvec remains the same

RIDvec_adnimerge = tbl_adnimerge.RID;
RIDvec_apoeres = tbl_apoeres.RID;
RIDvec_tau = tbl_tau.RID; VISCODEvec_tau = tbl_tau.VISCODE;
RIDvec_taupvc = tbl_taupvc.RID; VISCODEvec_taupvc = tbl_taupvc.VISCODE;
EXAMDATEvec_adnimerge = tbl_adnimerge.EXAMDATE;
VISCODEvec_adnimerge = tbl_adnimerge.VISCODE;
AGEvec_adnimerge = tbl_adnimerge.AGE;
SEXvec_adnimerge = tbl_adnimerge.SEX;
SEXvec_adnimerge(tbl_adnimerge.SEX==2) = 0 ;
ivec_CN = find(tbl_adnimerge.DIAGNOSIS == 1);
ivec_MCI = find(tbl_adnimerge.DIAGNOSIS == 2);
ivec_DEM = find(tbl_adnimerge.DIAGNOSIS == 3);
DXvec = NaN(size(RIDvec_adnimerge));
DXvec(ivec_CN) = 0;
DXvec(ivec_MCI) = 1;
DXvec(ivec_DEM) = 2;
RIDlist = unique(RIDvec_adnimerge);

EVENTCODEvec_adnimerge = cell(size(DXvec));
for i = 1:length(RIDvec_adnimerge)
  EVENTCODEvec_adnimerge{i} = sprintf('%04d_%s',RIDvec_adnimerge(i),VISCODEvec_adnimerge{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Read in other DSETS and Get ExamDate
% Datasets of interest testing with tbl_dx
%% xin changes the datasets 
% tbl_tau_csf: Xin changed the first var "CSF_PTAU_ABETA42" to "ABETA40" anf last var "CSF_PTAU_ABETA42" to "CSF_PTAU181_ABETA42"
% tbl_tau_plasma: Xin changed the first var 'pT217_AB42_F' to 'Plasma_pT217' and last var "pT217_AB42_F" to 'Plasma_PTAU217_ABETA42'
% tbl_tau: Xin added tau PET data
% tbl_taupvc: Xin added tau PET PVC data
datasets = {
  % The format is as follows:
  % tbl newVarName firstvar lastvar
  {[tbl_adnimerge] 'tbl_adnimerge'}
  {[tbl_mmse] 'MMSCORE' 'MMSCORE' 'MMSCORE'}
  {[tbl_moca] 'MOCA' 'MOCA' 'MOCA'}
  {[tbl_neurobat] 'RAVLT' 'AVTOT1' 'AVENDED'}      % AVTOT1 - AVENDED
  {[tbl_neurobat] 'RAVLT_DEL' 'AVDEL30MIN' 'AVDELERR2'} % AVDEL30MIN - AVDELERR2
  % {[tbl_tau_csf] 'CSF_PTAU_ABETA42' 'CSF_ABETA40' 'CSF_PTAU181_ABETA42'}  %
  {[tbl_tau_csf] 'CSF_PTAU181_ABETA42' 'CSF_PTAU181_ABETA42' 'CSF_PTAU181_ABETA42'}
  {[tbl_tau_plasma] 'Plasma_pT217_AB42_F' 'Plasma_pT217' 'Plasma_PTAU217_ABETA42'}
  {[tbl_tau] 'TAU_SUVR' 'META_TEMPORAL_SUVR' 'CTX_ENTORHINAL_SUVR'}
  {[tbl_taupvc] 'TAU_PVCSUVR' 'PVC_META_TEMPORAL_SUVR' 'PVC_CTX_ENTORHINAL_SUVR'}
  {[tbl_amyloid] 'AMYLOID_STATUS' 'AMYLOID_STATUS_COMPOSITE_REF' 'AMYLOID_STATUS_COMPOSITE_REF'}
  {[tbl_amysuvr] 'AMYLOID_SUVR' 'CTX_CAUDALANTERIORCINGULATE_SUVR' 'HIPPOCAMPUS_SUVR'}
  {[tbl_amyloid] 'AMYLOID_CENTI' 'CENTILOIDS' 'CENTILOIDS'}
  % {[tbl_faq] 'FAQTOTAL' 'FAQTOTAL' 'FAQTOTAL'}
  % {[tbl_cdr] 'CDRSB' 'CDRSB' 'CDRSB'}
  % {[tbl_neurobat] 'CLOCKSCOR' 'CLOCKSCOR' 'CLOCKSCOR'}
  % {[tbl_neurobat] 'CLOCK_ALL' 'CLOCKCIRC' 'CLOCKSCOR'} % CLOCKCIRC - CLOCKSCOR
  % {[tbl_neurobat] 'COPYSCOR' 'COPYSCOR' 'COPYSCOR'}
  % {[tbl_neurobat] 'COPY_ALL' 'COPYCIRC' 'COPYSCOR'} % COPYCIRC - COPYSCOR
  % {[tbl_neurobat] 'TRAASCOR' 'TRAASCOR' 'TRAASCOR'}
  % {[tbl_neurobat] 'TRABSCOR' 'TRABSCOR' 'TRABSCOR'}
  % {[tbl_neurobat] 'BNTTOTAL' 'BNTTOTAL' 'BNTTOTAL'}
  % {[tbl_neurobat] 'BNT_ALL' 'BNTSPONT' 'BNTTOTAL'}  % BNTSPONT - BNTTOTAL
};

%%%% Grab scores closest to EXAMDATE in DX table, recycled imaging code
% Changed: use 6 months, instead of 12 months

for dataset = 2:length(datasets)
  if isa(datasets{dataset}{1}, 'table')
    % Read in data
    tmp_tbl = datasets{dataset}{1};

    newVarName = datasets{dataset}{2};
    first_var = datasets{dataset}{3};
    last_var =  datasets{dataset}{4};
    first_var_idx = find(strcmp(first_var,tmp_tbl.Properties.VariableNames));
    last_var_idx = find(strcmp(last_var,tmp_tbl.Properties.VariableNames));
    if isempty(first_var_idx) | isempty(last_var_idx), error(); end
    fprintf(1, '%s: %s-%s (cols %i-%i)\n', newVarName, first_var, last_var, first_var_idx, last_var_idx);
    var = datasets{dataset}{1}.Properties.VariableNames(first_var_idx:last_var_idx);
    measmat = table2array(tmp_tbl(:,first_var_idx:last_var_idx));
    jvec = find(mean(isfinite(measmat))>0.5*max(mean(isfinite(measmat)))&nanstd(measmat)>0); measmat = measmat(:,jvec);var = var(:,jvec);

    % Generate outcome variable
    outcomeVar = NaN(size(RIDvec_adnimerge,1), size(measmat, 2));

    for RIDi = 1:length(RIDlist)
    RID = RIDlist(RIDi);
    ivec_adnimerge = find(RIDvec_adnimerge==RID);
    [sv si] = sort(AGEvec_adnimerge(ivec_adnimerge),'ascend');
    ivec_adnimerge = ivec_adnimerge(si); % DX table (RID, EXAMDATE, DX)
    ivec_tmp_tbl = find(tmp_tbl.RID==RID);% MMSE table (RID, EXAMDATE, DX)
    if length(ivec_tmp_tbl)>0
     for ii = 1:length(ivec_adnimerge)
       EXAMDATE = EXAMDATEvec_adnimerge(ivec_adnimerge(ii));
       [sv si] = sort(abs(tmp_tbl.EXAMDATE(ivec_tmp_tbl)-EXAMDATE));
       if sv(1) <= duration(365*24/2,0,0) % Within 6 months
         si = si(sv==sv(1));
         result = measmat(ivec_tmp_tbl(si), :);
         outcomeVar(ivec_adnimerge(ii), :) = nanmean(result, 1);
         if ivec_adnimerge(ii)==14000
    %          keyboard
         end
           else
             % disp(RID)
             % disp(EXAMDATE)
             % disp(EXAMDATEvec_nq(ivec_nq)')
           end
       end
     end
    end

    % Use eval to assign the table to the new name
    eval([newVarName ' = outcomeVar ;']);
    tbl_adnimerge(:, var) = array2table(outcomeVar, 'VariableNames', var); % Xin added this line to add variables above to tbl_adnimerge
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create PHSvec
% use PHSvec_EADB_noAPOE, excluding APOE effects
% similarly, use PHSvec_Desikan_noAPOE
% These models make sense only in conjunction with APOE4sumvec APOE2sumvec

RIDvec_phs = cellfun(@(x)str2num(getfield(regexp(x,'S_(?<rid>\d+)$','names'),'rid')),tbl_phs_eadb.IID);
PHSvec_EADB_noAPOE = NaN(size(RIDvec_adnimerge));
for subji = 1:length(RIDlist)
  ind = find(RIDvec_phs==RIDlist(subji));
  if ~isempty(ind)
    ivec_tmp = find(RIDvec_adnimerge==RIDlist(subji));
    PHSvec_EADB_noAPOE(ivec_tmp) = tbl_phs_eadb.EADB_PHS_excl_apoe(ind);
  end
end

% compute APOE-only effect of the Desikan model
desikan_e2_weight = -0.47;     desikan_e4_weight = 1.03;
tbl_phs_desikan.('Desikan_PHS_only_apoe')(strcmp('e2e2',tbl_phs_desikan.apoe_gt)) = desikan_e2_weight * 2;
tbl_phs_desikan.('Desikan_PHS_only_apoe')(strcmp('e2e3',tbl_phs_desikan.apoe_gt)) = desikan_e2_weight * 1;
tbl_phs_desikan.('Desikan_PHS_only_apoe')(strcmp('e3e3',tbl_phs_desikan.apoe_gt)) = 0;
tbl_phs_desikan.('Desikan_PHS_only_apoe')(strcmp('e2e4',tbl_phs_desikan.apoe_gt)) = desikan_e2_weight + desikan_e4_weight;
tbl_phs_desikan.('Desikan_PHS_only_apoe')(strcmp('e3e4',tbl_phs_desikan.apoe_gt)) = desikan_e4_weight;
tbl_phs_desikan.('Desikan_PHS_only_apoe')(strcmp('e4e4',tbl_phs_desikan.apoe_gt)) = desikan_e4_weight * 2;
tbl_phs_desikan.('Desikan_PHS_excl_apoe') = tbl_phs_desikan.('Desikan_PHS') - tbl_phs_desikan.('Desikan_PHS_only_apoe');

RIDvec_desikan = cellfun(@(x)str2num(getfield(regexp(x,'S_(?<rid>\d+)$','names'),'rid')),tbl_phs_desikan.PTID);
PHSvec_Desikan_noAPOE = NaN(size(RIDvec_adnimerge));
for subji = 1:length(RIDlist)
  ind = find(RIDvec_desikan==RIDlist(subji));
  if ~isempty(ind)
    ivec_tmp = find(RIDvec_adnimerge==RIDlist(subji));
    PHSvec_Desikan_noAPOE(ivec_tmp) = tbl_phs_desikan.Desikan_PHS_excl_apoe(ind);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CN, MCI, DEM counts
% counts in comments changed
disp(length(unique(RIDvec_adnimerge(DXvec==0)))) % 1504 subjects to have diagnosis CN at some point
disp(length(unique(RIDvec_adnimerge(DXvec==1)))) % 1611 subjects to have diagnosis MCI at some point
disp(length(unique(RIDvec_adnimerge(DXvec==2)))) % 969 subjects to have diagnosis DEM at some point

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Delete need to read in using CSV file {empty vector}
% no changes made
% Get neurocog measures
measmat_cog = table2array(tbl_adnimerge(:,find(strcmp('CDRSB',tbl_adnimerge.Properties.VariableNames)):find(strcmp('EcogSPTotal',tbl_adnimerge.Properties.VariableNames))));
jvec = find(mean(isfinite(measmat_cog))>0.5*max(mean(isfinite(measmat_cog)))&nanstd(measmat_cog)>0); measmat_cog = measmat_cog(:,jvec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Read in NQ/LQ measures
% changes removed RIDvec_nq(find(find_cell(tmp.IDvec_nq))) = cellfun(@(x)str2num(getfield(regexp(x,'^(?<rid>\d+)_','names'),'rid')),tmp.IDvec_nq(find(find_cell(tmp.IDvec_nq))));
% idx = find(strcmp('002_S_0295', tbl_adnimerge.PTID))
% tbl_adnimerge(idx, {'RID', 'PTID'})
% exam date to be same format as other datasets 'yyyy-MM-dd'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp = load(fname_nqlq); RIDvec_nq = NaN(size(tmp.IDvec_nq));

RIDvec_nq(find(find_cell(tmp.IDvec_nq))) = cellfun(@(x)str2num(getfield(regexp(x,'^(\d)+_S_(?<rid>\d+)','names'),'rid')),tmp.IDvec_nq(find(find_cell(tmp.IDvec_nq))));

tmp.datetimevec_mri_nq.Format = 'yyyy-MM-dd';
EXAMDATEvec_nq = tmp.datetimevec_mri_nq;
% tmp_nq = (tmp.measmat_rh+tmp.measmat_lh); % Use total volume across LH & RH per ROI
tmp_nq = cat(2,tmp.measmat_rh,tmp.measmat_lh); % Use concatenated LH and RH ROIs
rh_measnamevec_nq = cellfun(@(x) strcat('rh_',x), tmp.measnamevec_nq, 'UniformOutput', false);
lh_measnamevec_nq = cellfun(@(x) strcat('lh_',x), tmp.measnamevec_nq, 'UniformOutput', false);
measnamevec_nq = [rh_measnamevec_nq,lh_measnamevec_nq]

tmp.datetimevec_mri_lq.Format = 'yyyy-MM-dd';
EXAMDATEvec_lq = tmp.datetimevec_mri_lq;
tmp_lq = tmp.measmat_lq;
measnamevec_lq = tmp.measnamevec_lq;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NQLQ RID == DX RSID by EXAMDATE, within 6 months
% no changes made

measmat_nq = NaN(size(tbl_adnimerge,1),size(tmp_nq,2));
measmat_lq = NaN(size(tbl_adnimerge,1),size(tmp_lq,2));
for RIDi = 1:length(RIDlist)
  RID = RIDlist(RIDi);
  ivec_adnimerge = find(RIDvec_adnimerge==RID);
  [sv si] = sort(AGEvec_adnimerge(ivec_adnimerge),'ascend'); ivec_adnimerge = ivec_adnimerge(si);
  ivec_nq = find(RIDvec_nq==RID);
  if length(ivec_nq)>0
    for ii = 1:length(ivec_adnimerge)
      EXAMDATE = EXAMDATEvec_adnimerge(ivec_adnimerge(ii));
      [sv si] = sort(abs(datenum(EXAMDATEvec_nq(ivec_nq))-datenum(EXAMDATE)));
      if sv(1) <= 180 % Max half a year difference between imaging and adnimerge assessments
        fprintf(1,'RIDi=%d/%d RID=%d EXAMDATE=%s EXAMDATE_nq=%s,diff=%d\n',RIDi,length(RIDlist),RID,datestr(EXAMDATE,'yyyy-mm-dd'),datestr(EXAMDATEvec_nq(ivec_nq(si(1))),'yyyy-mm-dd'),round(sv(1)));
        si = si(sv==sv(1));
        measmat_nq(ivec_adnimerge(ii),:) = nanmean(tmp_nq(ivec_nq(si),:),1);
        measmat_lq(ivec_adnimerge(ii),:) = nanmean(tmp_lq(ivec_nq(si),:),1);
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get rid of missing columns on nq and lq
% Satisfy two conditions
% 1. At least 90% of the maximum proportion of finite values.
% 2. Nonzero standard deviation.
% Columns 71 and 72 in NQ are removed {'Unclassified', 'Unknown'}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Xin changed: to get the coresponding nq and lq name after getting rid missing columns
jvec = find(mean(isfinite(measmat_nq))>0.9*max(mean(isfinite(measmat_nq)))&nanstd(measmat_nq)>0); measmat_nq = measmat_nq(:,jvec); measnamevec_nq = measnamevec_nq(jvec);
jvec = find(mean(isfinite(measmat_lq))>0.9*max(mean(isfinite(measmat_lq)))&nanstd(measmat_lq)>0); measmat_lq = measmat_lq(:,jvec); measnamevec_lq = measnamevec_lq(:,jvec);

mean(isfinite(measmat_nq)) % 0.76  proportion of rows of measmat_nq with non-finite entries??
mean(isfinite(measmat_lq)) % 0.46  proportion of rows of measmat_nq with non-finite entries??
sum(mean(isfinite(measmat_lq),2)==0) % 7890 rows with NaN
sum(mean(isfinite(measmat_nq),2)==0) % 3379 rows with NaN

% Xin added these line below to add nqlq to tbl_adnimerge
measnamevec_nq = cellfun(@char, measnamevec_nq, 'UniformOutput', false);
measnamevec_lq = cellfun(@char, measnamevec_lq, 'UniformOutput', false);
tbl_nqlq = array2table([measmat_nq measmat_lq], 'VariableNames', [measnamevec_nq measnamevec_lq]);
tbl_adnimerge = [tbl_adnimerge,tbl_nqlq];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure
% no changes made
if 0
  figure(1); clf;
  subplot(1,2,1); imagesc(log10(measmat_nq));
  subplot(1,2,2); imagesc(log10(measmat_lq));
  colormap(hot);
end


% % % % % % % % % % % % % % APOE genotype sum
% no changes size(2758,1)
% % % % % % % % % % % % % %

APOE2sum = (double(tbl_apoeres.APGEN1==2)+double(tbl_apoeres.APGEN2==2));
APOE3sum = (double(tbl_apoeres.APGEN1==3)+double(tbl_apoeres.APGEN2==3));
APOE4sum = (double(tbl_apoeres.APGEN1==4)+double(tbl_apoeres.APGEN2==4));

% % % % % % % % % % % % % %  Classify subjects and compute longitudinal change
% Significant changes
% Output variables
% % % % % % % % % % % % % %

dmeasmat_nq = NaN(size(measmat_nq)); % measure of change

% AGEvec_DEM and AGEmaxvec_notDEM are used in non-dementia (CN or MCI)
% progression to dementia. 
% AGEvec_DEM is used to define time to event 
% AGEmaxvec_notDEM is used to define time to last observation (censored
% observations)
% Observations where patient was already observed with dementia needs to be
% excluded (this is done further down below, see 'includevec_fit' include "T>0")
AGEvec_DEM = NaN(size(AGEvec_adnimerge));
AGEmaxvec_notDEM = NaN(size(AGEvec_adnimerge));

% AGEvec_MCI  and AGEmaxvec_CN are not used, but can be used for CN to MCI
% progression. The meaning is similar to AGEvec_DEM and AGEmaxvec_notDEM
% above.
AGEvec_MCI = NaN(size(AGEvec_adnimerge));
AGEmaxvec_CN = NaN(size(AGEvec_adnimerge));

APOE2sumvec = NaN(size(AGEvec_adnimerge));
APOE3sumvec = NaN(size(AGEvec_adnimerge));
APOE4sumvec = NaN(size(AGEvec_adnimerge));

% not used for calculating
% RIDlist_CN = [];
% RIDlist_DEM = [];
% RIDlist_MCI = [];
% RIDlist_APOE44 = [];
% RIDlist_APOE34 = [];
% RIDlist_APOE33 = [];

RIDlist_CN_stab = []; RIDlist_DEM_stab = []; RIDlist_MCI_stab = []; RIDlist_CN_DEM = []; RIDlist_MCI_prog = []; RIDlist_DEM_prog = []; AGElist_DEM_prog = []; AGElist_MCI_prog = []; RIDlist_APOE44 = []; RIDlist_APOE34 = []; RIDlist_APOE33 = [];

for RIDi = 1:length(RIDlist)
  RID = RIDlist(RIDi);
  ind_apoeres = find(RIDvec_apoeres==RID);
  ivec = find(RIDvec_adnimerge==RID);
  ivec_bl = ivec(find(strcmp('bl',tbl_adnimerge.VISCODE(ivec)))); % Should perhaps allow for screening

  if length(ivec>0)
    [sv si] = sort(AGEvec_adnimerge(ivec),'ascend'); ivec = ivec(si);
    dxvec_tmp = DXvec(ivec); dxvec_tmp = dxvec_tmp(isfinite(dxvec_tmp));
                                
    % RIDlists for diagnosis at first visit
    % NOT USED for calculation
    % if length(dxvec_tmp)>0 && dxvec_tmp(1)==0
    % RIDlist_CN = cat(1,RIDlist_CN,RID);
    % elseif length(dxvec_tmp)>0 && dxvec_tmp(1)==1
    % RIDlist_MCI = cat(1,RIDlist_MCI,RID);
    % elseif length(dxvec_tmp)>0 && dxvec_tmp(1)==2
    % RIDlist_DEM = cat(1,RIDlist_DEM,RID);
    % end
                                
    % Age at highest diagnosis category 1) CN 2) CN/MCI but not Demented
    if max(DXvec(ivec))==0
      AGEmaxvec_CN(ivec) = AGEvec_adnimerge(ivec(end));
    end

    if max(DXvec(ivec))<2
      AGEmaxvec_notDEM(ivec) = AGEvec_adnimerge(ivec(end));
    end

    % Age at diagnosis category 1) convert to MCI 2) convert to DEM
    if min(DXvec(ivec))<1
    mi = min(find(DXvec(ivec)==1));
      if ~isempty(mi)
        AGEvec_MCI(ivec) = AGEvec_adnimerge(ivec(mi));
      end
    end

    if min(DXvec(ivec))<2
    mi = min(find(DXvec(ivec)==2));
      if ~isempty(mi)
      AGEvec_DEM(ivec) = AGEvec_adnimerge(ivec(mi));
      end
    end

    % Include APOE info for those who have it
    if length(ind_apoeres)==1
      APOE2sumvec(ivec) = APOE2sum(ind_apoeres);
      APOE4sumvec(ivec) = APOE4sum(ind_apoeres);
      APOE3sumvec(ivec) = APOE3sum(ind_apoeres);
    
      % not used below
      % if APOE4sum(ind_apoeres)==2
      % RIDlist_APOE44 = cat(1,RIDlist_APOE44,RID);
      % end
      % if APOE4sum(ind_apoeres)==1 && APOE3sum(ind_apoeres)==1
      % RIDlist_APOE34 = cat(1,RIDlist_APOE34,RID);
      % end
      % if APOE3sum(ind_apoeres)==2
      % RIDlist_APOE33 = cat(1,RIDlist_APOE33,RID);
      % end
    end
    
    % not used below
    % if length(ivec_bl)==1
    % AGE_bl = AGEvec_adnimerge(ivec_bl);
    % Years from baseline: assuming they have age at baseline
    % DTvec_bl = AGEvec_adnimerge(ivec)-AGE_bl;
    % end

    % change in NQ
    for ii = 2:length(ivec)
    dtvec = AGEvec_adnimerge(ivec)-AGEvec_adnimerge(ivec(ii));
    ind = find(dtvec<-0.25);
      if ~isempty(ind)
        dmeasmat_nq(ivec(ii),:) = (measmat_nq(ivec(ii),:) - measmat_nq(ivec(ind(end)),:))/(AGEvec_adnimerge(ivec(ii))-AGEvec_adnimerge(ivec(ind(end))));
      end
    end
  end
end

% Xin added these line below for ADNI_norm
RIDlist_CN_stab = unique(RIDvec_adnimerge(isfinite(AGEmaxvec_CN))); 
RIDlist_DEM_prog = unique(RIDvec_adnimerge(isfinite(AGEvec_DEM)));
RIDlist_MCI_prog = unique(RIDvec_adnimerge(isfinite(AGEvec_DEM) & DXvec == 1)); 
DTvec_DEM = NaN(size(AGEvec_adnimerge));
DTvec_DEM(isfinite(AGEvec_DEM)) = AGEvec_adnimerge(isfinite(AGEvec_DEM)) -  AGEvec_DEM(isfinite(AGEvec_DEM));



% AMD added the lines below for vertexwise data from MMPS and for neuropath data
addpath(genpath('/usr/pubsw/packages/cmig'));

tbl_path = readtable('/space/ceph/1/CMIG/ADNI/DATA/NEUROPATH/ADNI_NEUROPATH_Classfication_110.csv');
colnames_path = tbl_path.Properties.VariableNames(4:7);

%fpat_aseg = '/space/ceph/1/CMIG/ADNI/MMPS/MMPS_ADNI*/fsurf/*/mri/aseg.mgz';
fpat_aseg = '/space/ceph/1/CMIG/ADNI/MMPS/MMPS_ADNI*/fsurf/*_*_*/mri/aseg.mgz';
regexp_pat = '^FSURF_(?<site>...)_S_(?<subj>....)_(?<date>[^_]*)_(?<datetime>[^_]*)_1$';
outdir_cache = '~dale/ADNI/cache';

load_mmps_results;

%tmp_lh = load('~/ADNI/cache/thickness-sm256-lh.mat'); tmp_rh = load('~/ADNI/cache/thickness-sm256-rh.mat'); % Do these contain subjids, phases, and sites?
tmp_lh = load('~dale/ADNI/cache/thickness-sm1000-lh.mat'); tmp_rh = load('~dale/ADNI/cache/thickness-sm1000-rh.mat'); % Do these contain subjids, phases, and sites?
%tmp_lh = load('~/ADNI/cache/thickness-sm2819-lh.mat'); tmp_rh = load('~/ADNI/cache/thickness-sm2819-rh.mat'); % Do these contain subjids, phases, and sites?
measmat = cat(2,tmp_lh.measmat,tmp_rh.measmat);

[subjidlist IA IC] = unique(subjids);
sexvec = NaN(size(measmat,1),1); agevec = NaN(size(measmat,1),1); dxvec = NaN(size(measmat,1),1); dxmat_path = NaN(size(measmat,1),4);
for subji = 1:length(subjidlist)
  RID = str2num(subjidlist{subji});
  ivec_demo = find(tbl_demo.RID==RID);
  ivec_adnimerge = find(tbl_adnimerge.RID==RID);
  ivec_path = find(tbl_path.RID==RID);
  if length(ivec_demo)>0 & length(ivec_adnimerge)>0
    ivec_measmat = find(IC==subji);
    sexvec(ivec_measmat) = 2-tbl_demo.PTGENDER(ivec_demo(1));
    agevec(ivec_measmat) = (datenums(ivec_measmat)-datenum(tbl_demo.DOB(ivec_demo)))/365.25;
    dxvec(ivec_measmat) = nanmax(tbl_adnimerge.DIAGNOSIS(ivec_adnimerge));
    if ~isempty(ivec_path)
      dxmat_path(ivec_measmat,:) = repmat(strcmp('YES',table2cell(tbl_path(ivec_path,4:7))),[length(ivec_measmat) 1]);
    end
  else
    fprintf(1,'subji=%d/%d: no match found for %s\n',subji,length(subjidlist),subjidlist{subji});
  end
end

[phaselist IA IC_phase] = unique(phases);
[sitelist IA IC_site] = unique(sites);
X = [];
for phasei = 1:length(phaselist)
  for sitei = 1:length(sitelist)
    ivec_tmp = IC_phase==phasei&IC_site==sitei;
    if sum(ivec_tmp)>9
      X = cat(2,X,ivec_tmp);
    end
  end
end
[sv si] = sort(sum(X),'descend');
X_phase_site = X(:,si(2:end)); % Get rid of most common phase*site



% ToDos
%   Check why this is so slow
