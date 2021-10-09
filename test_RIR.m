%% test CCRMA classroom RIR

path = '../data/RIR/';
[rir, fs] = audioread([path,'K217-ir.wav']);
r = 100;  %downsampling factor
opt_flag = 0;   %optimization flag
room_flag = 1;  %RIR flag
[mode_params, rirhat] = frequency_warped_modal(rir,fs,[],opt_flag,room_flag);
% [mode_params, rirhat] = frequency_zoomed_modal(rir,fs,[],r,opt_flag,room_flag);
% [mode_params, rirhat] = frequency_zoomed_ARMA(rir,fs,[],r,opt_flag,room_flag);

RIR = struct('nmodes',size(mode_params,1),'mode_params',mode_params,'rirhat',rirhat);

%%

% hear results
soundsc(rir,fs);pause(2);
soundsc(rirhat,fs);

% plot spectrograms
ftgram([rir, rirhat],fs,'rir');

% save results
save([path,'FW_Modal_opt.mat'],'RIR');
audiowrite([path,'FW_Modal_opt.wav'],rirhat,fs);

% save([path,'FZ_Modal_opt_r=',num2str(r),'.mat'],'RIR');
% audiowrite([path,'FZ_Modal_opt_r=',num2str(r),'.wav'],rirhat,fs);

% save([path,'FZ_Arma_r=',num2str(r),'.mat'],'RIR');
% audiowrite([path,'FZ_Arma_r=',num2str(r),'.wav'],rirhat,fs);


%% plot results
close all;
path = '../data/RIR/';
method = {'FW_Modal','FW_Modal','FZ_Modal','FZ_Arma'};
ext = {'', '_alt', '_r=8', '_r=100'};

[rir,fs] = audioread([path,'K217-ir.wav']);
nmodes = zeros(length(method),1);

for i = 1:length(method)
    [rirhat,~] = audioread([path, method{i}, ext{i}, '.wav']);
    [rirhat_opt,~] = audioread([path, method{i}, '_opt', ext{i}, '.wav']);
    load([path, method{i}, ext{i}, '.mat']);
    nmodes(i) = size(RIR.mode_params,1);
    
    fig = figure('Units','inches', 'Position',[0 0 3.3 4.8],'PaperPositionMode','auto');
    % plot spectrograms
    [~,ax,~,~] = ftgram([rir, rirhat, rirhat_opt],fs,'rir','tanhflag',false);
    lgd = legend(ax(1),{'Measured', method{i}, [method{i},'_opt']},'FontSize',5,'Interpreter','none');
    lgd.ItemTokenSize = [10,30];
    set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');
    print(['../figures/rir_',method{i},ext{i},'.eps'], '-depsc');

end

    
    
    