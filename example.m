%% test on room impulse response

path = '../data/RIR/';
[rir, fs] = audioread([path,'K217-ir.wav']);
r = 100;  %downsampling factor
opt_flag = 0;   %optimization flag
room_flag = 1;  %RIR flag
[mode_params, rirhat] = frequency_warped_modal(rir,fs,[],opt_flag,room_flag);
% [mode_params, rirhat] = frequency_zoomed_modal(rir,fs,[],r,opt_flag,room_flag);


% hear results
soundsc(rir,fs);pause(2);
soundsc(rirhat,fs);

% plot spectrograms
ftgram([rir, rirhat],fs,'rir');

%% test on piano note

path = '../data/MIS/';
room_flag = false;

[ir,fs] = audioread([path,'Piano.mf.C1.aiff']);

% preprocessing
%convert to mono
if (size(ir,2)) > 1
    ir = ir(:,1);
end

% truncate signal
duration = 4.0;
start = round(0.5*fs); %to remove noisy transient and silence in beginning
finish = round(duration*fs);
ir = ir(start:finish);
f0 = 32.7;

[mode_params, irhat] = frequency_zoomed_modal(ir,fs,f0,r,opt_flag, room_flag);

% hear results
soundsc(ir,fs);pause(2);
soundsc(irhat,fs);

% plot spectrograms
ftgram([ir, irhat],fs,'music');

