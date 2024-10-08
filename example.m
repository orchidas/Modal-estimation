%% test on room impulse response
addpath('Toolbox/');
path = 'test_data/';
[rir, fs] = audioread([path,'K217-ir.wav']);
r = 20;  %downsampling factor
opt_flag = 0;   %optimization flag
room_flag = 1;  %RIR flag
plot_flag = 1;  % plot fitting in each subband

% run modal estimation
[mode_params, rirhat] = frequency_zoomed_modal(rir,fs,[],r,opt_flag,room_flag,plot_flag);

% resynthesize signal using parallel biquads
rirhat_biquad = resynthesize_signal(mode_params, length(rir), fs, 'use_parallel_biquads', true);

% hear results
soundsc(rir,fs);pause(2);
soundsc(rirhat_biquad,fs);


%% test on piano note

room_flag = false;
opt_flag = 1;

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

[mode_params, irhat] = frequency_warped_modal(ir, fs, f0, opt_flag, room_flag);
irhat_biquad = resynthesize_signal(mode_params, duration*fs, fs, 'use_parallel_biquads', true);

% hear results
soundsc(ir,fs);pause(4);
soundsc(irhat_biquad,fs);


