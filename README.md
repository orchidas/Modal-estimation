<h2>Modal estimation with frequency-warped ESPRIT</h2>

<p> Frequency-warped ESPRIT helps in estimating beating modes in the low frequency spectrum. A time-domain optimization scheme further tunes the estimated modal parameters. Subband ESPRIT or frequency-zoomed ESPRIT replaces warping with bandwise filtering and downsampling, but is slower than FW-ESPRIT. I'd recommend using FZ-ESPRIT when a large number of modes need to be estimated, such as in room impulse responses. 

See <a href = "test_RIR.m">test_RIR.m</a> for an example on how to use the toolbox. </p>