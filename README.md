<h2>Modal estimation with frequency-warped ESPRIT</h2>

<p> Frequency-warped ESPRIT helps in estimating beating modes in the low frequency spectrum. Subband ESPRIT or frequency-zoomed ESPRIT replaces warping with bandwise filtering and downsampling, but is slower than FW-ESPRIT. I'd recommend using FZ-ESPRIT when a large number of modes need to be estimated, such as in room impulse responses. </p>

<p> Time-domain sequential optimization, or frequency-domain constrained pole optimization further tunes the estimated modal parameters. The constrained pole optimizer has not been tested yet.</p>

See <a href = "example.m">example.m</a> for an example on how to use the toolbox. </p>