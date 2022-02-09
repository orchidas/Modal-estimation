<h2>Modal estimation toolbox</h2>

<p> Modes characterize resonant systems; estimating them accurately is tricky. Each mode of a system has a unique complex amplitude, a damping factor and a frequency. This toolbox contains algorithms to find modes from audio signals. Methods included are 

<ul>
	<li> <b>Frequency Zoomed-ARMA</b> - Karjalainen, Matti, Paulo AA Esquef, Poju Antsalo, Aki Mäkivirta, and Vesa Välimäki. "Frequency-zooming ARMA modeling of resonant and reverberant systems." Journal of the Audio Engineering Society 50, no. 12 (2002): 1012-1029.</li>
	<li> <b>Frequency Zoomed-ESPRIT </b>- Kereliuk, Corey, Woody Herman, Russel Wedelich, and Daniel J. Gillespie. "Modal analysis of room impulse responses using subband ESPRIT." In Proceedings of the International Conference on Digital Audio Effects. 2018. </li>
	<li> <b>Frequency Warped-ESPRIT </b>- Proposed by the author. </li>
</ul>
</p>


<p> Frequency-warped ESPRIT helps in estimating beating modes in the low frequency spectrum. Subband ESPRIT or frequency-zoomed ESPRIT replaces warping with bandwise filtering and downsampling, but is slower than FW-ESPRIT. FZ-ESPRIT is very similar to FZ-ARMA in spirit, but replaces ARMA modeling with ESPRIT, which is more robust. I'd recommend using FZ-ESPRIT when a large number of modes need to be estimated, such as in room impulse responses. </p>

<p> Time-domain sequential optimization, or frequency-domain constrained pole optimization further tunes the estimated modal parameters. The constrained pole optimizer has not been tested yet.</p>

See <a href = "example.m">example.m</a> for an example on how to use the toolbox. </p>