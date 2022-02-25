<h2>Modal estimation toolbox</h2>

<p> Modes characterize resonant systems; estimating them accurately is tricky. Each mode of a system has a unique complex amplitude, a damping factor and a frequency. This toolbox contains algorithms to find modes from audio signals. Methods included are 

<ul>
	<li> <b>Frequency Zoomed-ARMA</b> </li>
	<li> <b>Frequency Zoomed-ESPRIT </b> </li>
	<li> <b>Frequency Warped-ESPRIT </b> </li>
</ul>
</p>


<p> Frequency-warped ESPRIT helps in estimating beating modes in the low frequency spectrum. Subband ESPRIT or frequency-zoomed ESPRIT replaces warping with bandwise filtering and downsampling, but is slower than FW-ESPRIT. FZ-ESPRIT is very similar to FZ-ARMA in spirit, but replaces ARMA modeling with ESPRIT, which is more robust. I'd recommend using FZ-ESPRIT when a large number of modes need to be estimated, such as in room impulse responses. </p>

<p> Time-domain sequential optimization, or frequency-domain constrained pole optimization further tunes the estimated modal parameters. The constrained pole optimizer has not been tested yet.</p>

See <a href = "example.m">example.m</a> for an example on how to use the toolbox. </p>

<h3> References </h3>
<p>
<ul>
	<li> Das, Orchisama and Jomathan S. Abel. "Modal Estimation on a Warped Frequency Axis for Linear System Modeling." arXiv Preprint, 2022 <a href = "http://arxiv.org/abs/2202.11192
">http://arxiv.org/abs/2202.11192</a> </li>
	<li>  Karjalainen, Matti, Paulo AA Esquef, Poju Antsalo, Aki Mäkivirta, and Vesa Välimäki. "Frequency-zooming ARMA modeling of resonant and reverberant systems." Journal of the Audio Engineering Society 50, no. 12 (2002): 1012-1029. </li>
	<li> Kereliuk, Corey, Woody Herman, Russel Wedelich, and Daniel J. Gillespie. "Modal analysis of room impulse responses using subband ESPRIT." In Proceedings of the International Conference on Digital Audio Effects. 2018. </li>
	<li> Maestre, Esteban, Jonathan S. Abel, Julius O. Smith, and Gary P. Scavone. "Constrained pole optimization for modal reverberation." In Proc. of the 5th International Conference on Digital Audio Effects (DAFx). 2017. </li>

</ul>
</p>