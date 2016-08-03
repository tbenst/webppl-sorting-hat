// Based off of code written by Hubert Eichner
// Accessed At http://myselph.de/hodgkinHuxley.html
// Used under the MIT License

var Tensor = require('adnn/tensor');

function constructHHWaveform(samplingRate, C, GKMax, GNaMax, Gm,
						   EK, ENa, VRest, V0) {
	if (samplingRate % 200 !== 0) {
		throw "Unsupported sampling rate" 
	} else {
		var waveformLength = samplingRate / 200
		var stepsPerSample = 500 / waveformLength
	}

	var hh = hodgkinHuxley(C, GKMax, GNaMax, Gm, EK, ENa, VRest, V0)
	var voltages = Tensor.zeros(waveformLength)
	var j = 0
	for (var i = 0; i < 500; i++) { // 5ms @ .01ms step
		var v = hh.step()
		if (i % stepsPerSample === 0) {
			voltages[j] = v
			j += 1
		}
	}

	return voltages
}

function hodgkinHuxley(C, GKMax, GNaMax, Gm, EK, ENa, VRest, V0) {
	this.dt = 0.01; //ms
	this.C = C; //in muF/cm^2
	this.GKMax = GKMax;
	this.GNaMax = GNaMax;
	this.EK = EK; //mV
	this.ENa = ENa;
	this.Gm = Gm;
	this.VRest = VRest;
	this.V = V0;
	this.n = 0.32;
	this.m = 0.05;
	this.h = 0.60;
	this.INa = 0;
	this.IK = 0;
	this.Im = 0;
	this.Iinj = 0;
	   
	function alphaN(V) {
		if (V===10) return alphaN(V+0.001); // 0/0 -> NaN
		return (10-V) / (100*(Math.exp((10-V)/10)-1));
	};
	function betaN(V) {
		return 0.125 * Math.exp(-V/80);
	};

	function alphaM(V) {
		if (V===25) return alphaM(V+0.001);  // 0/0 -> NaN
		return (25-V) / (10 * (Math.exp((25-V)/10)-1));
	};
	function betaM(V) {
		return 4 * Math.exp(-V/18)
	};

	function alphaH(V) {
		return 0.07*Math.exp(-V/20);
	};
	function betaH(V) {
		return 1 / (Math.exp((30-V)/10)+1);
	};

	this.step = function () {
	var aN = alphaN(this.V);
	var bN = betaN(this.V);
	var aM = alphaM(this.V);
	var bM = betaM(this.V);
	var aH = alphaH(this.V);
	var bH = betaH(this.V);

	var tauN = 1 / (aN + bN);
	var tauM = 1 / (aM + bM);
	var tauH = 1 / (aH + bH);
	var nInf = aN * tauN;
	var mInf = aM * tauM;
	var hInf = aH * tauH;

	this.n += this.dt / tauN * (nInf - this.n);
	this.m += this.dt / tauM * (mInf - this.m);
	this.h += this.dt / tauH * (hInf - this.h);
	this.INa = this.GNaMax * this.m * this.m * this.m * this.h * (this.ENa - this.V);
	this.IK = this.GKMax * this.n * this.n * this.n * this.n * (this.EK - this.V);
	this.Im = this.Gm * (this.VRest - this.V);

	// this.V += this.dt / this.C * (this.INa + this.IK + this.Im + this.Iinj);
	return this.dt / this.C * (this.INa + this.IK + this.Im + this.Iinj);


module.exports = {
  // Adjust exports here
  constructHHWaveform: constructHHWaveform
}
