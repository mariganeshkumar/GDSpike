# GDSPIKE: AN ACCURATE SPIKE ESTIMATION ALGORITHM FROM NOISY CALCIUM FLUORESCENCE SIGNALS #

Accurate estimation of spike train from calcium (Ca 2+ ) fluorescence signals is challenging owing to significant fluctuations of fluorescence level. This method is a non-model-based approach for spike train inference using group delay (GD) analysis. It primarily exploits the property that changes in Ca 2+ fluorescence corresponding to a spike can be characterized by an onset, an attack, and a decay. 


### How do I get set up? ###

* Clone the repository with **"git clone https://mariganeshkumar@bitbucket.org/mariganeshkumar/gdspike.git"**
* Use the function GDspike.m to estimate the spike train from GCamp fluorescence traces.


### What are the inputs needed for GDspike.m? ###
You need to give the following two arguments, in order
* Calcium fluorescence signal - 1D time series.
* Sampling_rate - Sampling rate of calcium signal.

### Whom to contact? ###

* Write to jilt18 [at] gmail.com for any clarifications.