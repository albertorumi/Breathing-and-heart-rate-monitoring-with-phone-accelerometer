# Breathing and heart rate monitoring with phone accelerometer

In this project I've implemented in matlab different approaches and models to detect the breathing and heart rate of a human with the phone accelerometer. I've used a OnePlus 5T and the application Accelerometer analyzer to record the data. The experiments are taken on me while laing down over a period of 5 minutes with fixed different breathing periods.

The algorithms implemented for the breathing rate estimation are based on minimum detection and Lomb-Scargle power spectral density, as discussed [here](http://hdl.handle.net/11311/1045105).

While the algorithms implemented for the heart rate detection are based on thresholding for peak extraction ([algo1](https://ieeexplore.ieee.org/document/6091301), [algo2](https://pubmed.ncbi.nlm.nih.gov/19794234/)), with also the possibility to apply the Hilbert transform ([algo3](https://pubmed.ncbi.nlm.nih.gov/27681033/)).
