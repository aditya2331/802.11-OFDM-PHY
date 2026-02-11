# 802.11-OFDM-PHY

This repository contains the files required to simulate the Tx/Rx stack of 802.11 OFDM. It involves:
- LTS/STS based packet detection and sycnhronization equalization, demodulation.
- Coarse and fine carrier frequency offset (CFO) correction
- Sampling frequency offset (SFO) and residual phase offset correction.
- CP removal and channel estimation, equalization
- Symbol detection and demodulation

There is also SIMO and MIMO implementations using maximum ratio combining (SIMO) and SVD based precoding (MIMO), designed for single user testing.
