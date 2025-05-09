# ActiveGM
A Digital Twin of a Green Machine built with MZI's and internal modulators to compensate for errors.

A useful demo to see how the MZIs are affected by internal errors can be seen [here](https://www.desmos.com/calculator/ewmf8b2woe). a, b represent the input and ouput couplers, while p is a phase modulator inside of the MZI and q is a phase modulator on one output arm. The final model is composed of many of these MZIs.

Here's an example correction that was possible with this code, using biases of 40/60 coupling (a, b), random internal phase errors between each layer on the order of 0-pi/6, each layer losses of .2 dB, additional modulator loss of .3 dB, and small additional random errors on each bidirectional coupler (a, b) with a standard deviation of 1% coupling. The two images show power coupling between input and output ports before and after modulator corrections, respectively

![Pre Correction](https://github.com/BYUCamachoLab/ActiveGM/raw/main/Errors.png?raw=true)
![Post Correction](https://github.com/BYUCamachoLab/ActiveGM/raw/main/Corrected.png?raw=true)

### tst.py
Creates a digital twin with random imperfections. Uses an iterative optimizaiton method to correct these errors using the on chip modulators. Run this script (takes about an hour on my laptop)

### load_sparams.py
Creates the cirucit structure, defines MZI equations, etc.

### utils.py
Useful functions to build an programmable chip / green machine (using a Bokun or a Clements model, see [this](https://opg.optica.org/oe/fulltext.cfm?uri=oe-31-15-23851&id=532505) paper). Also includes quality functions to evaluate performance and functions to process s parameters from simphony.

