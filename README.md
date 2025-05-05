# ActiveGM
A Digital Twin of a Green Machine built with MZI's and internal modulators to compensate for errors.

### tst.py
Creates a digital twin with random imperfections. Uses an iterative optimizaiton method to correct these errors using the on chip modulators. Run this script (takes about an hour on my laptop)

### load_sparams.py
Creates the cirucit structure, defines MZI equations, etc.

### utils.py
Useful functions to build an programmable chip / green machine (using a Bokun or a Clements model, see [this](https://opg.optica.org/oe/fulltext.cfm?uri=oe-31-15-23851&id=532505) paper). Also includes quality functions to evaluate performance and functions to process s parameters from simphony.
