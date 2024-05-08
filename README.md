# Microbiome_Reversions

Code and data repository for the manuscript "Reversions mask the contribution of adaptive evolution in microbiomes" by Paul Torrillo and Tami Lieberman. Simulation code is meant to be run on a Mac operating system. On a system like Windows, it appears numpy.random.poisson and numpy.random.binomial only take up to 32 bit integers rather than 64 bits for the mean. I do not have an easy workaround but if you are on Windows and you have errors with those functions this is likely the cause. 

The current repository should contain all code and data used in the manuscript. Code repository may still be updated to improve readability.

In each figure file, one should find the script that plots it, and if relevant, the script that runs the simulation for the plot. There should also be pre generated simulations outputs to reproduce the runs from the manuscript.

Many parts of the code are very similar to one and other as they are twists and modifications on other simulations or figures. 

There is also a simulation folder that contains a more commented version of the main simulations used in the paper. The original simulations should be in the respective figure files. These may be helpful for those interested in simulating large population sizes. The main idea of what was done is also layed out in the paper's methods section.

Please inquire at torrillo@mit.edu for any questions or assistance
