# Hopfield models and derivatives
* You can find [here](https://github.com/mGolos/OldCodeSamples/blob/master/TVB/models.py#L3676) the continuous Hopfield model implementation (last version [here](https://github.com/the-virtual-brain/tvb-library/blob/trunk/tvb/simulator/models/hopfield.py)).  
It was generalized in order to also run my own model and two other derivatives ([equations](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004644#pcbi.1004644.e006)).  
* I had to add another coupling function [here](https://github.com/mGolos/OldCodeSamples/blob/master/TVB/coupling.py#L412) corresponding to an existing one because the simulator wasn't allowing me to use the sigmoid according to the equations' order (last version [here](https://github.com/the-virtual-brain/tvb-library/blob/trunk/tvb/simulator/coupling.py#L393)).
I made it configurable for the four kind of situations.
