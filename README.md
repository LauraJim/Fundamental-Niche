# Estimating the Fundamental Niche
Functions used to estimate a Fundamental Niche with the model presented by Jiménez et al. 2019:
 
Jiménez, L., Soberón, J., Christen, J. A., & Soto, D. (2019). On the problem of modeling a fundamental niche from occurrence data. Ecological Modelling, 397, 74-83.

There are two R files with code and two csv files with an example of the datasets needed to apply the model.
- Nf_Model_functions.R contains all the functions needed to apply our model and some functions for plotting the data and the results of the model. This file does not need to be changed except if you want to generalize the model to more than two dimensions, in this case, you may need to do some changes. All the code includes comments that describe what each function does.
- Nf_Model_application.R contains an example of how to use each of the functions. The code is adapted to the datasets included in this example (Background.csv, SpD_50_20), so, when using your own data, you will need to make some changes. In particular, if you want to include the tolerance ranges of the species for each environmental variable, change the parameters (a1,b1,a2,b2) defined in the section "A priori tolerance limits".
- Background.csv contains random points inside the area of study with their corresponding environmental values. It has four columns. The first two columns contain the coordinates of the points and the second two columns contains values of BIO1 (annual mean temperature) and BIO12 (annual precipitation) from WorldClim. The area of study for this example is the Americas and there are 20,000 random points representing this region.
- SpD_50_20.csv contains 50 presence points of a virtual species. It has the same columns than Background.csv.
