Supplementary Data Contents

All models are stored in .xml format because the models are coded in sbml which is based on .xml.

All python notebooks can be accessed and run using Jupiter Notebook. It may also be possible to use google collab to run the code, but the code was originally implemented and run using Jupiter Notebook.

S1 Dataset - contains the python script used to generate the models and an example of the input excel spreadsheets that the script takes as input. This also contains a text file with instructions on how to run the code. The example input is the full model.

S2 Dataset - contains the models for testing how increasing the number of reactions in the model influence the length of loading and simulating the model. These models have incrementing number of reactions. These are the models used in S3 Dataset when comparing simulation times of models with different numbers of reactions.

S3 Dataset - Python notebook for timing loading and simulation times of the models in S2 Dataset. The script also creates the S1 Fig from the data.

S4 Dataset - contains the model of each Vesicular Stomatitis Virus with eyeballed parameters. These are the models that are the starting point for S5 Dataset python notebook. These are .

S5 Dataset - contains the 3 notebooks used to fit the models. The models in the S4 Dataset are necessary to run the Fitting_1_ODE_Model_Fitting.ipynb notebook. This notebook outputs models that are then fed into the Fitting_2_Stoch_Model_1.ipynb python notebook. This then outputs models that are used by the Fitting_3_Stoch_Model_2_and_Data_Generation.ipynb to fit the model and get the final stochastic model used in this paper. These stochastic models are stored in the S6 Dataset file. This also generates and saves the dataset of trajectories. All of the output models from the 3 model fitting notebooks as well as this dataset are required run the python notebook in the S7 Dataset. 

S6 Dataset - Contains the final stochastic models that were fit using the notebooks in supplementary data 5.

S7 Dataset - The notebook used to generate the rest of the figures in the paper except for the parameter scan images. To run all of the code, this notebook requires all models generated in the S5 Dataset and the dataset of 10000 simulations generated using the final stochastic model. Creates S2 Fig and Fig 2. If you want to just generate the figures, it is possible to do so with only the dataset and the final stochastic model if you just do not run the code blocks that load and simulate the other models.

S8 Dataset - The notebook used to generate predictions for all 120 potential variants and then perform the parameter scan and generate the figures associated with it. This notebook just requires the final stochastic models in S6 Dataset. Creates Fig 3, Fig 4, and Fig 5.

