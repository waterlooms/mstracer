# MSTracer Workflow
 
### Code Setup
Clone & navigate to the repository.
```
git clone https://github.com/waterlooms/ms-tracer.git
```

Download scikit-learn(version 0.22) & tensorflow using python3.
```
pip3 install -U scikit-learn==0.22 scipy matplotlib
pip3 install tensorflow
```

### Usage 
Put your .mzXML or .mzML file in ``` data/ ``` or any designated folder.
Your output files will be generated in that folder as well.

Then you could simply run a bash file in terminal:
```
cd ms-tracer/
bash ms-tracer.bash data/file.mzXML
```
```file.mzXML_feature``` is the final output.

#### Explanation

You could also execute without using bash.
In details, the program does the following steps:

##### Step 1, initial detection of peptide features.
```
java -jar ms-tracer_detect.jar data/file.mzXML 
```
```file.mzXML_feature_all_z``` is the output in this step.

##### Step 2, SVR model prediction.
```
python3 SVR.py -feature data/file.mzXML
```
```file.mzXML_svr_score``` is the output in this step.

##### Step 3, charge state selection from the SVR model prediction.
```
java -jar ms-tracer_select.jar data/file.mzXML
```
```file.mzXML_feature_one_z``` is the output in this step.

##### Step 4, NN model prediction.
```
python3 NN.py -feature data/file.mzXML
```
```file.mzXML_nn_score``` is the output in this step.

##### Step 5, ranking from NN model prediction and finalization.
```
java -jar ms-tracer_final.jar data/file.mzXML 
```
```file.mzXML_feature``` is the final output.
