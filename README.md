# MSTracer Workflow
 
### Code Setup
Clone & navigate to the repository.
```
git clone https://github.com/waterlooms/ms-tracer.git
```
Make sure your python version is later than 3.6.
```
python3 --version # check version
```
If not downloaded, please visit https://python.org/download.

Download scikit-learn(version 0.22) & tensorflow(1.14.0) using python3.
```
pip3 install -U scikit-learn==0.22 scipy matplotlib
pip3 install tensorflow==1.14.0
pip3 install 'h5py==2.10.0' --force-reinstall
```

### Usage 
Put your .mzXML or .mzML file in ``` data/ ``` or any designated folder.
Your output files will be generated in that folder as well.

In the commands, replace ```your-file-name.mzXML``` with your file name, or try with a preloaded data ```toy.mzXML```.

Then you could simply run a bash file in terminal:
```
cd ms-tracer/
bash ms-tracer.bash data/your-file-name.mzXML
```
```your-file-name.mzXML_feature``` is the final output.

#### Explanation

You could also execute without using bash.
In details, the program does the following steps:

##### Step 1, initial detection of peptide features.
```
java -jar ms-tracer_detect.jar data/your-file-name.mzXML 
```
```your-file-name.mzXML_feature_all_z``` is the output in this step.

##### Step 2, SVR model prediction.
```
python3 SVR.py -feature data/your-file-name.mzXML
```
```your-file-name.mzXML_svr_score``` is the output in this step.

##### Step 3, charge state selection from the SVR model prediction.
```
java -jar ms-tracer_select.jar data/your-file-name.mzXML
```
```your-file-name.mzXML_feature_one_z``` is the output in this step.

##### Step 4, NN model prediction.
```
python3 NN.py -feature data/your-file-name.mzXML
```
```your-file-name.mzXML_nn_score``` is the output in this step.

##### Step 5, ranking from NN model prediction and finalization.
```
java -jar ms-tracer_final.jar data/your-file-name.mzXML 
```
```your-file-name.mzXML_feature``` is the final output.

#### Trouble shooting
##### Note1
Should there be any issue with Step 4, another way is to run ```NN.py``` using PyCharm. 

Right click the "tensorflow" at the third line, choose "Show Context Actions"; choose "install package tensorflow". This sets up a virtual environment that runs Tensorflow.

In the configuration, (1) choose Python from Templates. (2)Set the script path to that of ```NN.py``` (e.g.```home/Desktop/ms-tracer/NN.py```). (3) Parameters should be "-feature data/your-file-name.mzXML" 

##### Note2
If you encounter any issues please feel free to reach by emailing x25zeng@uwaterloo.ca.
