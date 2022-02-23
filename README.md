# MSTracer -- Now Available on Docker Hub!

MSTracer version 1.1.0
 
 ### Quick Start
* [Docker](https://docs.docker.com/engine/install/ubuntu/) download and install docker
* Run application:
```
docker run --rm -v $PWD:/app -w /app trackerrr/mstracer java -jar /mstracer/mstracer.jar -mzXML [PATH_TO_MZXML]
```
Note: Put your .mzXML file in ``` data/ ``` or any designated folder.
Your output files will be generated in that folder as well.

In the commands, replace ```[PATH_TO_MZXML]``` with the path to your file, or try with a preloaded data ```toy.mzXML```.

The output file is```your-file-name_precursors.tsv```.

### Developer Usage 

Clone & navigate to the repository.
```
git clone https://github.com/waterlooms/mstracer.git
cd mstracer
```
####Code Setup
Structure at dia_data_reading/ should look like the following:  

```
dia_data_reading/
└───.github/workflows/
└───data/
└───scripts/
└───src/
│   └───main/
│      └───java/edu/uw/waterlooms/
│      └───python/
│      └───resources/
│   └───test/
.dockerignore
.gitignore
Dockerfile
Makefile
docker-compose.yml
pom.xml
README.md
output_file_instruction.txt
```

Your data files should be located in the *data/* folder.  
By default, the program will execute on the *toy.mzXML* file located in *mstracer/data/*.  
JAVA source files should go under *mstracer/src/main/java/edu/uw/waterlooms/*.  
PYTHON source files should go under *mstracer/src/main/python/*.



#### Troubleshooting
##### Note1
Should there be any issue with Step 4, another way is to run ```NN.py``` using PyCharm. 

Right click the "tensorflow" at the third line, choose "Show Context Actions"; choose "install package tensorflow". This sets up a virtual environment that runs Tensorflow.

In the configuration, (1) choose Python from Templates. (2)Set the script path to that of ```NN.py``` (e.g.```home/Desktop/ms-tracer/src/main/python/NN.py```). (3) Parameters should be "-feature data/your-file-name.mzXML" 

##### Note2
If you encounter any issues please feel free to reach by emailing x25zeng@uwaterloo.ca.
