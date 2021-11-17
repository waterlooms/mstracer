# MSTracer Workflow

ms-tracer version 1.1.0
 
 ### Quick Start
* [Docker](https://docs.docker.com/engine/install/ubuntu/) & [Docker Compose](https://docs.docker.com/compose/install/)  
Download and install Docker.    
```
sudo apt-get update
sudo apt-get install -y apt-transport-https ca-certificates curl software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add - 
sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"
sudo apt-get update
sudo apt-get install docker-ce

sudo groupadd docker
sudo usermod -aG docker $USER
```  
Reboot your computer to apply the new group permissions at this point.  
Install Docker-Compose.  
```
sudo curl -L "https://github.com/docker/compose/releases/download/1.27.4/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose

sudo chmod +x /usr/local/bin/docker-compose
```
### Code Setup
Clone & navigate to the repository.
```
git clone https://github.com/waterlooms/mstracer.git
cd ms-tracer
```

### Usage 
Put your .mzXML or .mzML file in ``` data/ ``` or any designated folder.
Your output files will be generated in that folder as well.

In the commands, replace ```your-file-name.mzXML``` with your file name, or try with a preloaded data ```sample.mzXML```.

The output file is```your-file-name_mstracer.tsv```.

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
By default, the program will execute on the *toy.mzXML* file located in *dia_data_reading/data/*.  
JAVA source files should go under *dia_data_reading/src/main/java/edu/uw/waterlooms/*.  
PYTHON source files should go under *dia_data_reading/src/main/python/*.



#### Troubleshooting
##### Note1
Should there be any issue with Step 4, another way is to run ```NN.py``` using PyCharm. 

Right click the "tensorflow" at the third line, choose "Show Context Actions"; choose "install package tensorflow". This sets up a virtual environment that runs Tensorflow.

In the configuration, (1) choose Python from Templates. (2)Set the script path to that of ```NN.py``` (e.g.```home/Desktop/ms-tracer/src/main/python/NN.py```). (3) Parameters should be "-feature data/your-file-name.mzXML" 

##### Note2
If you encounter any issues please feel free to reach by emailing x25zeng@uwaterloo.ca.
