#!bin bash
java -jar ms-tracer_detect.jar $1
python3 SVR.py -feature $1
java -jar ms-tracer_select.jar $1
python3 NN.py -feature $1
java -jar ms-tracer_final.jar $1
