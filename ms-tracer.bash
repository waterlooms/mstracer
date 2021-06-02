#!bin bash
java -jar ms-tracer_detect.jar $1
python3 src/main/python/SVR.py -feature $1
java -jar ms-tracer_select.jar $1
python3 src/main/python/NN.py -feature $1
java -jar ms-tracer_final.jar $1

filepath=$(dirname $(readlink -f $1))
filename=$(basename -- "$1")
cp $1_feature "${filepath}/${filename%.*}_mstracer".tsv