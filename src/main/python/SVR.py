#!/usr/bin/python
import csv
import pickle
import sys
import os
sys.path.append("/mstracer/src/main/python")

import argparse

class SVR:
    def __init__(self):
        self.ID = 0
        self.mz = 1
        self.rt = 2
        self.z = 3
        self.isonum = 4
        self.int_shape = 5
        self.iso_distr = 6
        self.intensity_area_percentage = 7
        self.scan_num = 8
        self.quantification_peaks_sum = 9
        self.quantification_peaks_area = 10
        self.svr_score = 11
        self.quality_score = 12

    # TO-DO: change to run(self, feature_file, parameter_dict)
    def run(self, feature_file):
        start = self.isonum
        end = self.iso_distr
        filepath = os.path.splitext(feature_file)[0]
        with open(filepath + "_featureAllZ.tsv", 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            test_data = list(reader)
        X_svr = test_data[1: len(test_data)]  # jump header
        for i in range(0, len(X_svr)):
            X_svr[i] = X_svr[i][start:end+1]

        for i in range(0, len(X_svr)):
            for j in range(0, len(X_svr[i])):
                X_svr[i][j] = float(X_svr[i][j])

        path = "/mstracer/src/main/python"
        clf = pickle.load(open(path + "/model/SVR_z_selection", 'rb'))
        predict_svr = clf.predict(X_svr)

        with open(filepath + "_SVRScore.tsv", "w+") as outfile:
            outfile.write("id\tmz\trt\tz\tisotope_num\tintensity_shape_score\tisotope_distribution_score\tintensity_area_percentage\tscan_num\tquantification_peaks_sum\tquantification_peaks_area\tsvr_score\tquality_score\tmzs,rts,ints\n")
            for i in range(len(X_svr)):
                scannum = int(test_data[i + 1][self.scan_num])
                mzsStr = test_data[i + 1][self.quality_score + 1]
                rtsStr =  test_data[i + 1][self.quality_score + 1 + scannum]
                intsStr =  test_data[i + 1][self.quality_score + 1 + scannum * 2]
                for j in range(1, scannum):
                    mzsStr += "\t%s" % test_data[i + 1][self.quality_score + 1 + j]
                    rtsStr += "\t%s" % test_data[i + 1][self.quality_score + 1 + scannum + j]
                    intsStr += "\t%s" % test_data[i + 1][self.quality_score + 1 + scannum * 2 + j]

                outfile.write("%s\t" % test_data[i + 1][self.ID])
                outfile.write("%s\t" % test_data[i + 1][self.mz])
                outfile.write("%s\t" % test_data[i + 1][self.rt])
                outfile.write("%s\t" % test_data[i + 1][self.z])
                outfile.write("%s\t" % test_data[i + 1][self.isonum])
                outfile.write("%s\t" % test_data[i + 1][self.int_shape])
                outfile.write("%s\t" % test_data[i + 1][self.iso_distr])
                outfile.write("%s\t" % test_data[i + 1][self.intensity_area_percentage])
                outfile.write("%s\t" % scannum)
                outfile.write("%s\t" % test_data[i + 1][self.quantification_peaks_sum])
                outfile.write("%s\t" % test_data[i + 1][self.quantification_peaks_area])
                outfile.write("%f\t" % predict_svr[i])
                outfile.write("0\t")
                outfile.write(mzsStr + "\t")
                outfile.write(rtsStr + "\t")
                outfile.write(intsStr + "\n")
        print("Done!")



# Read in arguments
import argparse # add to beginning of file
parser = argparse.ArgumentParser(description='Train on features extracted from mzXML.')
parser.add_argument("-features",
    help="Features computed from FeatureDetect."
    )
parser.add_argument("-parameters",
    default="/data/featuredetect.params",
    help="FeatureDetect parameter file."
    )
args = parser.parse_args()



# Parse the .env like parameter file
# TO-DO uncomment this
# from dotenv import dotenv_values
# parsed = dotenv_values(stream=args.parameters)


# Invoke the training
training = SVR()
# TO-DO: change to training.run(args.features, parsed)
training.run(args.features)