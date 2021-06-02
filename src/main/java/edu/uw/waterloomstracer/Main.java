package edu.uw.waterlooms;

import edu.uw.waterlooms.entity.SVRScore;
import edu.uw.waterlooms.entity.XICEntry;
import edu.uw.waterlooms.match.*;
import edu.uw.waterlooms.ms1.FeatureDetect;
import edu.uw.waterlooms.ms1.FeatureSelect;
import edu.uw.waterlooms.ms2.MS2FeatureDetection;
import edu.uw.waterlooms.msutil.OpenMzxml;
import edu.uw.waterlooms.service.PeptideScoreService;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.util.*;

import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.stream.Collectors;

import org.apache.commons.io.FilenameUtils;

/**
 * MSTracer Main Class
 *
 * @author Xiangyuan Zeng, Jia Wu
 */
public class Main {

  public static void main(String[] args) throws IOException {

    Path path = FileSystems.getDefault().getPath("").toAbsolutePath();
    final String LOCAL_WORKING_DIR = path + "/data/";

    ArgumentParser parser =
            ArgumentParsers.newFor("java -jar ms-tracer.jar")
                    .build()
                    .defaultHelp(true)
                    .description("MS-Tracer Peptide Detection.");
    parser
            .addArgument("-mzXML")
            .metavar("FILE")
            .type(Arguments.fileType().acceptSystemIn().verifyExists())
            .help("mzXML File");
    parser
            .addArgument("-detectionParams")
            .metavar("FILE")
            .type(Arguments.fileType().acceptSystemIn().verifyExists())
            .help("FeatureDetect Parameters File");
    parser
            .addArgument("-selectionParams")
            .metavar("FILE")
            .type(Arguments.fileType().acceptSystemIn().verifyExists())
            .help("FeatureSelect Parameters File");

    // Parse Arguments
    Namespace ns = null;
    try {
      ns = parser.parseArgs(args);
    } catch (ArgumentParserException $exception) {
      parser.handleError($exception);
      System.exit(1);
    }

    String mzXMLInFile = ns.getString("mzXML");
    String selectionParamsInFile = ns.getString("selectionParams");
    String detectionParamsInFile = ns.getString("detectionParams");

    String mzXMLFile =
            (mzXMLInFile != null && !mzXMLInFile.isEmpty())
                    ? mzXMLInFile
                    : LOCAL_WORKING_DIR + "sample.mzXML";
    String selectionParams =
            (selectionParamsInFile != null && !selectionParamsInFile.isEmpty())
                    ? selectionParamsInFile
                    : LOCAL_WORKING_DIR + "featureselect.params";
    String detectionParams =
            (detectionParamsInFile != null && !detectionParamsInFile.isEmpty())
                    ? detectionParamsInFile
                    : LOCAL_WORKING_DIR + "featuredetect.params";

    // Sanity check for file existence(s)
    File mzXMLIn = new File(mzXMLFile);
    if (!mzXMLIn.exists()) {
      parser.printHelp();
      System.exit(1);
    }
    String mzXMLFname = FilenameUtils.getName(mzXMLFile);

    // Read in Parameters for Feature Detection & Feature Selection
    Properties featureDetectParams = new Properties(), featureSelectParams = new Properties();
    try {
      FileReader fin = new FileReader(detectionParams);
      featureDetectParams.load(fin);
      fin = new FileReader(selectionParams);
      featureSelectParams.load(fin);
    } catch (IOException $e) {
      System.err.printf($e.getMessage());
      System.exit(1);
    }
    // Step 1 : write to mzXMLFile_feature_all_z
    OpenMzxml of = new OpenMzxml(mzXMLFile);
    FeatureDetect featureDetect = new FeatureDetect(of, featureDetectParams);
    featureDetect.detectFeatures(mzXMLFile);
    // Step 2 : write to mzXMLFile_svr_score
    ProcessBuilder processBuilder = new ProcessBuilder();
    ArrayList<String> svrCommand = new ArrayList<>(processBuilderCommand);
    svrCommand.add("/waterloomstracer/src/main/python/SVR.py");
    svrCommand.add("-features");
    svrCommand.add(LOCAL_WORKING_DIR + mzXMLFname + "_feature_all_z");

    processBuilder.command(svrCommand);
    try {
      // TODO: Suppress output if necessary for python
      Process process = processBuilder.inheritIO().start();
      int exitCode = process.waitFor();
      System.out.println("\nSVR.py exited with code : " + exitCode);
    } catch (InterruptedException $e) {
      $e.printStackTrace();
    }

    // Step 3 : write to mzXMLFile_feature_one_z
    FeatureSelect featureSelect = new FeatureSelect();
    featureSelect.selectFeature(
            LOCAL_WORKING_DIR + mzXMLFname + "_svr_score", featureSelectParams);

    // Step 4 : write to mzXMLFile_nn_score
    ArrayList<String> nnCommand = new ArrayList<>(processBuilderCommand);
    nnCommand.add("/waterloomstracer/src/main/python/NN.py");
    nnCommand.add("-features");
    nnCommand.add(LOCAL_WORKING_DIR + mzXMLFname + "_feature_one_z");
    processBuilder.command(nnCommand);
    try {
      Process process = processBuilder.start();
      int exitCode = process.waitFor();
      System.out.println("\nNN.py exited with code : " + exitCode);
    } catch (InterruptedException e) {
      e.printStackTrace();
    }

    // Step 5 : write to mzXMLFile_feature
    featureSelect.finalizeFeature(
            LOCAL_WORKING_DIR + mzXMLFname + "_nn_score", featureSelectParams);
  }
}


