package edu.uw.waterlooms;
import edu.uw.waterlooms.entity.*;
import edu.uw.waterlooms.ms1.*;

import java.io.*;
import java.nio.file.FileSystems;
import java.nio.file.Path;
import java.util.*;

import edu.uw.waterlooms.msutil.OpenMzxml;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

import java.util.stream.Collectors;

import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.util.Pair;
import org.json.JSONArray;
import org.json.JSONObject;

/**
 * MS-TRACER Main Class
 *
 * @author Xiangyuan ZengF
 */
public class Main {

  public static void main(String[] args) throws IOException {

    // Parse environment variable set in Dockerfile to validate if the executable is running within
    // the container
    boolean containerExecution = Boolean.parseBoolean(System.getenv("RUN_WITHIN_CONTAINER"));
    ArrayList<String> processBuilderCommand = new ArrayList<String>();
    if (!containerExecution) {
      processBuilderCommand.add("docker-compose");
      processBuilderCommand.add("exec");
      processBuilderCommand.add("-T");
      processBuilderCommand.add("mstracer");
    }
    processBuilderCommand.add("python");

    Path path = FileSystems.getDefault().getPath("").toAbsolutePath();
    final String DOCKER_WORKING_DIR = "/data/";
    final String LOCAL_WORKING_DIR = containerExecution ? "/data/" : path + "/data/";

    ArgumentParser parser =
        ArgumentParsers.newFor("java -jar mstracer.jar")
            .build()
            .defaultHelp(true)
            .description("DIA Analysis.");
    parser
        .addArgument("-mzXML")
        .metavar("FILE")
        .type(Arguments.fileType().acceptSystemIn().verifyExists())
        .help("mzXML File");
    parser.addArgument("-parameters").type(String.class).help("JSON String Dict of Parameters");
    parser
        .addArgument("-outputDir")
        .type(String.class)
        .help("Directory for outputting the result CSV.");

    // Parse Arguments
    Namespace ns = null;
    try {
      ns = parser.parseArgs(args);
    } catch (ArgumentParserException $exception) {
      parser.handleError($exception);
      System.exit(1);
    }

    String mzXMLInFile = ns.getString("mzXML");
    String fastaInFile = ns.getString("fasta");
    String outputDir = ns.getString("outputDir");
    String serializedParameters = ns.getString("parameters");

    //String rawFileName = "r01_dia_data.mzXML";
    String rawFileName = "toy.mzXML";
    if ((mzXMLInFile != null && !mzXMLInFile.isEmpty())){
      // TODO: DIA-WEBAPP submits mzXMLInFile as a path/file.mzXML
      // TODO: Need to strip this and set rawFileName
      rawFileName = mzXMLInFile;
    }

    // Set defaults from main/resources/ folder if arguments not specified (development purposes)
    String mzXMLFile =
            (mzXMLInFile != null && !mzXMLInFile.isEmpty())
                    ? mzXMLInFile
                    : LOCAL_WORKING_DIR + rawFileName;
    String outDir =
            (outputDir != null && !outputDir.isEmpty())
                    ? outputDir
                    : LOCAL_WORKING_DIR + "output/";

    // Sanity check for file existence(s)
    File mzXMLIn = new File(mzXMLFile);
    if (!mzXMLIn.exists()) {
      parser.printHelp();
      System.exit(1);
    }

    StopWatch stopWatch = new StopWatch();
    stopWatch.start();

    ///* MS1 Feature Detection
    // Step 1: Identify possible precursors using Xiangyuan's MS1 component
    OpenMzxml of = new OpenMzxml(mzXMLFile);
    FeatureDetect featureDetect = new FeatureDetect(of, FeatureDetect.DetectionType.MS1);
    featureDetect.detectFeatures(mzXMLFile);

    // write to mzXMLFile_svr_score
    System.out.println("Running SVR.py ...");
    ProcessBuilder processBuilder = new ProcessBuilder();
    ArrayList<String> svrCommand = new ArrayList<>(processBuilderCommand);
    svrCommand.add("/mstracer/src/main/python/SVR.py");
    svrCommand.add("-features");
    svrCommand.add(DOCKER_WORKING_DIR + rawFileName);
    processBuilder.command(svrCommand);
    try {
      // TODO: Suppress output if necessary for python
      Process process = processBuilder.inheritIO().start();
      int exitCode = process.waitFor();
      System.out.println("\nSVR.py exited with code : " + exitCode);
    } catch (InterruptedException $e) {
      $e.printStackTrace();
    }

    // write to mzXMLFile_feature_one_z
    System.out.println("Selecting charge with FeatureSelect ...");
    FeatureSelect featureSelect = new FeatureSelect();
    featureSelect.selectFeature(
            LOCAL_WORKING_DIR + rawFileName);

    // Step 4 : write to mzXMLFile_nn_score
    ArrayList<String> nnCommand = new ArrayList<>(processBuilderCommand);
    nnCommand.add("/mstracer/src/main/python/NN.py");
    nnCommand.add("-features");
    nnCommand.add(DOCKER_WORKING_DIR + rawFileName);

    processBuilder.command(nnCommand);
    try {
      Process process = processBuilder.start();
      int exitCode = process.waitFor();
      System.out.println("\nNN.py exited with code : " + exitCode);
    } catch (InterruptedException e) {
      e.printStackTrace();
    }

    // Step 5 : write to mzXMLFile_feature
    // Updated to just save as a list of SVRScores
    List<SVRScore> svrScores;
    svrScores =
            featureSelect.finalizeFeature(
                    LOCAL_WORKING_DIR + rawFileName,
                    LOCAL_WORKING_DIR, rawFileName
            );
    System.out.println("Completed MSTracer Precursor Detection for " + rawFileName + "...");

    stopWatch.stop();
    System.out.println("Elapsed Time in Minutes: " + stopWatch.getTime() + " ...");
    System.exit(0);
  }
}
