package edu.uw.waterlooms.msutil;

import org.json.JSONArray;
import org.json.JSONObject;
import umich.ms.datatypes.scan.IScan;
import umich.ms.datatypes.spectrum.ISpectrum;
import umich.ms.fileio.exceptions.FileParsingException;
import java.nio.ByteBuffer;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

public class Scans {
  public OpenMzxml openFile = null;
  public Set<Map.Entry<Integer, IScan>> scanEntries = null;
  public TreeMap<Integer, IScan> num2scan = null;

  public Scans(OpenMzxml of, Integer mslevel) {
    if (mslevel == 1) {
      openFile = of;
      scanEntries = of.scanEntries;
      num2scan = of.num2scan;
    } else {
      openFile = of;
      scanEntries = of.scanEntries2;
      num2scan = of.num2scan2;
    }
  }

  public Scans(OpenMzxml of, Set<Map.Entry<Integer, IScan>> se, TreeMap<Integer, IScan> ns) {
    openFile = of;
    scanEntries = se;
    num2scan = ns;
  }

  // write lcms scans as Json strings
  private int totoalpeaks = 0;
  //    public Set<JSONObject> getLCMS() {
  //        Set<JSONObject> specs = new HashSet<JSONObject>();
  //        for (Map.Entry<Integer, IScan> scanEntry : scanEntries) {
  //            IScan scan = scanEntry.getValue();
  //            JSONObject spectrumJObj = scan2Json(scan);
  //            specs.add(spectrumJObj);
  //        }
  //        return specs;
  //    }

  public JSONObject scan2Json(IScan scan) {
    JSONObject spectrumJObj = new JSONObject();
    spectrumJObj.put("sN", scan.getNum()); // scanNum
    // childScans
    JSONArray childScanJArr = new JSONArray();
    if (scan.getChildScans() != null && !scan.getChildScans().isEmpty()) {
      for (int child : scan.getChildScans()) {
        JSONArray childScan = new JSONArray();
        childScan.put(child);
        childScan.put(openFile.scans.getScanByNum(child).getPrecursor().getMzTarget());
        childScanJArr.put(childScan);
      }
    }
    spectrumJObj.put("cS", childScanJArr); // child scan num and targetmz
    spectrumJObj.put("rt", scan.getRt()); // rt
    // basepeak
    spectrumJObj.put("bpmz", scan.getBasePeakMz());
    spectrumJObj.put("bpint", scan.getBasePeakIntensity());
    spectrumJObj.put("tic", scan.getTic()); // tic
    spectrumJObj.put("wL", scan.getScanMzWindowLower()); // scanMzWindowLower
    spectrumJObj.put("wU", scan.getScanMzWindowUpper()); // scanMzWindowUpper
    if (scan.getMsLevel() == 2) {
      spectrumJObj.put("pcsmz", scan.getPrecursor().getMzTarget()); // precursor
    }
    ISpectrum spectrum = null;
    try {
      spectrum = scan.fetchSpectrum();
    } catch (FileParsingException e) {
      e.printStackTrace();
    }
    double[] mz = spectrum.getMZs();
    double[] intensities = spectrum.getIntensities();
    int peakNum = mz.length;
    totoalpeaks += peakNum;
    ByteBuffer byteBufferMZ = ByteBuffer.allocate(Double.BYTES * peakNum);
    ByteBuffer byteBufferInt = ByteBuffer.allocate(Double.BYTES * peakNum);
    for (int i = 0; i < peakNum; i++) {
      byteBufferMZ.putDouble(mz[i]);
      byteBufferInt.putDouble(intensities[i]);
    }
    //            byteBuffer.position(0);
    String mzstr = Utils.byte2Base64StringFun(byteBufferMZ.array());
    String intstr = Utils.byte2Base64StringFun(byteBufferInt.array());
    //            byte[] backr = Utils.base64String2ByteFun(mzstr);//test
    spectrumJObj.put("mz", mzstr);
    spectrumJObj.put("int", intstr);
    return spectrumJObj;
  }
}
