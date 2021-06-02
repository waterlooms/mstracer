package edu.uw.waterlooms.msutil;

import org.apache.commons.codec.binary.Base64;

public class Utils {
  // base64 String to byte[]
  public static byte[] base64String2ByteFun(String base64Str) {
    return Base64.decodeBase64(base64Str);
  }
  // byte[] to base64
  public static String byte2Base64StringFun(byte[] b) {
    return Base64.encodeBase64String(b);
  }
}
