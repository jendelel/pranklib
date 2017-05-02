package cz.siret.prank.lib.utils;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.awt.Color;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.Scanner;
import java.util.zip.GZIPInputStream;


public class Utils {
    private final transient Logger logger = LoggerFactory.getLogger(getClass());

    public static void copyStream(InputStream in, OutputStream out) throws IOException {
        byte[] buffer = new byte[1024];
        int len = in.read(buffer);
        while (len != -1) {
            out.write(buffer, 0, len);
            len = in.read(buffer);
        }
    }

    public static InputStream readFile(File file) throws IOException {
        if (file.getName().endsWith(".gz")) {
            return new GZIPInputStream(new FileInputStream(file));
        } else {
            return new FileInputStream(file);
        }
    }

    // http://stackoverflow.com/questions/309424/read-convert-an-inputstream-to-a-string
    // The stream is not closed.
    public static String convertStreamToString(InputStream is) {
        Scanner s = new Scanner(is).useDelimiter("\\A");
        return s.hasNext() ? s.next() : "";
    }

    public static void stringToFile(String content, File destination) throws FileNotFoundException {
        PrintWriter writer = new PrintWriter(destination);
        writer.print(content);
        writer.close();
    }

    // http://stackoverflow.com/questions/340209/generate-colors-between-red-and-green-for-a-power-meter
    private static Color getColor (double power) {
        double H = power * 0.5;
        return Color.getHSBColor((float)H, (float)0.9, (float)0.9);
    }

}
