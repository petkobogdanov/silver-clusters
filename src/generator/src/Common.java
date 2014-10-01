import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.TreeMap;

public class Common {
	public static String TMP_DIR = "/dev/shm/";
	public static String PLOT_DIR="/home/petko/Dropbox/plot/";
	public static void plotHist(TreeMap<Double,Integer> h){
		try {
			writeHist(h, Common.TMP_DIR + "hist");
			Runtime run = Runtime.getRuntime() ;
			System.err.print(Common.PLOT_DIR + "plothist " + Common.TMP_DIR + "hist\n");
			Process pr = run.exec(Common.PLOT_DIR + "plothist " + Common.TMP_DIR + "hist");
			pr.waitFor();
		} catch (Exception e) {
			System.err.print("Error Plotting a histogram\n");
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	public static void writeHist(TreeMap<Double,Integer> h, String fn) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(fn));
		for (double i:h.keySet()) {
			bw.write(i + " " + h.get(i) + "\n");
		}	
		bw.close();
	}

}
