import java.util.ArrayList;


public class Motif implements Comparable<Motif>{
	public Motif(String nm, int nreg) {
		name = nm;
		isPos = name.endsWith("+");
//		regex = name.substring(0, name.length()-1).replaceAll("_", "[ACGT]");
		regex = name.substring(0, name.length()-1);
		prob = new double[nreg];
	}
	
	public String name;
	public boolean isPos;
	public String regex;
	// probabilities for region inclusion with gaps 
	public double[] prob;
	
	@Override
	public int compareTo(Motif o) {
		return name.compareTo(o.name);
	}
	
	@Override
	public String toString() {
		return name;
	}
}
