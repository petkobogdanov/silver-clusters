import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/*
From the distribution of brights (biased to higher Integrated Intensity sequences) 
REPEAT:
1) sample number of dark/bright motifs in parts 1,2,3  
2) average stickiness in parts 1,2,3
// so far we have a "recipe" for a sequence next we use motifs
3) sample motifs according to numbers in step 1) and such we satisfy stickiness from step 2)
REPEAT:
   4) fill in the gaps (multiple times) randomly
   5) ensure the generated sequence is different from our current population
   6) run through classifier and accept if Prob(bright) > Some threshold
*/

public class Generator {
	public static HashMap<String,Double> BASE_STICKY = null;
	static{
		BASE_STICKY = new HashMap<String, Double>();
		BASE_STICKY.put("A", 0.5);
		BASE_STICKY.put("T", 0.0);
		BASE_STICKY.put("C", 1.0);
		BASE_STICKY.put("G", 1.0);
	}
	
	private Random r = new Random(System.currentTimeMillis());
	private int n = 0, nreg = 0;
	private String[] seqs = null;
	private int[][] bcnt;
	private double[][] bsize;
	private double[][] bsticky;
	private int[][] dcnt;
	private String[] classes;
	private ArrayList<Motif> motifs;
	private double[] intensity;
	
	// distributions for sampling
	private ArrayList<TreeMap<Double,Integer>> distBnum = new ArrayList<TreeMap<Double, Integer>>();
	private ArrayList<TreeMap<Double,Integer>> distDnum = new ArrayList<TreeMap<Double, Integer>>();
	private ArrayList<TreeMap<Double,Double>> distSticky = new ArrayList<TreeMap<Double, Double>>();
	private ArrayList<TreeMap<Double,Motif>> distMotif = new ArrayList<TreeMap<Double, Motif>>();
	
	public Generator(String inSeqFn) throws IOException{
		readSequences(inSeqFn);
		computeDists("bright");
	}

	private ArrayList<String> generate(Integer numNew) {
		int cnt = 0;
		int atMostNgaps = 1;
		ArrayList<String> generated = new ArrayList<String>();
		while (cnt < numNew) {
			String seq = "__________";
			seq = placeMotifs(seq, atMostNgaps);
			seq = fillGapsBasedOnStickyness(seq);
			generated.add(seq);
			cnt++;
		}
		return generated;
	}

	private String fillGapsBasedOnStickyness(String seq) {
		// Fill in gaps according to sticky[]
		int ind = -1;
		while ((ind = seq.indexOf('_')) >= 0) {
			ArrayList<Integer> regs = getRegions(ind, seq.length(), nreg);
			TreeMap<String,Double> base2score = new TreeMap<String, Double>();
			double totScore = 0.0;
			for (String base: BASE_STICKY.keySet()) base2score.put(base, 0.0);
			for (Integer reg: regs) {
				int s = reg*(seq.length()/nreg);
				int e = (reg+1)*(seq.length()/nreg) + seq.length()%nreg;
				String str = seq.substring(s,e);
				for (String base: base2score.keySet()) {
					int replace = ind-s;
					assert(str.charAt(replace)=='_');
					String strnew = str.substring(0,replace) + 
							base + str.substring(replace+1,str.length());
					Double avgSticky = computeAvgStickyness(strnew);
					double sc = 0;
					if (null == distSticky.get(reg).floorKey(avgSticky)) {
						sc = distSticky.get(reg).firstKey();
					} else {
						sc = ((null==distSticky.get(reg).ceilingKey(avgSticky))
								?distSticky.get(reg).lastKey()
								:distSticky.get(reg).ceilingKey(avgSticky)) - 
							((null==distSticky.get(reg).floorKey(avgSticky - 0.000000001))?0:
								distSticky.get(reg).floorKey(avgSticky - 0.000000001));
					}
					base2score.put(base, base2score.get(base) + sc);
					totScore += sc;
				}
			}
			TreeMap<Double,String> prob2base = new TreeMap<Double, String>();
			double acc = 0.0;
			for(String base: base2score.keySet()) {
				prob2base.put(base2score.get(base)/totScore + acc, base);
				acc += base2score.get(base)/totScore;
			}
			String base = prob2base.get(prob2base.ceilingKey(r.nextDouble()));
			seq = seq.substring(0,ind) + base + seq.substring(ind+1,seq.length());
		}
		return seq;
	}

	private String placeMotifs(String seq, int atMostNgaps) {
		int[] pos = getShuffledIndices();
		int numGaps = seq.length();
		// sample motifs in 1 then 2 and 3 positions
		while (numGaps > atMostNgaps) {
			for (int j=0;j<nreg;j++) {
				int i = pos[j];
				Motif m = distMotif.get(i).get(distMotif.get(i).ceilingKey(r.nextDouble()));
				ArrayList<Integer> possible_pos = getPossiblePositions(m,seq,i);
				if(possible_pos.size()>0) {
					int position = possible_pos.get(r.nextInt(possible_pos.size()));
					seq = insert(m,seq,position);
				}
			}
			numGaps = 0;
			for (int i=0; i < seq.length(); i++)
				if (seq.charAt(i) == '_') 
					numGaps++;
		}
		return seq;
	}

	private static Double computeAvgStickyness(String strnew) {
		int cnt = 0;
		double sticky = 0.0;
		for (int i = 0; i < strnew.length(); i++) {
			if (strnew.charAt(i)!='_') {
				cnt++;
				sticky += BASE_STICKY.get(strnew.substring(i, i+1));
			}
		}
		if (cnt>0) sticky /= cnt;
		return sticky;
	}

	private static ArrayList<Integer> getRegions(int ind, int length, int numreg) {
		int chunk = length / numreg;
		int rem = length % numreg;
		ArrayList<Integer> res = new ArrayList<Integer>();
		for (int pos = 0; pos < numreg; pos++) {
			if(ind>=pos*chunk && ind < (pos+1)*chunk + rem)
				res.add(pos);
		}
		return res;
	}

	private int[] getShuffledIndices() {
		int[] pos = new int[nreg];
		for (int i=0; i<nreg; i++) pos[i] = i;
		for (int i=0; i<10; i++) {
			int a = r.nextInt(nreg);
			int b = r.nextInt(nreg);
			int sub = pos[a];
			pos[a] = pos[b];
			pos[b] = sub;
		}
		return pos;
	}

	private static int editDist(String seq1, String seq2) {
		int d = 0;
		assert(seq1.length() == seq2.length());
		for (int i=0; i < seq1.length(); i++) {
			if(seq1.charAt(i)!='_' && seq2.charAt(i)!='_' 
					&& seq1.charAt(i) != seq2.charAt(i)) d++;
		}
		return d;
	}

	private String insert(Motif m, String seq, int pos) {
		assert(fits(m,seq,pos)); 
		String res = "";
		for (int i=0; i < seq.length(); i++) {
			if (i<pos || i >= pos+m.regex.length()) res += seq.charAt(i);
			else if (seq.charAt(i) == '_') res += m.regex.charAt(i-pos);
			else if (m.regex.charAt(i-pos) == '_') res += seq.charAt(i);
			else {
				assert(seq.charAt(i) == m.regex.charAt(i-pos));
				res += seq.charAt(i);
			}
		}
		return res;
	}

	private ArrayList<Integer> getPossiblePositions(Motif m, String seq, int region) {
		ArrayList<Integer> resArrayList = new ArrayList<Integer>();
		double half = m.regex.length()/2.0;
		for (int st = 0; st < seq.length()-m.regex.length()+1; st++) {
			if (Math.floor((st+half)*nreg*1.0/seq.length()) == region && fits(m, seq, st)) 
				resArrayList.add(st);
		}
		return resArrayList;
	}

	public boolean fits(Motif m, String seq, int position) {
		for (int i=0; i< m.regex.length(); i++) {
			if (m.regex.charAt(i)!='_' && seq.charAt(i+position) != '_' && 
				m.regex.charAt(i) != seq.charAt(i+position)) {
				return false;
			}
		}
		return true;
	}
	
	private static boolean matchInReg(Motif m, String seq, int reg, int numreg) {
		assert(!seq.contains("_"));
		String regex = m.regex.replace("_", "[ACTG]{1}");
		Pattern p = Pattern.compile(regex);
		Matcher matcher = p.matcher(seq);
		while(matcher.find()) {
			double middle = (matcher.start() + matcher.end())/2.0;
			if(Math.floor(middle*numreg*1.0/seq.length()) == reg)
				return true;
		}
		return false;
	}

	// STATS #######################################################
	
	private void plotNeighborStats(ArrayList<String> generated, int maxdist) {
		TreeMap<Double,Integer> stats_close = new TreeMap<Double, Integer>();
		for (int i =0; i < seqs.length; i++) stats_close.put(i*1.0, 0);
		for (String seq: generated) {
			countNeighboringSequences(stats_close, seq, maxdist);
		}
		int tot = 0, b = 0, d = 0;
		for (Double ind: stats_close.keySet()) {
			if("bright".equals(classes[ind.intValue()])) b+= stats_close.get(ind);
			else d+=stats_close.get(ind);
			tot+=stats_close.get(ind);
		}
		Double perc = b*100.0/tot;
		System.out.print( "bright: " + perc.intValue() + "%\n");
		Common.plotHist(stats_close);
	}
	
	private void countNeighboringSequences(
			TreeMap<Double, Integer> stats_close, String seq, int radius) {
		for (int i = 0 ; i < seqs.length; i++) {
			int d = editDist(seqs[i],seq);
			if (d<=radius) {
				stats_close.put(i*1.0, stats_close.get(i*1.0) + 1);
			}
		}
	}
	
	private boolean isClosestKnownWithin(int minD, int maxD, String s) {
		int closestD = s.length(), d = 0;
		assert(s.length() == seqs[0].length());
		for (int i = 0 ; i < seqs.length; i++) {
			d = editDist(seqs[i],s);
			if( d < closestD)
				closestD = d;
		}
		return (minD<=closestD) && (maxD >= closestD);
	}
	
	// DISTRIBUTIONS ###################################################
	
	private void computeDists(String clas) {
		// compute the total II in bright
		double sumII = getTotalIntensity(clas);
		for (int pos = 0; pos < nreg; pos++) {
			distBnum.add(computeDist(clas,sumII,pos,bcnt));
			distDnum.add(computeDist(clas,sumII,pos,dcnt));
			distSticky.add(computeDist(clas, sumII, pos, bsticky));
		}
		distMotif = computeDistMotif(clas);
	}

	private ArrayList<TreeMap<Double,Motif>> computeDistMotif(String clas) {
		ArrayList<TreeMap<Double,Motif>> res = new ArrayList<TreeMap<Double,Motif>>();
		double[] totalII = new double[nreg];
		for (int reg = 0; reg < nreg; reg++) {
			for (Motif m: motifs ) {
				for (int i = 0; i < n; i++) {
					if (classes[i].equals(clas)) {
						if (matchInReg(m,seqs[i],reg,nreg)){
							m.prob[reg] += intensity[i];
							totalII[reg] += intensity[i];
						}
					}
				}
			}
			if (totalII[reg]>0)
				for (Motif m: motifs) 
					m.prob[reg] /= totalII[reg];
			TreeMap<Double, Motif> prob2mot = new TreeMap<Double, Motif>();
			double acc = 0.0;
			for (Motif m: motifs) {
				prob2mot.put(m.prob[reg]+acc, m);
				acc += m.prob[reg];
			}
			res.add(prob2mot);
		}
		return res;
	}

	public double getTotalIntensity(String clas) {
		double sumII = 0.0;
		for (int i = 0; i < n; i++) {
			if (classes[i].equals(clas)) {
				sumII += intensity[i];
			}
		}
		return sumII;
	}

	public TreeMap<Double, Integer> computeDist(String clas, double sumII, int pos, int[][] vals) {
		TreeMap<Integer, Double> m = new TreeMap<Integer, Double>();	
		for (int i = 0; i < n; i++) {
			if (classes[i].equals(clas)) {
				if (!m.containsKey(vals[i][pos])) m.put(vals[i][pos], 0.0);
				m.put(vals[i][pos], m.get(vals[i][pos]) + intensity[i]/sumII);
			}
		}
		TreeMap<Double, Integer> prob2num = new TreeMap<Double, Integer>();
		double acc = 0.0;
		for (Integer num: m.keySet()) {
			prob2num.put(m.get(num)+acc, num);
			acc += m.get(num);
		}
		return prob2num;
	}
	
	public TreeMap<Double, Double> computeDist(String clas, double sumII, int pos, double[][] vals) {
		TreeMap<Double, Double> m = new TreeMap<Double, Double>();	
		for (int i = 0; i < n; i++) {
			if (classes[i].equals(clas)) {
				if (!m.containsKey(vals[i][pos])) m.put(vals[i][pos], 0.0);
				m.put(vals[i][pos], m.get(vals[i][pos]) + intensity[i]/sumII);
			}
		}
		TreeMap<Double, Double> prob2dbl = new TreeMap<Double, Double>();
		double acc = 0.0;
		for (Double num: m.keySet()) {
			prob2dbl.put(m.get(num)+acc, num);
			acc += m.get(num);
		}
		return prob2dbl;
	}
	
	// READ/WRITE ###################################################
	
	private void readSequences(String fn) throws IOException {
		n= getNumLines(fn) -1;
		BufferedReader br = new BufferedReader(new FileReader(fn));
		String line = br.readLine();
		String[] sline = line.split(",");
		int colFeat = 5, colSeq = 1, colII = 4, colBPos = -1, colDPos = -1, 
				colBsize = -1, colBSticky = -1, colClass = -1;
		for (int i = 0; i < sline.length; i++) {
			if (sline[i].equalsIgnoreCase("BrightPos0")) colBPos = i;
			if (sline[i].equalsIgnoreCase("DarkPos0")) colDPos = i;
			if (sline[i].equalsIgnoreCase("BaseSize0")) colBsize = i;
			if (sline[i].equalsIgnoreCase("BaseSticky0")) colBSticky = i;
			if (sline[i].equalsIgnoreCase("class")) colClass = i;
		}
		nreg=colDPos-colBPos;
		motifs = new ArrayList<Motif>();
		for (int i = colFeat; i < colBPos; i++) {
			motifs.add(new Motif(sline[i],nreg)); 
			if (sline[i].contains("_")) 
				motifs.add(new Motif(sline[i].replace("_", ""),nreg)); 
		}
		
		seqs = new String[n];
		bcnt = new int[n][nreg];
		dcnt = new int[n][nreg];
		bsize = new double[n][nreg];
		bsticky = new double[n][nreg];
		classes = new String[n];
		intensity = new double[n];
		
		int nr = 0;
		while((line = br.readLine())!=null) {
			sline = line.split(",");
			seqs[nr] = sline[colSeq];
			classes[nr] = sline[colClass];
			intensity[nr] = Double.parseDouble(sline[colII]);
			for(int i = 0; i < nreg; i++) {
				bcnt[nr][i] = Integer.parseInt(sline[colBPos + i]);
				dcnt[nr][i] = Integer.parseInt(sline[colDPos + i]);
				bsize[nr][i] = Double.parseDouble(sline[colBsize + i]);
				bsticky[nr][i] = Double.parseDouble(sline[colBSticky + i]);
			}	
			nr++;
		}
	}
	
	private int getNumLines(String fn) throws FileNotFoundException, IOException {
		LineNumberReader lnr = new LineNumberReader(new FileReader(fn));
		lnr.skip(Long.MAX_VALUE);
		return lnr.getLineNumber();
	}
	
	private void writeSequences(String fn, ArrayList<String> sequences, int MinD, int MaxD) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(fn));
		int id =0;
		TreeSet<String> nseq = new TreeSet<String>();
		for(String s: sequences) nseq.add(s);
		bw.write("Name,Sequence\n");
		for(String s:nseq) {
			if (isClosestKnownWithin(MinD,MaxD,s)) {
				bw.write("10basegen"+id + "," + s + "\n");
				id++;
			}
		}
		System.out.print("Wrote " + id + " new sequences\n");
		bw.close();
	}
	
	// TEST #####################################

	public static void testMatchInRegion() throws IOException {
		String seq = "ATCCCCAAAC";
		Motif m = new Motif("AA_C+",3);
		assert(matchInReg(m, seq, 2, 3));
	}
	
	// MAIN #####################################
	public static void main(String[] args) throws IOException{
//		testMatchInRegion();
//		System.exit(0);
		String inSeqFn = args[0];
		Integer numNew = Integer.parseInt(args[1]);
		int Tmin = Integer.parseInt(args[2]);
		int Tmax = Integer.parseInt(args[3]);
		
		Generator gen = new Generator(inSeqFn);
		ArrayList<String> newseqs = gen.generate(numNew);
		gen.plotNeighborStats(newseqs, 3);
		gen.writeSequences(inSeqFn.substring(0,inSeqFn.indexOf(".")) + "-new.csv", 
				newseqs, Tmin, Tmax);
	}

}
