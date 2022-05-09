package main;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

public class InParam {
	
	private static InParam ARGS = new InParam();
	private List<String> input_sam = new ArrayList<>();
	private List<String> treat_files = new ArrayList<>();
	private List<String> control_files = new ArrayList<>();
	private List<String> treatFiles2 = new ArrayList<>();
	private List<String> controlFiles2 = new ArrayList<>();
	private String genome_file = null;
	private String exon_file = null;
	private String peak_file = null;
	private String mut_file = null;
	private String circ_file = null;
	private String gene_file = null;
	private String scaf_file = null;
	private String base_quality = "J";
	private String combineP = "log";
	private Set<Character> bases = null;
	private double read_depth = -1.0;
	private double total_reads = 10E6;
	private double enrich_num = -1.0;
	private double back_scripts = 0.5;
	private double fix_allele_fre = -1.0;
	private int method = -1;
	private int exam = Integer.MAX_VALUE;
	private int window_size = 50;
	private int read_length = 300;
	private int alignment_length = 150;
	private int min_align_len = 30;
	private int rand_read_len = 100;
	private int min_back_read = 0;
	private int permutate = 1000;
	private long random_seed = 20210625;
	private double mutation_prop = 0.5;
	private double background_num = 0.01;
	private double linear_prop = 50.0;
	private String out_pre = null;
	private boolean sam_format = false;
	private boolean chr_sep = false;
	private boolean out_zip = true;
	private boolean sim_with_peak = false;
	private boolean linear_peak = true;
	private boolean mix_peak = false;
	private boolean rand_depth = false;
	private boolean annote = false;
	private boolean rtMut = false;
	private boolean replicate = false;
	private boolean negative = false;
	private boolean help = false;
	
	private InParam() {
	}
	
	public List<String> getInputSams() {
		return input_sam;
	}

	public String getTreatFile() {
		return treat_files.isEmpty() ? null : treat_files.get(0);
	}
	
	public List<String> getTreatFiles() {
		return treat_files;
	}
	
	public List<String> getTreatFiles2() {
		return treatFiles2;
	}
	
	public String getControlFile() {
		return control_files.isEmpty() ? null : control_files.get(0);
	}
	
	public List<String> getControlFiles() {
		return control_files;
	}
	
	public List<String> getControlFiles2() {
		return controlFiles2;
	}
	
	public String getGenomeFile() {
		return genome_file;
	}

	public String getExonFile() {
		return exon_file;
	}

	public String getGeneFile() {
		return gene_file;
	}
	
	public String getPeakFile() {
		return peak_file;
	}

	public String getMutFile() {
		return mut_file;
	}

	public String getCircFile() {
		return circ_file;
	}
	
	public String getScafFile() {
		return scaf_file;
	}

	public String getQuality() {
		return base_quality;
	}

	public Set<Character> getBases() {
		return bases;
	}

	public boolean isVaildBase(char c) {
		return bases == null || c == 0 || bases.contains(c);
	}
	
	public String getCombineP() {
		return combineP;
	}
	
	public double getReadDepth() {
		return read_depth;
	}
	
	public double getTotalReads() {
		return total_reads;
	}

	public double getEnrich() {
		return enrich_num;
	}

	public double getNonPeakBack() {
		return background_num;
	}
	
	public double getLinearProp() {
		return linear_prop;
	}
	
	public double getFixAlleleFre() {
		return fix_allele_fre;
	}
	
	public double getBackScripts() {
		return back_scripts;
	}
	
	public int getReadLength() {
		return read_length;
	}

	public int getPermutate() {
		return permutate;
	}
	
	public int getWindowSize() {
		return window_size;
	}
	
	public int getMethod() {
		return method;
	}
	
	public int getAlignmentLength() {
		return alignment_length;
	}

	public int getMinAlignmentLength() {
		return min_align_len;
	}

	public int getMinBackRead() {
		return min_back_read;
	}
	
	public int getRandReadLength() {
		return rand_read_len;
	}

	public double getMutationProportion() {
		return mutation_prop;
	}

	public int getExam() {
		return exam;
	}

	public String getOutPrefix() {
		return out_pre;
	}
	
	public boolean isSamFormat() {
		return sam_format;
	}
	
	public boolean isChrSep() {
		return chr_sep;
	}
	
	public boolean isOutZip() {
		return out_zip;
	}
	
	public boolean isPeakSim() {
		return sim_with_peak;
	}
	
	public boolean isLinearPeak() {
		return linear_peak;
	}
	
	public boolean isSharedPeak() {
		return mix_peak;
	}
	
	public boolean isRandDepth() {
		return rand_depth;
	}
	
	public boolean isAnnote() {
		return annote;
	}
	
	public boolean isRtMut() {
		return rtMut;
	}
	
	public boolean isReplicate() {
		return replicate;
	}
	
	public boolean isNegative() {
		return negative;
	}
	
	public boolean isPairEnd() {
		return read_length != alignment_length;
	}
	
	public boolean setParams(String[] args) {
		return set(args);
	}
	
	private int setMethod(String m_string) {
		List<String> methods = new ArrayList<>();
		methods.add("sim");
		methods.add("maf");
		methods.add("seq");
		methods.add("stat");
		methods.add("cover");
		methods.add("ribosnitch");
		methods.add("split");
		methods.add("twosample");
		return methods.indexOf(m_string);
	}
	
	private boolean set(String[] args) {
		if (args.length <= 1) {
			help = true;
		}
		if (args.length > 1) {
			method = setMethod(args[0].toLowerCase());
		}
		for (int i = 1; i < args.length; i++) {
			switch (args[i].toLowerCase()) {
			case "-is":
				while(++i < args.length && args[i].charAt(0) != '-') {
					input_sam.add(args[i]);
				}
				--i;
				break;
			
			case "-o":
				if (++i < args.length) {
					out_pre = args[i];
				}
				break;
				
			case "-tf":
				while(++i < args.length && args[i].charAt(0) != '-') {
					treat_files.add(args[i]);
				}
				--i;
				break;
				
			case "-tf2":
				while(++i < args.length && args[i].charAt(0) != '-') {
					treatFiles2.add(args[i]);
				}
				--i;
				break;
				
			case "-cf":
				while(++i < args.length && args[i].charAt(0) != '-') {
					control_files.add(args[i]);
				}
				--i;
				break;
				
			case "-cf2":
				while(++i < args.length && args[i].charAt(0) != '-') {
					controlFiles2.add(args[i]);
				}
				--i;
				break;
				
			case "-gf":
				if (++i < args.length) {
					genome_file = args[i];
				}
				break;
				
			case "-ef":
				if (++i < args.length) {
					exon_file = args[i];
				}
				break;
				
			case "-bases":
				bases = new HashSet<>();
				if (++i < args.length && args[i].charAt(0) != '-') {
					for (int j = 0; j < args[i].length(); ++j) {
						char c = Character.toUpperCase(args[i].charAt(j));
						if (c >= 'A' && c <= 'Z') {
							bases.add(c);
						}
					}
				}
				else {
					--i;
				}

				break;
				
			case "-genef":
				if (++i < args.length) {
					gene_file = args[i];
				}
				break;
				
			case "-mutf":
				if (++i < args.length) {
					mut_file = args[i];
				}
				break;
				
			case "-peakf":
				if (++i < args.length) {
					peak_file = args[i];
				}
				break;
				
			case "-circf":
				if (++i < args.length) {
					circ_file = args[i];
				}
				break;
				
			case "-sf":
				if (++i < args.length) {
					scaf_file = args[i];
				}
				break;
				
			case "-baseq":
				if (++i < args.length) {
					base_quality = args[i];
				}
				break;
				
			case "-depth":
				if (++i < args.length) {
					read_depth = Math.abs(Double.parseDouble(args[i]));
				}
				break;
				
			case "-reads":
				if (++i < args.length) {
					total_reads = Math.abs(Double.parseDouble(args[i]));
				}
				break;
				
			case "-enrich":
				if (++i < args.length) {
					enrich_num = Math.abs(Double.parseDouble(args[i]));
				}
				break;
				
			case "-npback":
				if (++i < args.length) {
					background_num = Math.abs(Double.parseDouble(args[i]));
				}
				break;
				
			case "-lnp":
				if (++i < args.length) {
					linear_prop = Math.abs(Double.parseDouble(args[i]));
				}
				break;
				
			case "-faf":
				if (++i < args.length) {
					fix_allele_fre = Math.abs(Double.parseDouble(args[i]));
				}
				break;
				
			case "-exam":
				if (++i < args.length) {
					exam = Math.abs(Integer.parseInt(args[i]));
				}
				break;
				
			case "-rlen":
				if (++i < args.length) {
					read_length = Math.abs(Integer.parseInt(args[i]));
				}
				break;
				
			case "-win":
				if (++i < args.length) {
					window_size = Math.abs(Integer.parseInt(args[i]));
				}
				break;
				
			case "-alen":
				if (++i < args.length) {
					alignment_length = Math.abs(Integer.parseInt(args[i]));
				}
				break;
				
			case "-permu":
				if (++i < args.length) {
					permutate = Math.abs(Integer.parseInt(args[i]));
				}
				break;
				
			case "-minalen":
				if (++i < args.length) {
					min_align_len = Math.abs(Integer.parseInt(args[i]));
				}
				break;
				
			case "-minbr":
				if (++i < args.length) {
					min_back_read = Math.abs(Integer.parseInt(args[i]));
				}
				break;
				
			case "-randrlen":
				if (++i < args.length) {
					rand_read_len = Math.abs(Integer.parseInt(args[i]));
				}
				break;
				
			case "-mutp":
				if (++i < args.length) {
					mutation_prop = Math.abs(Double.parseDouble(args[i]));
				}
				break;
				
			case "-seed":
				if (++i < args.length && args[i].charAt(0) != '-') {
					random_seed = Long.parseLong(args[i]);
				}
				else {
					--i;
				}
				CommonMethod.rand = new Random(random_seed);
				break;
				
			case "-sam":
				sam_format = true;
				break;
				
			case "-peak":
				sim_with_peak = true;
				break;
				
			case "-chrs":
				chr_sep = true;
				break;
				
			case "-nolinear":
				linear_peak = false;
				break;
				
			case "-share":
				mix_peak = true;
				break;
				
			case "-anno":
				annote = true;
				break;
				
			case "-mut":
				rtMut = true;
				break;
				
			case "-neg":
				negative = true;
				break;
				
			case "-repli":
				replicate = true;
				break;
				
			case "-h":
				help = true;
				break;
				
			default:
				System.err.println("Warning: Unknown param " + args[i]);
				break;
			}
		}
		return checkParams();
	}

	public static InParam getParams() {
		return ARGS;
	}
	
	private boolean checkParams() {
		try {
			if (help) {
				help();
				return false;
			}
			if (method < 0) {
				helpMethod();
			}

//			if (genome_file == null || exon_file == null) {
//				throw new FileNotFoundException("Lack parameter(s)");
//			}
//			FileReader fr = new FileReader(new File(genome_file));
//			fr.close();
//			fr = new FileReader(new File(exon_file));
//			fr.close();
//			
//			if (mutation_prop > 1.0) {
//				throw new NumberFormatException("Mutation proportion should between 0.0 and 1.0");
//			}
//			if (read_length <= alignment_length) {
//				read_length = alignment_length;
//				System.out.println("Warning: Segment length is no greater than read length, simulating single-end fastq using read length");
//			}
		} catch (Exception e) {
			System.err.println(e.getClass().getSimpleName() + ": " + e.getMessage());
			help();
			return false;
		}
		return true;
	}

	private void helpMethod() {
		System.out.println("To Use:\njava [-Xmx16g] -jar RTstop.jar <method> [options]\nChoose from a method: sim, ribosnitch, twosample");
	}
	
	private void help() {
		InputStream is = InParam.class.getResourceAsStream("/Help");
		try {
			InputStreamReader isr = new InputStreamReader(is, "UTF-8");
			BufferedReader reader = new BufferedReader(isr);
			String line = null;
			while ((line = reader.readLine()) != null) {
				System.out.println(line);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void setEnrich(double enrich_num) {
		this.enrich_num = enrich_num;
	}
}
