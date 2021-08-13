package sim;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.lang.reflect.AccessibleObject;
import java.lang.reflect.Array;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.distribution.*;

import genome.*;
import bed.*;
import exonseq.SeqMethod;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import main.BamMethod;
import main.CommonMethod;
import main.InParam;
import main.Method;
import sim.Scaffolds.Scaffold;

public class SimMethod {
	
	private final static NormalDistribution mut_per_dis = new NormalDistribution(0.0, 0.15);
	private final static NormalDistribution OR_dis = new NormalDistribution(0.0, 0.9);
	
	private final static char[] base_units = new char[] {'A', 'G', 'C', 'T'};
	private final static HashMap<Character, Character> BASE_MAP = new HashMap<>();
	static {
		BASE_MAP.put('A', 'T');
		BASE_MAP.put('T', 'A');
		BASE_MAP.put('G', 'C');
		BASE_MAP.put('C', 'G');
		BASE_MAP.put('t', 'a');
		BASE_MAP.put('a', 't');
		BASE_MAP.put('g', 'c');
		BASE_MAP.put('c', 'g');
	}
	private static int id_count = 0;
	private int base_num = 0;
	private double read_num = 0.0;
	private double read_num_ip = 0.0;
	private double linear_ratio = 0.0;
	private String FixMapQ = null;
	private Map<String, IntervalTree<Mutation>> muts = new HashMap<>();
	private Map<String, IntervalTree<Peak>> peak_table = new HashMap<>();
	private Map<String, IntervalTree<CircleClip>> circ_table = new HashMap<>();
	private Map<String, IntervalTree<Scaffold[]>> scaf_table = new HashMap<>();
	private Genome genome = null;
	private NormalDistribution enrich_dis = new NormalDistribution(0.5, 0.3);
//	private int[] counts = new int[149];
	private InParam args = null;
	
	private static SimMethod m = new SimMethod();
	
	
	private SimMethod() {
	}
	
	public static void run() {
		m.args = InParam.getParams();
		try {
			m.genome = Method.loadGenomeInfo(m.args.getExonFile(), m.args.getGenomeFile());
			Method.printNow("Genome complete!");
			m.loadInfo();
			Method.printNow("Load complete!");
			m.scaf_table = m.getAlleleScaffold(Method.readScaffolds(m.args.getScafFile(), null));
			m.ensureRelation(null);
			Method.printNow("SNV complete!");
			m.FixMapQ = buildQualiString(m.args.getAlignmentLength() + 16, "J");
			ensureMutDes(m.peak_table);
			
			m.prepareSim();
			m.sim();
			m.simPeak();
			
			if (m.args.getMutFile() != null) {
				Method.getInstance().calMAF(m.muts, m.genome);
				Method.writeFile(m.args.getOutPrefix() + "_mut.bed", m.muts, true, false, Mutation.getHeader());
			}
			if (m.args.getCircFile() != null) {
				Method.writeFile(m.args.getOutPrefix() + "_circ.bed", m.circ_table, true, false, CircleClip.getHeader());
			}
			if (m.args.getPeakFile() != null) {
				Method.writeFile(m.args.getOutPrefix() + "_peak.bed", m.peak_table, true, false, Peak.getHeader());
			}
			if (m.args.getScafFile() != null) {
				m.writeMutExonSeq();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
//	public static void runScaf() {
//		m.args = InParam.getParams();
//		try {
//			m.genome = Method.loadGenomeInfo(m.args.getExonFile(), m.args.getGenomeFile());
//			Method.printNow("Genome complete!");
//			m.muts = Method.readMutation(m.args.getMutFile(), m.muts);
//			
//			List<String> chrRemove = new ArrayList<>();
////			String geneName = "gene-MTMR12";
//			m.genome.getChrMap().forEach((chrString, chr) -> {
//				if (!m.muts.containsKey(chrString)) {
//					chrRemove.add(chrString);
//					return;
//				}
//				for (Node<Gene> geneNode : chr.getGeneTree()) {
////					if (!geneName.equals(geneNode.getValue().getGene_id())) {
//					if (!m.muts.get(chrString).overlappers(geneNode.getStart(), geneNode.getEnd()).hasNext()) {
//						chr.getGeneTree().remove(geneNode.getStart(), geneNode.getEnd());
//					}
//				}
//			});
//			for (String chr : chrRemove) {
//				m.genome.getChrMap().remove(chr);
//			}
//			
//
//			m.scaf_table = m.getAlleleScaffold(Method.readScaffolds(m.args.getScafFile(), null));
//			m.ensureRelation(null);
//			Method.printNow("SNV complete!");
//			m.writeMutExonSeq();
//			m.FixMapQ = buildQualiString(m.args.getAlignmentLength() + 16, "J");
//			
//			m.prepareSim();
//			m.simScaf();
//			
//			if (m.args.getMutFile() != null) {
//				Method.getInstance().calMAF(m.muts, m.genome);
//				Method.writeFile(m.args.getOutPrefix() + "_mut.bed", m.muts, true, false, Mutation.getHeader());
//			}
//			if (m.args.getCircFile() != null) {
//				Method.writeFile(m.args.getOutPrefix() + "_circ.bed", m.circ_table, true, false, CircleClip.getHeader());
//			}
//			if (m.args.getPeakFile() != null) {
//				Method.writeFile(m.args.getOutPrefix() + "_peak.bed", m.peak_table, true, false, Peak.getHeader());
//			}
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//	}

	private void loadInfo() throws IOException {
		boolean peak_flag = args.getPeakFile() != null && Files.exists(Paths.get(args.getPeakFile()));
		boolean mut_flag = args.getMutFile() != null && Files.exists(Paths.get(args.getMutFile()));
		boolean circ_flag = args.getCircFile() != null && Files.exists(Paths.get(args.getCircFile()));
		peak_table = peak_flag ? Method.getInstance().readPeak(args.getPeakFile(), peak_table) : peak_table;
		Method.printNow("Peak complete!");
		muts = mut_flag ? Method.readMutation(args.getMutFile(), muts) : muts;
		Method.printNow("SNV complete!");
		circ_table = circ_flag ? Method.getInstance().readCircFile(args.getCircFile(), circ_table, genome) : circ_table;
		Method.printNow("Circ complete!");
//		List<List<Gene>> gene_enrich = Method.getInstance().readGeneFile(args.getGeneFile());
		ensureRelation(null);
	}

//	private HashMap<String, IntervalTree<Peak>> simPeak(HashMap<String, IntervalTree<Peak>> peak_table) throws IOException {
//		final HashMap<String, IntervalTree<Peak>> peak_tree = peak_table == null ? new HashMap<>() : peak_table;
//		if (args.getPeakFile() == null) {
//			return peak_tree;
//		}
//		String[] des = {"Peak", "Gene", "None", "Exon"};
//		genome.getChrMap().forEach((chr_string, chr) -> {
//			UniformIntegerDistribution uid = new UniformIntegerDistribution(1, chr.getSeq().length());
//			int count = chr.getSeq().length() / 50000;
//			int none_count = count / 50;
//			if (!peak_tree.containsKey(chr_string)) {
//				peak_tree.put(chr_string, new IntervalTree<>());
//			}
//			while (count > 0) {
//				int s = uid.sample();
//				Iterator<Node<Gene>> nodes = chr.getGeneTree().overlappers(s, s);
//				if (nodes.hasNext()) {
//					if (!peak_tree.get(chr_string).overlappers(s - args.getAlignmentLength(),
//							s + args.getAlignmentLength()).hasNext() && --count >= 0) {
//						peak_tree.get(chr_string).put(s, s, new Peak(s, chr.getExonTree().overlappers(s, s).hasNext() ? des[3] : des[1]));
//					}
//				}
//				else {
//					if (!peak_tree.get(chr_string).overlappers(s - args.getAlignmentLength(),
//							s + args.getAlignmentLength()).hasNext() && --none_count >= 0) {
//						peak_tree.get(chr_string).put(s, s, new Peak(s, des[2]));
//					}
//				}
//			}
//		});
//		Method.writeFile(args.getPeakFile(), peak_tree, true, true, Peak.getHeader());
//		return peak_tree;
//	}
	
//	private HashMap<String, IntervalTree<Mutation>> simMut(HashMap<String, IntervalTree<Mutation>> mut_table,
//			HashMap<String, IntervalTree<Peak>> peak_table) {
//		final HashMap<String, IntervalTree<Mutation>> mut_tree = mut_table == null ? new HashMap<>() : mut_table;
//		if (args.getMutFile() == null) {
//			return mut_tree;
//		}
//		final char[] base_units = new char[] {'A', 'G', 'C', 'T'};
//		UniformIntegerDistribution buid = new UniformIntegerDistribution(0, base_units.length - 1);
//		genome.getChrMap().forEach((chr_string, chr) -> {
//			if (!mut_tree.containsKey(chr_string)) {
//				mut_tree.put(chr_string, new IntervalTree<>());
//			}
//			if (peak_table.containsKey(chr_string)) {
//				for (Node<Peak> peak_node : peak_table.get(chr_string)) {
//					for (int pos : randomNoRepeat(peak_node.getStart() - args.getReadLength(), 
//							peak_node.getStart() + args.getReadLength(), 3)) {
//						char c = base_units[buid.sample()];
//						while (c == chr.getSeq().charAt(pos - 1)) {
//							c = base_units[buid.sample()];
//						}
//						mut_tree.get(chr_string).put(pos, pos, new Mutation(pos, chr.getSeq().charAt(pos - 1), c, "Peak"));
//					}
//				}
//			}
//			int mut_count = chr.getSeq().length() / 50000;
//			UniformIntegerDistribution uid = new UniformIntegerDistribution(1, chr.getSeq().length());
//			while (--mut_count >= 0) {
//				int s = uid.sample();
//				if (!mut_tree.get(chr_string).overlappers(s - args.getReadLength(), s + args.getReadLength()).hasNext()) {
//					char c = base_units[buid.sample()];
//					while (c == chr.getSeq().charAt(s - 1)) {
//						c = base_units[buid.sample()];
//					}
//					mut_tree.get(chr_string).put(s, s, new Mutation(s, chr.getSeq().charAt(s - 1), c, 
//							chr.getGeneTree().overlappers(s, s).hasNext() ? "Gene" : "None"));
//				}
//				
//			}
//		});
//		return mut_tree;
//	}
	
//	private HashMap<String, IntervalTree<CircleClip>> simCirc(HashMap<String, IntervalTree<CircleClip>> circ_table) {
//		final HashMap<String, IntervalTree<CircleClip>> circ_tree = circ_table == null ? new HashMap<>() : circ_table;
//		if (args.getCircFile() == null) {
//			return circ_tree;
//		}
//		genome.getChrMap().forEach((chr_string, chr) -> {
//			if (!circ_tree.containsKey(chr_string)) {
//				circ_tree.put(chr_string, new IntervalTree<>());
//			}
//			for (int i = chr.getScrpitTree().size() / 100 + 1; i > 0; i--) {
//				CircleClip circle = simCircPeak(chr, true);
//				if (circle != null) {
//					circ_tree.get(chr_string).put(circle.getStart() + 1, circle.getEnd(), circle);
//				}
//			}
//		});
//		return circ_tree;
//	}

	
//	private CircleClip simCircPeak(Chromosome chr, boolean in_exon) {
//		Node<Exon> node = chr.getExonTree().findByIndex((int) (Math.random() * chr.getExonTree().size()));
//		while (node.getValue().getScript().isCirc()) {
//			node = chr.getExonTree().findByIndex((int) (Math.random() * chr.getExonTree().size()));
//		}
//		return simCircExon(node.getValue(), chr.getSeq());
//	}
	
//	private CircleClip simCircExon(Exon exon, String bases){
//		Transcript script = exon.getScript();
//		int index = 0;
//		for (Exon start_exon = script.getExons(); start_exon != null; start_exon = start_exon.getDownStream()) {
//			++index;
//		}
//		int left_index = (int) (Math.random() * (index / 2 + 1));
//		int right_index = (int) (Math.random() * (index - left_index)) + left_index;
//		int start = 0;
//		int end = 0;
//		int len = 0;
//		for (Exon start_exon = script.getExons(); start_exon != null && right_index >= 0; start_exon = start_exon.getDownStream()) {
//			start = left_index == 0 ? start_exon.getStart() : start;
//			end = right_index == 0 ? start_exon.getEnd() : end;
//			len += left_index <= 0 && right_index >= 0 ? start_exon.getLength() : 0;
//			--left_index;
//			--right_index;
//		}
//		if (len <= args.getAlignmentLength()) {
//			return null;
//		}
////		script.setCirc(true);
//		return new CircleClip(start, end, script, isGTAG(bases.substring(end, end + 2), bases.substring(start - 2, start)));
//	}
	
	private boolean isGTAG(String gt, String ag) {
		String[] AG = {"AG", "AC", "AG", "AC", "AT", "GC"};
		String[] GT = {"GT", "AT", "GC", "CT", "GT", "CT"};
		for (int i = 0; i < AG.length; i++) {
			if (AG[i].equals(ag) && GT[i].equals(gt)) {
				return true;
			}
		}
		return false;
	}
	
	private Map<String, IntervalTree<Scaffold[]>> getAlleleScaffold(Map<String, IntervalTree<Scaffolds>> scafTable) {
		Map<String, IntervalTree<Scaffold[]>> out = new HashMap<>();
 		scafTable.forEach((chr, scafTree) -> {
			out.put(chr, new IntervalTree<>());
			for (Node<Scaffolds> scafNode : scafTree) {
//				out.get(chr).put(scafNode.getStart(), scafNode.getEnd(), scafNode.getValue().getRandScafPair());
				out.get(chr).put(scafNode.getStart(), scafNode.getEnd(), scafNode.getValue().getDiffScafPair());
			}
		});
 		return out;
	}
	
	private void ensureRelation(List<List<Gene>> gene_enrich) {
		final int[] peak_num = {0, 0};
		final boolean[] sim_flag = {args.getPeakFile() != null && !Files.exists(Paths.get(args.getPeakFile())),
				args.getCircFile() != null && !Files.exists(Paths.get(args.getCircFile())),
				args.getMutFile() != null && !Files.exists(Paths.get(args.getMutFile()))};
		genome.getChrMap().forEach((chr_string, chr) -> {
			if (sim_flag[0]) {
				peak_table.put(chr_string, new IntervalTree<>());
			}
			IntervalTree<Peak> peak_tree = peak_table.getOrDefault(chr_string, new IntervalTree<>());
			IntervalTree<CircleClip> circ_tree = circ_table.getOrDefault(chr_string, new IntervalTree<>());
			IntervalTree<Exon> exon_tree = chr.getExonTree();
			for (Node<Gene> gene_node : chr.getGeneTree()) {
				for (Gene gene = gene_node.getValue(); gene != null; gene = gene.getAnotherGene()) {
					gene.setSimMaf(0.5 + (args.getFixAlleleFre() < 0.0 ? 
							sampleCutInRange(mut_per_dis, -0.5, 0.5) : 0.0));
				}
			}
			for (Node<CircleClip> circ_node : circ_tree) {
//				if (circ_node.getEnd() == 94171496) {
//					System.out.println("debug");
//				}
				HashSet<Transcript> peak_scripts = new HashSet<>();
				HashSet<Transcript> circ_exon = new HashSet<>();
				int circ_exon_stat = 0;
				HashSet<Transcript> exon_bound = new HashSet<>();
				Iterator<Node<Exon>> exon_nodes = exon_tree.overlappers(circ_node.getStart(), circ_node.getEnd());
				while (exon_nodes.hasNext()) {
					Node<Exon> exon_node = exon_nodes.next();
					circ_exon_stat = 1;
					if (peak_tree.overlappers(exon_node.getStart(), exon_node.getEnd()).hasNext()) {
						for (Exon exon = exon_node.getValue(); exon != null; exon = exon.getAnotherExon()) {
							peak_scripts.add(exon.getScript());
						}
					}
					if (circ_node.getStart() == exon_node.getStart()) {
						for (Exon exon = exon_node.getValue(); exon != null; exon = exon.getAnotherExon()) {
							exon_bound.add(exon.getScript());
						}
					}
					if (circ_node.getEnd() == exon_node.getEnd()) {
						for (Exon exon = exon_node.getValue(); exon != null; exon = exon.getAnotherExon()) {
							circ_exon.add(exon.getScript());
						}
					}
				}
				circ_exon.retainAll(exon_bound);
				if (circ_exon_stat == 0 && chr.getScrpitTree().overlappers(circ_node.getStart(), circ_node.getEnd()).hasNext()) {
					circ_exon_stat = 2;
				}
				if (circ_exon_stat == 0 && chr.getGeneTree().overlappers(circ_node.getStart(), circ_node.getEnd()).hasNext()) {
					circ_exon_stat = 3;
				}
				Transcript script = circ_node.getValue().getScript();
				if (script == null || !circ_exon.contains(script)) {
					HashSet<Transcript> set = new HashSet<>(peak_scripts);
					set.retainAll(circ_exon);
					script = biggestScript(set);
				}
				if (script != null) {
					setCircle(circ_node.getValue(), script, chr);
					script = circ_node.getValue().getScript();
					int len = 0;
					int index = script.getExonLength() / 2;
					Peak peak = null;
					for (Exon exon = script.getExons(); exon != null; exon = exon.getDownStream()) {
						Iterator<Node<Peak>> peak_nodes = peak_tree.overlappers(exon.getStart() + 1, exon.getEnd());
						while (peak_nodes.hasNext()) {
							Node<Peak> peak_node = peak_nodes.next();
							if (Math.abs(index - script.getExonLength() / 2.0) <= Math.abs(len + peak_node.getValue().getStart() - exon.getStart() - script.getExonLength() / 2.0)) {
								index = len + peak_node.getValue().getStart() - exon.getStart();
								peak = peak_node.getValue();
							}
							peak_node.getValue().setCirc(len + peak_node.getValue().getStart() - exon.getStart() <= args.getReadLength()
									|| script.getExonLength() - len - peak_node.getValue().getStart() + exon.getStart() <= args.getReadLength());
						}
						len += exon.getLength();
					}
					circ_node.getValue().setPeak(peak, index);
				}
				else {
					script = biggestScript(circ_exon);
					if (script == null) {
						Gene gene = new Gene(chr, null, null, null, '.', circ_node.getStart() - 1, circ_node.getEnd());
						for (Iterator<Node<Gene>> gene_nodes = chr.getGeneTree().overlappers(circ_node.getStart(), circ_node.getEnd());
								gene_nodes.hasNext();) {
							Node<Gene> gene_node = gene_nodes.next();
							if (gene_node.getStart() <= circ_node.getStart() && gene_node.getEnd() >= circ_node.getEnd()) {
								gene = gene_node.getValue();
								break;
							}
							if (gene.getGene_id() == null) {
								gene = gene_node.getValue();
							}
							else {
								gene = moreNear(circ_node.getValue(), gene, gene_node.getValue());
							}
						}
						if (gene.getGene_id() == null) {
							gene.addAnotherGene(chr.getGeneTree().put(gene.getStart() + 1, gene.getEnd(), gene));
						}
						script = new Transcript(gene, "None", "None", circ_node.getStart() - 1, circ_node.getEnd());
						gene.addScript(script);
						Exon exon = new Exon(script, circ_node.getStart() - 1, circ_node.getEnd());
						script.addLastExon("exon", exon);
						int index = script.getExonLength() / 2;
						Peak peak = null;
						for (Iterator<Node<Peak>> peak_nodes = peak_tree.overlappers(exon.getStart() + 1, exon.getEnd());
								peak_nodes.hasNext();) {
							Node<Peak> peak_node = peak_nodes.next();
							if (Math.abs(index - script.getExonLength() / 2.0) < 
							Math.abs(peak_node.getValue().getStart() - exon.getStart() - script.getExonLength() / 2.0)) {
								index = peak_node.getValue().getStart() - exon.getStart();
								peak = peak_node.getValue();
							}
							peak_node.getValue().setCirc(peak_node.getValue().getStart() - exon.getStart() <= args.getReadLength()
									|| script.getExonLength() - peak_node.getValue().getStart() + exon.getStart() <= args.getReadLength());
						}
						circ_node.getValue().setPeak(peak, index);
						script.setCirc(circ_node.getValue());
						script.setCircTrans(true);
						circ_node.getValue().setExon(false);
						circ_node.getValue().setGTAG(chr != null && isGTAG(chr.getSeq().substring(circ_node.getValue().getEnd(), circ_node.getValue().getEnd() + 2), 
								chr.getSeq().substring(circ_node.getValue().getStart() - 2, circ_node.getValue().getStart())));
						circ_node.getValue().setScript(script);
					}
					else {
						setCircle(circ_node.getValue(), script, chr);
					}
				}
				PoissonDistribution expression_dis = new PoissonDistribution(1.0);
				for (Node<Gene> gene_node : chr.getGeneTree()) {
					for (Transcript s = gene_node.getValue().getScript(); s != null; s = s.getGeneScript()) {
//						double exp = expression_dis.sample();
//						while (exp < 1E-6) {
//							exp = expression_dis.sample();
//						}
//						s.setExpression(exp);
						s.setExpression(Math.exp(expression_dis.sample() - 1.0));
					}
				}
			}
			for (Node<Mutation> mut_node : muts.getOrDefault(chr_string, new IntervalTree<>())) {
				Transcript script = null;
				for (Iterator<Node<Exon>> exon_nodes = chr.getExonTree().overlappers(mut_node.getStart(), mut_node.getStart());
						exon_nodes.hasNext();) {
					Node<Exon> exon_node = exon_nodes.next();
					for (Exon exon = exon_node.getValue(); exon != null; exon = exon.getAnotherExon()) {
						if (script == null || script.getExonLength() < exon.getScript().getExonLength()) {
							script = exon.getScript();
						}
						exon.getScript().setOnceMut(mut_node.getValue());
					}
				}
				if (mut_node.getValue().getScript() == null) {
					mut_node.getValue().setScript(script);
				}
				if (mut_node.getValue().getScript() == null) {
					muts.get(chr_string).remove(mut_node.getStart(), mut_node.getEnd());
				}
				else {
					mut_node.getValue().setMaf(sampleDiscreteMut(mut_node.getValue().getScript().getGene().getSimMaf()));
				}
			}
			
			if (!muts.containsKey(chr_string)) {
				muts.put(chr_string, new IntervalTree<>());
			}
			IntervalTree<Mutation> mut_tree = muts.get(chr_string);
			for (Node<Peak> peak_node : peak_tree) {
				Peak peak = peak_node.getValue();
				HashSet<Transcript> linear_scripts = new HashSet<>();
				HashSet<Gene> used_genes = new HashSet<>();
				for (Iterator<Node<Exon>> exon_nodes = chr.getExonTree().overlappers(peak_node.getEnd(), peak_node.getEnd());
						exon_nodes.hasNext();) {
					Node<Exon> exon_node = exon_nodes.next();
					for (Exon exon = exon_node.getValue(); exon != null; exon = exon.getAnotherExon()) {
						if (exon.getScript().isCirc()) {
							peak_num[0] = exon.getScript().isPeak() ? peak_num[0] : (peak_num[0] + 1);
							Transcript script = uniqPeakScript(exon.getScript(), peak, peak_tree, chr);
							peak.addScript(script);
							if (sim_flag[2] && script != null) {
								mutSim(mut_tree, exon.getScript(), peak);
							}
						}
						else {
							linear_scripts.add(exon.getScript());
						}
						exon.getScript().setPeak(true);
						if (used_genes.add(exon.getScript().getGene())) {
							exon.getScript().getGene().incPeakCount();
						}
					}
				}
				if (args.isSharedPeak()) {
					HashSet<Transcript> used_scripts = new HashSet<>();
					for (Transcript script = peak.getScript(); script != null; script = script.getGeneScript()) {
						if (!script.isCirc()) {
							continue;
						}
						for (Transcript s : linear_scripts) {
							if (script.getScript_id().equals(s.getScript_id())) {
								if (used_scripts.add(script)) {
									Transcript peak_script = getPeakScript(s, peak, chr);
									peak.addScript(peak_script);
									peak_script.setMut(peak.getScript().getMut());
									++peak_num[1];
								}
								break;
							}
						}
					}
				}
				if (peak.getScript() == null && args.isLinearPeak()) {
					Transcript script = biggestScript(linear_scripts);
					peak.addScript(getPeakScript(script, peak, chr));
					if (sim_flag[2] && script != null) {
						mutSim(mut_tree, script, peak);
					}
					++peak_num[1];
				}
			}
//			for (Node<Gene> gene_node : chr.getGeneTree()) {
//				Gene gene = gene_node.getValue();
//				int enrichIndex = gene.getPeakCount() > 0 ? 0 : 1;
//				Gene enrich = CommonMethod.randElement(gene_enrich == null ? null : gene_enrich.get(enrichIndex));
//				gene.setIp_enrich(enrich == null ? 0.01 / (Math.random() + 0.01) : enrich.getIp_enrich());
//			}
			if (sim_flag[2]) {
				for (Node<Mutation> mut_node : muts.getOrDefault(chr_string, new IntervalTree<>())) {
					for (Iterator<Node<Exon>> exon_nodes = chr.getExonTree().overlappers(mut_node.getStart(), mut_node.getStart());
							exon_nodes.hasNext();) {
						Node<Exon> exon_node = exon_nodes.next();
						for (Exon exon = exon_node.getValue(); exon != null; exon = exon.getAnotherExon()) {
							exon.getScript().setOnceMut(mut_node.getValue());
						}
					}
				}
			}
		});
		if (sim_flag[2]) {
			try {
				Method.writeFile(args.getMutFile(), muts, true, false, Mutation.getHeader());
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		Method.printNow("Circle peak scripts: " + peak_num[0]);
		Method.printNow("Linear peak scripts: " + peak_num[1]);
	}
	
//	private CircleClip circSim(Transcript script) {
//		return null;
//	}
	
	private void mutSim(IntervalTree<Mutation> mut_tree, Transcript script, Peak peak) {
		if (script == null) {
			return;
		}
		int read_len = args.getReadLength();
		Chromosome chr = genome.getChr(script.getChr());
		if (chr == null || chr.getSeq() == null) {
			return;
		}
		if (peak == null) {
			int index = (int) (script.getExonLength() * Math.random());
			int genome_index = SeqMethod.getGenomeIndex(script, index) + 1;
			Mutation mut = buildMut(genome_index - 1, chr);
			mut.setDescription("Exon");
			mut_tree.put(genome_index, genome_index, mut);
			mut.setMaf(sampleDiscreteMut(script.getGene().getSimMaf()));
			mut.setScript(script);
			script.setOnceMut(mut);
		}
		else if (peak.isCirc()) {
			if (script.getExonLength() < read_len * 2) {
				int index = (int) (script.getExonLength() * Math.random());
				int genome_index = SeqMethod.getGenomeIndex(script, index) + 1;
				Mutation mut = buildMut(genome_index - 1, chr);
				mut.setDescription("Peak");
				mut_tree.put(genome_index, genome_index, mut);
				mut.setMaf(sampleDiscreteMut(script.getGene().getSimMaf()));
				mut.setScript(script);
				script.setOnceMut(mut);
				for (Transcript s = peak.getScript(); s != null; s = s.getGeneScript()) {
					s.setOnceMut(mut);
				}
			}
			else {
				int seq_index = SeqMethod.getSeqIndex(script, peak.getStart());
				if (seq_index - read_len + 1 < 0) {
					int seq_start = seq_index + read_len;
					int index = (int) ((script.getExonLength() - 2 * read_len + 1) * Math.random()) + seq_start;
					int genome_index = SeqMethod.getGenomeIndex(script, index) + 1;
					Mutation mut = buildMut(genome_index - 1, chr);
					mut.setDescription("Exon");
					mut_tree.put(genome_index, genome_index, mut);
					mut.setMaf(sampleDiscreteMut(script.getGene().getSimMaf()));
					mut.setScript(script);
					
					index = (int) ((2 * read_len - 1) * Math.random());
					index = index >= seq_start ? script.getExonLength() - 2 * read_len + 1 + index : index;
					genome_index = SeqMethod.getGenomeIndex(script, index) + 1;
					mut = buildMut(genome_index - 1, chr);
					mut.setDescription("Peak");
					mut_tree.put(genome_index, genome_index, mut);
					mut.setMaf(sampleDiscreteMut(script.getGene().getSimMaf()));
					mut.setScript(script);
					script.setOnceMut(mut);
					for (Transcript s = peak.getScript(); s != null; s = s.getGeneScript()) {
						s.setOnceMut(mut);
					}
				}
				else if (seq_index + read_len > script.getExonLength()) {
					int seq_start = seq_index + read_len - script.getExonLength();
					int index = (int) ((script.getExonLength() - 2 * read_len + 1) * Math.random()) + seq_start;
					int genome_index = SeqMethod.getGenomeIndex(script, index) + 1;
					Mutation mut = buildMut(genome_index - 1, chr);
					mut.setDescription("Exon");
					mut_tree.put(genome_index, genome_index, mut);
					mut.setMaf(sampleDiscreteMut(script.getGene().getSimMaf()));
					mut.setScript(script);
					
					index = (int) ((2 * read_len - 1) * Math.random());
					index = index >= seq_start ? script.getExonLength() - 2 * read_len + 1 + index : index;
					genome_index = SeqMethod.getGenomeIndex(script, index) + 1;
					mut = buildMut(genome_index - 1, chr);
					mut.setDescription("Peak");
					mut_tree.put(genome_index, genome_index, mut);
					mut.setMaf(sampleDiscreteMut(script.getGene().getSimMaf()));
					mut.setScript(script);
					script.setOnceMut(mut);
					for (Transcript s = peak.getScript(); s != null; s = s.getGeneScript()) {
						s.setOnceMut(mut);
					}
				}
				else {
					int seq_start = seq_index - read_len + 1;
					int index = (int) ((2 * read_len - 1) * Math.random()) + seq_start;
					int genome_index = SeqMethod.getGenomeIndex(script, index) + 1;
					Mutation mut = buildMut(genome_index - 1, chr);
					mut.setDescription("Peak");
					mut_tree.put(genome_index, genome_index, mut);
					mut.setMaf(sampleDiscreteMut(script.getGene().getSimMaf()));
					mut.setScript(script);
					script.setOnceMut(mut);
					for (Transcript s = peak.getScript(); s != null; s = s.getGeneScript()) {
						s.setOnceMut(mut);
					}
					
					index = (int) ((script.getExonLength() - 2 * read_len + 1) * Math.random());
					index = index >= seq_start ? 2 * read_len - 1 + index : index;
					genome_index = SeqMethod.getGenomeIndex(script, index) + 1;
					mut = buildMut(genome_index - 1, chr);
					mut.setDescription("Exon");
					mut_tree.put(genome_index, genome_index, mut);
					mut.setMaf(sampleDiscreteMut(script.getGene().getSimMaf()));
					mut.setScript(script);
				}
			}
		}
		else {
			int seq_index = SeqMethod.getSeqIndex(script, peak.getStart());
			int seq_start = Math.max(0, seq_index - read_len + 1);
			int seq_end = Math.min(script.getExonLength(), seq_index + read_len);
			int index = (int) ((seq_end - seq_start) * Math.random()) + seq_start;
			int genome_index = SeqMethod.getGenomeIndex(script, index) + 1;
			Mutation mut = buildMut(genome_index - 1, chr);
			mut.setDescription("Peak");
			mut_tree.put(genome_index, genome_index, mut);
			mut.setMaf(sampleDiscreteMut(script.getGene().getSimMaf()));
			mut.setScript(script);
			script.setOnceMut(mut);
			for (Transcript s = peak.getScript(); s != null; s = s.getGeneScript()) {
				s.setOnceMut(mut);
			}
			
			if (seq_end - seq_start < script.getExonLength()) {
				index = (int) ((script.getExonLength() - seq_end + seq_start) * Math.random());
				index = index >= seq_start ? seq_end - seq_start + index : index;
				genome_index = SeqMethod.getGenomeIndex(script, index) + 1;
				mut = buildMut(genome_index - 1, chr);
				mut.setDescription("Exon");
				mut_tree.put(genome_index, genome_index, mut);
				mut.setMaf(sampleDiscreteMut(script.getGene().getSimMaf()));
				mut.setScript(script);
			}
		}
	}
	
	private Mutation buildMut(int index, Chromosome chr) {
		char ref = Character.toUpperCase(chr.getSeq().charAt(index));
		char alt = ref;
		while (alt == ref) {
			alt = base_units[(int) (Math.random() * base_units.length)];
		}
		return new Mutation(ref, new char[] {alt}, null);
	}
	
	private Transcript biggestScript(Collection<Transcript> scripts) {
		if (scripts == null || scripts.size() < 1) {
			return null;
		}
		Transcript script = null;
		for (Transcript s : scripts) {
			script = script == null || script.getExonLength() < s.getExonLength() ? s : script;
		}
		return script;
	}
	
	private Transcript uniqPeakScript(Transcript script, Peak peak, IntervalTree<Peak> peak_tree, Chromosome chr) {
		if (!args.isPeakSim() && script.isPeak()) {
			if (peak.getScript() == null) {
				for (Exon exon = script.getExons(); exon != null; exon = exon.getDownStream()) {
					for (Iterator<Node<Peak>> peak_nodes = peak_tree.overlappers(exon.getStart() + 1, exon.getEnd()); peak_nodes.hasNext();) {
						Node<Peak> peak_node = peak_nodes.next();
						if (peak_node.getValue().getScript() != null) {
							return peak_node.getValue().getScript();
						}
					}
				}
			}
			else {
				return null;
			}
		}
		script.setPeak(true);
		return getPeakScript(script, peak, chr);
	}
	
	private void setCircle(CircleClip circ, Transcript script, Chromosome chr) {
		boolean exon_flag = false;
		if (script != null) {
//			if (script.getScript_id().equals("ENST00000627315.2")) {
//				System.out.println("debug");
//			}
			Transcript circ_script = new Transcript(script.getGene(), script.getScript_id(), script.getScript_type(), 0, 0);
			circ_script.setCirc(circ);
			circ_script.setCircTrans(true);
			boolean add_flag = false;
			for (Exon start_exon = script.getExons(); start_exon != null; start_exon = start_exon.getDownStream()) {
				if (add_flag || (start_exon.getStart() <= circ.getStart() && (start_exon.getDownStream() == null || start_exon.getDownStream().getStart() > circ.getStart()))) {
					circ_script.addLastExon("exon", new Exon(circ_script, start_exon.getStart(), start_exon.getEnd()));
					exon_flag = add_flag ? exon_flag : start_exon.getStart() == circ.getStart();
					add_flag = true;
				}
				if (start_exon.getEnd() >= circ.getEnd()) {
//					if (!add_flag) {
//						System.out.printf("SID:%s CSID:%s CIRC:%s %s\n", script.getScript_id(), circ_script.getScript_id(), chr.getChr(), circ.toString());
//					}
					circ_script.sort();
					exon_flag = exon_flag && start_exon.getEnd() == circ.getEnd();
					break;
				}
			}
			circ.setScript(circ_script);
			script.addGeneScript(circ_script);
		}
//		else {
//			System.out.println("debug");
//		}
		circ.setExon(exon_flag);
		circ.setGTAG(chr != null && isGTAG(chr.getSeq().substring(circ.getEnd(), circ.getEnd() + 2), 
				chr.getSeq().substring(circ.getStart() - 2, circ.getStart())));
	}
	
	private double sampleDiscreteMut(double maf) {
		if (Math.random() < 0.5) {
			return maf;
		}
		double f = args.getFixAlleleFre();
		if (f >= 0.0 && f <= 0.5) {
			return maf + (Math.random() < 0.5 ? f : -f);
		}
		int s = (int) Math.floor(maf / 0.1);
		UniformIntegerDistribution uid = new UniformIntegerDistribution(-s, 9 - s);
		return maf + 0.1 * uid.sample();
	}
	
	private Transcript getPeakScript(Transcript script, Peak peak, Chromosome chr) {
		Transcript out = null;
		if (script != null) {
			out = new Transcript(script.getGene(), script.getScript_id(), script.getScript_type(), 0, 0);
			script.setExpression(1.0);
			out.setExpression(script.getExpression());
			out.setMut(script.getMut());
			out.setPeak(true);
			int index = 0;
			for (Exon exon = script.getExons(); exon != null; exon = exon.getDownStream()) {
				if (exon.getStart() <= peak.getStart() && exon.getEnd() >= peak.getEnd()) {
					index += peak.getStart() - exon.getStart();
					break;
				}
				index += exon.getLength();
			}
			if (index < script.getExonLength()) {
				if (!args.isPeakSim()) {
					for (Exon exon = script.getExons(); exon != null; exon = exon.getDownStream()) {
						out.addLastExon("exon", new Exon(out, exon.getStart(), exon.getEnd()));
					}
					out.setCircTrans(script.isCircTrans());
					out.setCirc(script.getCirc());
					return out;
				}
				int len = 0;
				if (script.isCirc()) {
					if (script.getExonLength() > args.getReadLength()) {
						if (index - args.getReadLength() < 0) {
							Exon head_exon = new Exon(null, 0, 0);
							for (Exon exon = script.getExons(); exon != null; exon = exon.getDownStream()) {
								if (len <= index + args.getReadLength()) {
									out.addLastExon("exon", new Exon(out, exon.getStart(),
										exon.getStart() + Math.min(exon.getLength(), index + args.getReadLength() - len)));
								}
								if (len + exon.getLength() > index - args.getReadLength() + script.getExonLength()) {
									head_exon.addLastExon(new Exon(out, exon.getStart() + Math.max(0, index + script.getExonLength() - args.getReadLength() - len),
											exon.getEnd()));
								}
								len += exon.getLength();
							}
							out.addFirstExon("exon", head_exon.getDownStream());
							out.setCirc(script.getCirc());
							return out;
						}
						else if (index + args.getReadLength() > script.getExonLength()) {
							Exon head_exon = new Exon(null, 0, 0);
							for (Exon exon = script.getExons(); exon != null; exon = exon.getDownStream()) {
								if (len <= index + args.getReadLength() - script.getExonLength()) {
									out.addLastExon("exon", new Exon(out, exon.getStart(),
										exon.getStart() + Math.min(exon.getLength(), index + args.getReadLength() - script.getExonLength() - len)));
								}
								if (len + exon.getLength() > index - args.getReadLength()) {
									head_exon.addLastExon(new Exon(out, exon.getStart() + Math.max(0, index - args.getReadLength() - len),
											exon.getEnd()));
								}
								len += exon.getLength();
							}
							out.addFirstExon("exon", head_exon.getDownStream());
							out.setCirc(script.getCirc());
							return out;
						}
					}
					else {
						for (Exon exon = script.getExons(); exon != null; exon = exon.getDownStream()) {
							out.addLastExon("exon", new Exon(out, exon.getStart(), exon.getEnd()));
						}
						out.setCircTrans(script.isCircTrans());
						out.setCirc(script.getCirc());
						return out;
					}
				}
				for (Exon exon = script.getExons(); exon != null; exon = exon.getDownStream()) {
					if (len + exon.getLength() > index - args.getReadLength() 
							&& len <= index + args.getReadLength()) {
						out.addLastExon("exon", new Exon(out, exon.getStart() + Math.max(index - args.getReadLength() - len, 0),
								exon.getStart() + Math.min(exon.getLength(), index + args.getReadLength() - len)));
					}
					len += exon.getLength();
				}
				out.setCirc(script.getCirc());
				return out;
			}
		}
		Gene gene = null;
		for (Iterator<Node<Gene>> gene_nodes = chr.getGeneTree().overlappers(peak.getEnd(), peak.getEnd()); gene_nodes.hasNext();) {
			Node<Gene> gene_node = gene_nodes.next();
			gene = gene == null ? gene_node.getValue() : (gene.getLength() < gene_node.getLength() ? gene_node.getValue() : gene);
		}
		gene = gene == null ? new Gene(chr, "None", "None", "None", '.', peak.getStart() - args.getReadLength(), peak.getEnd() + args.getReadLength()) : gene;
		out = new Transcript(gene, "None", "None", peak.getStart() - args.getReadLength(), peak.getEnd() + args.getReadLength());
		out.addLastExon("exon", new Exon(out, out.getStart(), out.getEnd()));
		out.setPeak(true);
		return out;
	}
	
	private void prepareSim() {
		final int[] script_num = {0, 0, 0};
		base_num = 0;
		read_num = 0.0;
		read_num_ip = 0.0;
		if (args.isLinearPeak()) {
			genome.getChrMap().forEach((chr_string, chr) -> {
				for (Node<Gene> gene_node : chr.getGeneTree()) {
					Transcript script = gene_node.getValue().getScript();
					for (Transcript s = gene_node.getValue().getScript(); s != null; s = s.getGeneScript()) {
						if (s.isCirc()) {
							script_num[1] += Math.round(s.getExonLength() * s.getExpression());
							script_num[0]++;
						}
						else if (script.getExonLength() < s.getExonLength()) {
							script = s;
						}
					}
					if (script != null) {
						script_num[2] += Math.round(script.getExonLength() * script.getExpression());
						script_num[0]++;
					}
					if (args.isPeakSim()) {
						
						gene_node.getValue().setIp_enrich(calEnrich(2 * args.getReadLength(),
								script.getExonLength(), enrich_dis.cumulativeProbability(gene_node.getValue().getIp_enrich())));
					}
//					gene_node.getValue().setIp_enrich(enrich_dis.cumulativeProbability(gene_node.getValue().getIp_enrich()));
//					gene_node.getValue().setIp_enrich(enrich_dis.inverseCumulativeProbability(gene_node.getValue().getIp_enrich()));
					read_num += gene_node.getValue().getInputReadCount();
					read_num_ip += gene_node.getValue().getIpReadCount();
				}
			});
			linear_ratio = args.getLinearProp() / script_num[2] * script_num[1];
			Method.printNow("Linear Ratio: " + linear_ratio);
			base_num = linear_ratio > 0 ? (int) (script_num[1] / linear_ratio) : 0;
			base_num += script_num[2];
		}
		else {
			circ_table.forEach((chr, circ_tree) -> {
				for (Node<CircleClip> circ_node : circ_tree) {
					base_num += circ_node.getValue().getScript() != null ? circ_node.getValue().getScript().getExonLength() : 0;
				}
			});
		}
		
//		peak_table.forEach((chr_string, peak_tree) -> {
//			HashSet<Transcript> used_scripts = new HashSet<>();
//			HashMap<Gene, Integer> gene_map = new HashMap<>();
//			for (Node<Peak> peak_node : peak_tree) {
//				for (Transcript script = peak_node.getValue().getScript(); script != null; script = script.getGeneScript()) {
//					if (used_scripts.add(script)) {
//						gene_map.put(script.getGene(), gene_map.getOrDefault(script.getGene(), 0) + 1);
//					}
//				}
//			}
//			gene_map.forEach((gene, peak_num) -> {
//				gene.setIp_enrich(gene.getIp_enrich() / peak_num);
//				read_num_ip += gene.getIp_enrich() * gene.getReadCount();
//			});
//		});
		Method.printNow("Total Bases: " + base_num + " in " + script_num[0]);
		Method.printNow("Total Circle Bases: " + script_num[1]);
		Method.printNow("Total Linear Bases: " + script_num[2]);
	}
	
	private double calEnrich(int peak, int script, double enrich) {
		return enrich;
	}
	
	private void sim() {
		try (final BufferedWriter writer_input1 = buildWriter(args.getOutPrefix() + "_C_L1.fastq", false, args.isOutZip());
				final BufferedWriter writer_input2 = buildWriter(args.getOutPrefix() + "_C_L2.fastq", false, args.isOutZip());
				final BufferedWriter writer_ip1 = buildWriter(args.getOutPrefix() + "_T_L1.fastq", false, args.isOutZip());
				final BufferedWriter writer_ip2 = buildWriter(args.getOutPrefix() + "_T_L2.fastq", false, args.isOutZip());
				final BufferedWriter writer_annote = args.isAnnote() ? 
						buildWriter(args.getOutPrefix() + "_annote.txt", false, args.isOutZip()) : null) {
			
			if (args.getEnrich() < 0.0) {
				final double[] total_count = {0.0}; 
				peak_table.entrySet().forEach(entry -> {
					HashSet<Transcript> used_scripts = new HashSet<>();
					for (Node<Peak> peak_node : entry.getValue()) {
						for (Transcript script = peak_node.getValue().getScript(); script != null; script = script.getGeneScript()) {
							if (used_scripts.add(script)) {
								total_count[0] += script.getExonLength() / (script.isCirc() ? linear_ratio : 1.0);
							}
						}
					}
				});
				args.setEnrich((double) base_num / (double) total_count[0]);
			}
			final int[] total_count = {0, 0, 0};
			if (args.isLinearPeak()) {
//				genome.getChrMap().entrySet().parallelStream().forEach(entry -> {
				genome.getChrMap().entrySet().forEach(entry -> {
					for (Node<Gene> gene_node : entry.getValue().getGeneTree()) {
						Transcript script = gene_node.getValue().getScript();
						for (Transcript s = gene_node.getValue().getScript();s != null; s = s.getGeneScript()) {
							if (s.isCirc()) {
								double count = read_num > 10.0 ? getReadGeneCount(gene_node.getValue(), false, false) 
										: (getBaseScriptCount(s) / linear_ratio);							
								int read_count = writeWithAnnote(writer_input1, writer_input2, writer_annote, buildFastqs(s, count, false));
								synchronized (total_count) {
									total_count[0] += read_count;
								}
								if (!s.isPeak()) {
									count = read_num_ip > 10.0 ? getReadGeneCount(gene_node.getValue(), true, false)
											: Math.max(args.getNonPeakBack() * count, args.getMinBackRead());
									read_count = writeWithAnnote(writer_ip1, writer_ip2, writer_annote, buildFastqs(s, count, true));
									synchronized (total_count) {
										total_count[1] += read_count;
									}
								}
							}
							else if (script.getExonLength() < s.getExonLength()) {
								script = s;
							}
						}
						if (script != null) {
							double count = read_num > 10.0 ? getReadGeneCount(gene_node.getValue(), false, false) : getBaseScriptCount(script);
							int read_count = writeWithAnnote(writer_input1, writer_input2, writer_annote, buildFastqs(script, count, false));
							synchronized (total_count) {
								total_count[0] += read_count;
							}
							count = read_num_ip > 10.0 ? getReadGeneCount(gene_node.getValue(), true, false)
									: Math.max(args.getNonPeakBack() * count, args.getMinBackRead());
							read_count = writeWithAnnote(writer_ip1, writer_ip2, writer_annote, buildFastqs(script, count, true));
							synchronized (total_count) {
								total_count[2] += read_count;
							}
						}
					}
				});
			}
			else {
				circ_table.entrySet().parallelStream().forEach(entry -> {
					for (Node<CircleClip> circ_node : entry.getValue()) {
						for (Transcript script = circ_node.getValue().getScript(); script != null; script = null) {
							double count = getBaseScriptCount(script);
							int read_count = writeWithAnnote(writer_input1, writer_input2, writer_annote, buildFastqs(script, count, false));
							synchronized (total_count) {
								total_count[0] += read_count;
							}
							if (!script.isPeak()) {
								count = args.getNonPeakBack() * count;
								read_count = writeWithAnnote(writer_ip1, writer_ip2, writer_annote, buildFastqs(script, count, true));
								synchronized (total_count) {
									total_count[1] += read_count;
								}
							}
						}
					}
				});
			}
			Method.printNow("Total Input Count: " + total_count[0]);
			Method.printNow("Total Circle Count: " + total_count[1]);
			Method.printNow("Total Linear Count: " + total_count[2]);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void simPeak() {
		try (final BufferedWriter writer_ip1 = buildWriter(args.getOutPrefix() + "_T_L1.fastq", true, args.isOutZip());
				final BufferedWriter writer_ip2 = buildWriter(args.getOutPrefix() + "_T_L2.fastq", true, args.isOutZip());
				final BufferedWriter writer_annote = args.isAnnote() ? 
						buildWriter(args.getOutPrefix() + "_annote.txt", true, args.isOutZip()) : null){
			Method.printNow("Peak start");
			System.out.println("Enrichment is " + args.getEnrich());
			final int[] total_count = {0, 0, 0};
//			peak_table.entrySet().parallelStream().forEach(entry -> {
			peak_table.entrySet().forEach(entry -> {
				HashSet<Transcript> used_scripts = new HashSet<>();
				for (Node<Peak> peak_node : entry.getValue()) {
					for (Transcript script = peak_node.getValue().getScript(); script != null; script = script.getGeneScript()) {
						if (used_scripts.add(script)) {
							double count = (read_num_ip > 10.0 ? getReadGeneCount(script.getGene(), true, true) : 
								(getBaseScriptCount(script) * args.getEnrich())) / (script.isCirc() ? linear_ratio : 1.0);
							int read_count = writeWithAnnote(writer_ip1, writer_ip2, writer_annote, buildFastqs(script, count, true));
							synchronized (total_count) {
								if (script.isCirc()) {
									total_count[0] += script.getExonLength();
								}
								else {
									total_count[1] += script.getExonLength();
								}
								total_count[2] += read_count;
							}
						}
					}
				}
			});
			Method.printNow("Total IP Circle Bases: " + total_count[0]);
			Method.printNow("Total IP Linear Bases: " + total_count[1]);
			Method.printNow("Total Count: " + total_count[2]);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void simScaf() throws IOException {
		try (final BufferedWriter writer_input1 = buildWriter(args.getOutPrefix() + "_C_L1.fastq", false, args.isOutZip());
				final BufferedWriter writer_input2 = buildWriter(args.getOutPrefix() + "_C_L2.fastq", false, args.isOutZip());
				final BufferedWriter writer_ip1 = buildWriter(args.getOutPrefix() + "_T_L1.fastq", false, args.isOutZip());
				final BufferedWriter writer_ip2 = buildWriter(args.getOutPrefix() + "_T_L2.fastq", false, args.isOutZip());
				final BufferedWriter writer_annote = args.isAnnote() ? 
						buildWriter(args.getOutPrefix() + "_annote.txt", false, args.isOutZip()) : null) {
			
			if (args.getEnrich() < 0.0) {
				final double[] total_count = {0.0}; 
				peak_table.entrySet().forEach(entry -> {
					HashSet<Transcript> used_scripts = new HashSet<>();
					for (Node<Peak> peak_node : entry.getValue()) {
						for (Transcript script = peak_node.getValue().getScript(); script != null; script = script.getGeneScript()) {
							if (used_scripts.add(script)) {
								total_count[0] += script.getExonLength() / (script.isCirc() ? linear_ratio : 1.0);
							}
						}
					}
				});
				args.setEnrich((double) base_num / (double) total_count[0]);
			}
			final int[] total_count = {0, 0, 0};
//			genome.getChrMap().entrySet().parallelStream().forEach(entry -> {
			genome.getChrMap().entrySet().forEach(entry -> {
				for (Node<Gene> gene_node : entry.getValue().getGeneTree()) {
					Transcript script = gene_node.getValue().getScript();
					for (Transcript s = gene_node.getValue().getScript();s != null; s = s.getGeneScript()) {
						script = s.getExonLength() > script.getExonLength() ? s : script;
					}
					double count = read_num > 10.0 ? getReadGeneCount(gene_node.getValue(), false, false) : getBaseScriptCount(script);
					int read_count = writeWithAnnote(writer_input1, writer_input2, writer_annote, buildFastqs(script, count, false));
					synchronized (total_count) {
						total_count[0] += read_count;
					}
					count = read_num_ip > 10.0 ? getReadGeneCount(gene_node.getValue(), true, false)
							: Math.max(args.getNonPeakBack() * count, args.getMinBackRead());
					read_count = writeWithAnnote(writer_ip1, writer_ip2, writer_annote, buildFastqs(script, count, true));
					synchronized (total_count) {
						total_count[2] += read_count;
					}
				}
			});
			Method.printNow("Total Input Count: " + total_count[0]);
			Method.printNow("Total Circle Count: " + total_count[1]);
			Method.printNow("Total Linear Count: " + total_count[2]);
		}
	}
	
	private double getBaseScriptCount(Transcript script) {
		return getDepth(false) * script.getExonLength() * script.getExpression();
	}
	
	private double getReadGeneCount(Gene gene, boolean ip_flag, boolean peak) {
		return getDepth(true) * (ip_flag ? (peak ? gene.getIp_enrich() : 1.0 - gene.getIp_enrich())
				* gene.getIpReadCount() / read_num_ip : (gene.getInputReadCount() / read_num));
	}
	
	private double getDepth(boolean read_flag) {
		return read_flag ? args.getTotalReads() : args.getReadDepth() < 0.0 ? 
				args.getTotalReads() / base_num : (args.getReadDepth() / args.getReadLength());
	}
	
	public List<Fastq> buildFastqs(Transcript script, double count, boolean ip) {
		List<Fastq> out = new ArrayList<>();
		if (script == null || script.getGene().getChromosome().getSeq() == null) {
			return out;
		}
		List<String> seqs = getRTmutExonBases(script);// {ref, alt}
		UniformIntegerDistribution uid = new UniformIntegerDistribution(0, script.getExonLength() - 1);
		BinomialDistribution ald = new BinomialDistribution((int) (args.getAlignmentLength() / 0.95 + 0.99), 0.95);
		BinomialDistribution rld = new BinomialDistribution(args.getReadLength() * 2, 0.5);
		boolean[] flags = new boolean[] {args.isPairEnd(), false, script.isPeak(), 
				false, false, ip, false};
		if (args.isRandDepth()) {
			count = sampleReadCount(count);
		}
		count -= 0.5;
//		counts = new int[149];
		for (; count >= 0.0; count--) {
			int start = uid.sample();
			int align_len = Math.min(args.getAlignmentLength(), ald.sample());
			int read_len = Math.max(rld.sample(), align_len);
			int end = start + read_len + 1;
			flags[0] = args.isPairEnd();
			flags[1] = false;
			flags[3] = false;
			flags[4] = false;
			flags[6] = false;
			if (!scaf_table.isEmpty() && !args.isRtMut() && ip) {
				flags[4] = CommonMethod.randInt(2) == 1;
				flags[6] = true;
				int[] region = getRTstopRegion(script, start, end, seqs, !flags[4] && !args.isNegative());
				start = region[0];
				end = region[1];
			}
			if (script.isCirc() && script.isCircTrans()) {
				int[] regions = flags[0] ? new int[]{start, start + align_len, end - align_len, end} : new int[]{start, end};
				String annote = buildAnnote(script, regions, flags);
				String id = buildUniqID(script.getChr(), annote);
				String exon_seq = flags[3] && flags[4] ? seqs.get(1) : seqs.get(0);
				annote = args.isAnnote() ? id + "\t" + annote : null;
				for (int i = 0; i < regions.length; i += 2) {
					int round = (regions[i + 1] + exon_seq.length() - 1) / exon_seq.length() - 
							(regions[i] + exon_seq.length()) / exon_seq.length() ;
					if (round == 0) {
						out.add(new Fastq(id, reverseFliq(getRtmutBases(exon_seq.substring(regions[i] % exon_seq.length(),
								(regions[i + 1] - 1) % exon_seq.length()), ip), i / 2), annote));
					}
					else {
						StringBuilder sb = new StringBuilder();
						sb.append(exon_seq.substring((regions[i] + exon_seq.length()) % exon_seq.length()));
						while (--round > 0) {
							sb.append(exon_seq);
						}
						sb.append(exon_seq.substring(0, (regions[i + 1] + exon_seq.length() - 1) % exon_seq.length() + 1));
						out.add(new Fastq(id, reverseFliq(getRtmutBases(sb.toString(), ip), i / 2), annote));
					}
				}
			}
			else {
				if (end > script.getExonLength()) {
					if (script.getExonLength() <= align_len) {
						int[] regions = flags[0] ? new int[]{0, script.getExonLength(), 0, script.getExonLength()} : new int[]{0, script.getExonLength()};
						String annote = buildAnnote(script, regions, flags);
						String id = buildUniqID(script.getChr(), annote);
						String exon_seq = flags[3] && flags[4] ? seqs.get(1) : seqs.get(0);
						annote = args.isAnnote() ? id + "\t" + annote : null;
						out.add(new Fastq(id, getRtmutBases(exon_seq.substring(regions[0], regions[1]), ip), annote));
						if (args.isPairEnd()) {
							out.add(new Fastq(id, reverseFliq(getRtmutBases(exon_seq.substring(regions[2], regions[3]), ip)), annote));
						}
					}
					else {
						if (end - script.getExonLength() < script.getExonLength() - start) {
							end = script.getExonLength();
						}
						else {
							end = Math.min(script.getExonLength(), end - script.getExonLength());
							start = 0;
						}
						int[] regions = flags[0] ? new int[]{start, Math.min(start + align_len, script.getExonLength()), 
								Math.max(0, end - align_len), end} : new int[]{start, Math.min(end, script.getExonLength())};
						String annote = buildAnnote(script, regions, flags);
						String id = buildUniqID(script.getChr(), annote);
						String exon_seq = flags[3] && flags[4] ? seqs.get(1) : seqs.get(0);
						annote = args.isAnnote() ? id + "\t" + annote : null;
						out.add(new Fastq(id, getRtmutBases(exon_seq.substring(regions[0], regions[1]), ip), annote));
						if (args.isPairEnd()) {
							out.add(new Fastq(id, reverseFliq(getRtmutBases(exon_seq.substring(regions[2], regions[3]), ip)), annote));
						}
					}
//					end = end > 2 * exon_seq.length() ? exon_seq.length() : end - exon_seq.length();
//					if (exon_seq.length() - start < end) {
//						start = Math.max(args.getMinAlignmentLength(), end - read_len + align_len);
//						int[] regions = flags[0] ? new int[]{0, start, Math.max(0, end - align_len), end} : new int[]{0, start};
//						String annote = buildAnnote(script, regions, flags);
//						String id = buildUniqID(script.getChr(), false, flags[1], flags[2], ip);
//						out.add(new Fastq(id, exon_seq.substring(regions[0], regions[1]), id + "\t" + annote));
//						if (args.isPairEnd()) {
//							out.add(new Fastq(id, reverseFliq(exon_seq.substring(regions[2], regions[3])), id + "\t" + annote));
//						}
//					}
//					else {
//						end = Math.min(exon_seq.length() - args.getMinAlignmentLength(), end + exon_seq.length() - align_len);
//						int[] regions = flags[0] ? new int[]{start, Math.min(exon_seq.length(), start + align_len), end, exon_seq.length()} :
//							new int[]{start, Math.min(exon_seq.length(), start + align_len)};
//						String annote = buildAnnote(script, regions, flags);
//						String id = buildUniqID(script.getChr(), false, flags[1], flags[2], ip);
//						out.add(new Fastq(id, exon_seq.substring(regions[0], regions[1]), id + "\t" + annote));
//						if (args.isPairEnd()) {
//							out.add(new Fastq(id, reverseFliq(exon_seq.substring(regions[2], regions[3])), id + "\t" + annote));
//						}
//					}
				}
				else {
					int[] regions = flags[0] ? new int[]{start, Math.min(start + align_len, script.getExonLength()), Math.min(start, end - align_len), end} :
						new int[]{start, end};
					String annote = buildAnnote(script, regions, flags);
					String id = buildUniqID(script.getChr(), annote);
					String exon_seq = flags[3] && flags[4] ? seqs.get(1) : seqs.get(0);
					annote = args.isAnnote() ? id + "\t" + annote : null;
					out.add(new Fastq(id, getRtmutBases(exon_seq.substring(regions[0], regions[1]), ip), annote));
					if (args.isPairEnd()) {
						out.add(new Fastq(id, reverseFliq(getRtmutBases(exon_seq.substring(regions[2], regions[3]), ip)), annote));
					}
				}
			}
		}
//		Method.logGlobal(Arrays.toString(counts));
		return out;
	}
	
	private String buildAnnote(Transcript script, int[] regions, boolean[] flags) {
		if (script == null) {
			return null;
		}
		StringBuilder sb = new StringBuilder();
		char sep = '\t';
		sb.append(script.getChr());
		sb.append(sep);
		sb.append(script.getScript_id());
		sb.append(sep);
		if (script.isCirc() && script.isCircTrans()) {
			for (int i = 0; i < regions.length; i += 2) {
				if (regions[i] < 0) {
					if (regions[i + 1] > script.getExonLength()) {
						flags[1] = true;
						sb.append(buildExonAnnote(script, script.getExonLength() + regions[i], script.getExonLength(), flags));
						sb.append(':');
						sb.append(buildExonAnnote(script, 0, script.getExonLength(), flags));
						sb.append(':');
						sb.append(buildExonAnnote(script, 0, regions[i + 1] - script.getExonLength(), flags));
						sb.append(sep);
					}
					else {
						flags[1] = true;
						sb.append(buildExonAnnote(script, script.getExonLength() + regions[i], script.getExonLength(), flags));
						sb.append(':');
						sb.append(buildExonAnnote(script, 0, regions[i + 1], flags));
						sb.append(sep);
					}
				}
				else {
					if (regions[i + 1] > script.getExonLength()) {
						if (regions[i] > script.getExonLength()) {
							sb.append(buildExonAnnote(script, regions[i] - script.getExonLength(), regions[i + 1] - script.getExonLength(), flags));
							sb.append(sep);
						}
						else {
							flags[1] = true;
							sb.append(buildExonAnnote(script, regions[i], script.getExonLength(), flags));
							sb.append(':');
							sb.append(buildExonAnnote(script, 0, regions[i + 1] - script.getExonLength(), flags));
							sb.append(sep);
						}
					}
					else {
						sb.append(buildExonAnnote(script, regions[i], regions[i + 1], flags));
						sb.append(sep);
					}
				}
			}
		}
		else {
			for (int i = 0; i < regions.length; i += 2) {
				sb.append(buildExonAnnote(script, regions[i], regions[i + 1], flags));
				sb.append(sep);
			}
			flags[1] = flags[0];
		}
		if (script.isCirc()) {
			script.getCirc().incCircNum(flags[1], flags[5]);
		}
		sb.append(script.isCirc() ? 'C' : 'L');
		sb.append(flags[1] ? 'C' : 'L');
		sb.append(flags[2] ? 'P' : 'U');
		sb.append(flags[3] ? 'A' : 'R');
		sb.append(flags[4] ? 'A' : 'R');
		sb.append(flags[5] ? ":T" : ":C");
		return sb.toString();
	}

	private String buildExonAnnote(Transcript script, int start, int end, boolean[] flag) {
		if (script == null || script.getExons() == null || flag == null || flag.length < 6) {
			return null;
		}
		StringBuilder sb = new StringBuilder();
		flag[0] = false;
		for (Exon exon = script.getExons(); exon != null && end > 0; exon = exon.getDownStream()) {
			if (exon.getLength() > start) {
				int exon_start = exon.getStart() + start + 1;
				int exon_end = Math.min(end + exon.getStart(), exon.getEnd());
				flag[0] = flag[0] || (script.isCirc() ? (end > exon.getLength() && exon.getEnd() == script.getCirc().getEnd())
						: isCirc(script.getChr(), exon_start, exon_end));
				flag[3] = flag[3] || (script.isMut() && muts.get(script.getChr()) != null && 
						muts.get(script.getChr()).overlappers(exon_start, exon_end).hasNext());
				if (flag[3] && !flag[6]) {
					flag[4] = flag[5] ? Math.random() < muts.get(script.getChr()).overlappers(exon_start, exon_end).next().getValue().getMaf()
							: Math.random() < muts.get(script.getChr()).overlappers(exon_start, exon_end).next().getValue().getScript().getGene().getSimMaf();
					flag[6] = true;
				}
				sb.append(exon_start);
				sb.append('-');
				sb.append(exon_end);
				sb.append(':');
				if (flag[3]) {
					for (Iterator<Node<Mutation>> mut_nodes = muts.get(script.getChr()).overlappers(exon_start, exon_end);
							mut_nodes.hasNext();) {
						Node<Mutation> mut_node = mut_nodes.next();
						if (flag[4]) {
							mut_node.getValue().incCount(mut_node.getValue().getSNP()[0], flag[5]);
						}
						else {
							mut_node.getValue().incCount(mut_node.getValue().getRef(), flag[5]);
						}
					}
				}
			}
			start = Math.max(0, start - exon.getLength());
			end = Math.max(0, end - exon.getLength());
		}
		if (sb.length() > 0) {
			sb.setLength(sb.length() - 1);
		}
		return sb.toString();
	}
	
	private String getMutExonBases(Transcript script) {
		if (script == null) {
			return null;
		}
		if (!script.isMut() || muts == null || muts.get(script.getChr()) == null) {
			return script.getExonBases();
		}
		StringBuilder sb = new StringBuilder();
		int len = 0;
		
		for (Exon e = script.getExons(); e != null; e = e.getDownStream()) {
			sb.append(script.getGene().getChromosome().getSeq().substring(e.getStart(), e.getEnd()));
			for (Iterator<Node<Mutation>> mut_nodes = muts.get(script.getChr()).overlappers(e.getStart() + 1, e.getEnd());
					mut_nodes.hasNext();) {
				Node<Mutation> mut_node = mut_nodes.next();
				sb.setCharAt(len + mut_node.getStart() - 1 - e.getStart(), mut_node.getValue().getSNP()[0]);
			}
			len += e.getLength();
		}
		return sb.toString();
	}
	
	private List<String> getRTmutExonBases(Transcript script) {
		if (script == null) {
			return null;
		}
		List<String> seqs = new ArrayList<>();
		if (scaf_table == null || scaf_table.get(script.getChr()) == null) {
			seqs.add(script.getExonBases().toUpperCase());
			seqs.add(getMutExonBases(script).toUpperCase());
			return seqs;
		}
		StringBuffer sb1 = new StringBuffer();
		StringBuffer sb2 = new StringBuffer();
		StringBuilder str1 = new StringBuilder();
		StringBuilder str2 = new StringBuilder();
		boolean structFlag = false;
		String bases = script.getGene().getChromosome().getSeq();
		int len = 0;
		
		for (Exon e = script.getExons(); e != null; e = e.getDownStream()) {
			String exonSeq = bases.substring(e.getStart(), e.getEnd()).toUpperCase();
			sb1.append(exonSeq);
			sb2.append(exonSeq);
			for (int i = e.getStart(); i < e.getEnd(); ++i) {
				str1.append('(');
				str2.append('(');
			}
			for (Iterator<Node<Scaffold[]>> scaf_nodes = scaf_table.get(script.getChr()).overlappers(e.getStart() + 1, e.getEnd());
					scaf_nodes.hasNext();) {
				structFlag = true;
				Node<Scaffold[]> scaf_node = scaf_nodes.next();
				List<Integer> sites = scaf_node.getValue()[0].getUnpairSites();
				for (Integer site : sites) {
					int index = site - e.getStart() - 1;
					if (index >= 0 && index < e.getLength() && args.isVaildBase(sb1.charAt(len + index))) {
						sb1.setCharAt(index + len, Character.toLowerCase(sb1.charAt(index + len)));
						str1.setCharAt(index + len, '.');
					}
				}
				sites = scaf_node.getValue()[1].getUnpairSites();
				for (Integer site : sites) {
					int index = site - e.getStart() - 1;
					if (index >= 0 && index < e.getLength() && args.isVaildBase(sb2.charAt(len + index))) {
						sb2.setCharAt(index + len, Character.toLowerCase(sb1.charAt(index + len)));
						str2.setCharAt(index + len, '.');
					}
				}
			}
			for (Iterator<Node<Mutation>> mut_nodes = muts.get(script.getChr()).overlappers(e.getStart() + 1, e.getEnd());
					mut_nodes.hasNext();) {
				Node<Mutation> mut_node = mut_nodes.next();
				sb2.setCharAt(len + mut_node.getStart() - 1 - e.getStart(), mut_node.getValue().getSNP()[0]);
			}
			len += e.getLength();
		}
		if (args.isRtMut()) {
			seqs.add(sb1.toString());
			seqs.add(sb2.toString());
		}
		else {
			seqs.add(script.getExonBases().toUpperCase());
			seqs.add(getMutExonBases(script).toUpperCase());
		}
		if (structFlag) {
			seqs.add(str1.toString());
			seqs.add(str2.toString());
		}
		return seqs;
	}

	private String getRtmutBases(String bases, boolean ip) {
		if (!args.isRtMut() || bases.equals(bases.toUpperCase())) {
			return bases;
		}
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < bases.length(); ++i) {
			char c = bases.charAt(i);
			c = ip && Character.isLowerCase(c) && CommonMethod.randDouble() < 0.03 ? BamMethod.getMutBase(Character.toUpperCase(c)) : Character.toUpperCase(c);
			c = ip || CommonMethod.randDouble() >= 0.003 ? c : BamMethod.getMutBase(c);
			sb.append(c);
		}
		return sb.toString();
	}
	
	private int[] getRTstopRegion(Transcript script, int start, int end, List<String> seqs, boolean refFlag) {
		int[] region = new int[] {start,end};
		if (scaf_table == null || script == null || scaf_table.get(script.getChr()) == null || seqs.size() <= 2) {
			return region;
		}
		if (end > script.getExonLength()) {
			if (end - script.getExonLength() > script.getExonLength() - start) {
				start = 0;
				end = Math.min(script.getExonLength(), end - script.getExonLength());
			}
			else {
				end = script.getExonLength();
			}
		}
		region[0] = script.getGene().getStrand() == '-' ? start + args.getMinAlignmentLength() : start;
		region[1] = script.getGene().getStrand() == '-' ? end : end - args.getMinAlignmentLength();
		String struct = refFlag ? seqs.get(2) : seqs.get(3);
		String seq = refFlag ? seqs.get(0) : seqs.get(1);
		List<Integer> unpairSites = new ArrayList<>();
		for (int i = region[0]; i < region[1]; ++i) {
			if (struct.charAt(i) == '.' && args.isVaildBase(Character.toUpperCase(seq.charAt(i)))) {
				unpairSites.add(i + (script.getGene().getStrand() == '-' ? -1 : 1));
			}
		}
		if (unpairSites.isEmpty()) {
			region[0] = start;
			region[1] = end;
			return region;
		}
		unpairSites.add(script.getGene().getStrand() == '-' ? end : start);

//		int scriptStart = script.getExonSeqSite(start);
//		int scriptEnd = script.getExonSeqSite(end - 1);
//		int scafStart = start;
//		for (int i = 0; i < pairSites.size(); ++i) {
//			int scafSite = pairSites.get(i) > scriptStart && pairSites.get(i) < scriptEnd && CommonMethod.randDouble() < 0.01 ? script.inExons(pairSites.get(i)) : -1; 
//			scafStart = scafSite > 0 ? scafSite : scafStart;
//		}
		
		if (script.getGene().getStrand() == '-') {
//			region[1] = CommonMethod.randInt(2) == 0 ? CommonMethod.randElement(unpairSites) : end;
			region[1] = CommonMethod.randElement(unpairSites);
			region[0] = start;
		}
		else {
//			region[0] = CommonMethod.randInt(2) == 0 ? CommonMethod.randElement(unpairSites) : start;
			region[0] = CommonMethod.randElement(unpairSites);
			region[1] = end;
		}
		return region;
	}

	private BufferedWriter buildWriter(String file_name, boolean append, boolean zip_flag) throws IOException {
		return zip_flag ? new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(new File(file_name + ".gz"), append)))) :
			new BufferedWriter(new FileWriter(new File(file_name), append));
	}
	
	private boolean isCirc(String chr, int start, int end) {
		if (circ_table == null || circ_table.get(chr) == null) {
			return false;
		}
		Iterator<Node<CircleClip>> circ_nodes = circ_table.get(chr).overlappers(start, end);
		while (circ_nodes.hasNext()) {
			Node<CircleClip> node = circ_nodes.next();
			if ((start <= node.getStart() && end >= node.getStart())
					|| (start <= node.getEnd() && end >= node.getEnd())) {
				return true;
			}
			
		}
		return false;
	}
	
	private void writeMutExonSeq() throws IOException {
//		List<String> out = new ArrayList<>();
//		muts.forEach((chr, mutTree) -> {
//			for (Iterator<Node<Mutation>> mutNodes = mutTree.overlappers(0, Integer.MAX_VALUE); mutNodes.hasNext();) {
//				Node<Mutation> mutNode = mutNodes.next();
//				Transcript script = mutNode.getValue().getScript();
//				int start = Math.max(script.getStart(), mutNode.getStart() - args.getAlignmentLength() - 1);
//				int end = Math.min(script.getEnd(), mutNode.getEnd() + args.getAlignmentLength());
//				int start = script.inExons(mutNode.getEnd());
//				int end = start < 0 ? Math.min(mutNode.getEnd() + args.getAlignmentLength(), script.getEnd()) : 
//					script.getExonSeqSite(Math.min(start + args.getAlignmentLength(), script.getExonLength() - 1));
//				start = start < 0 ? Math.max(mutNode.getEnd() - args.getAlignmentLength() - 1, script.getStart()) :
//					script.getExonSeqSite(Math.max(start - args.getAlignmentLength() - 1, 0));
//				StringBuilder sb = new StringBuilder();
//				sb.append(chr);
//				sb.append('\t');
//				sb.append(start);
//				sb.append('\t');
//				sb.append(end);
//				sb.append('\t');
//				sb.append(genome.getChr(chr).getSeq().substring(start, end));
//				out.add(sb.toString());
//			}
//		});
//		Method.writeFile(args.getMutFile() + ".seq.txt", out, null);
		List<String> structList = new ArrayList<>();
		m.scaf_table.forEach((chr, scafTree) -> {
			for (Iterator<Node<Scaffold[]>> scafNodes = scafTree.overlappers(0, Integer.MAX_VALUE); scafNodes.hasNext();) {
				Node<Scaffold[]> scafNode = scafNodes.next();
				StringBuilder sb = new StringBuilder();
				sb.append(chr);
				sb.append('\t');
				sb.append(scafNode.getStart() - 1);
				sb.append('\t');
				sb.append(scafNode.getEnd());
				sb.append('\t');
				sb.append(scafNode.getValue()[0].getSitePair());
				sb.append('\t');
				sb.append(scafNode.getValue()[1].getSitePair());
				sb.append('\t');
				sb.append(CommonMethod.getHammingDistance(scafNode.getValue()[0].getSitePair(), scafNode.getValue()[1].getSitePair()));
				structList.add(sb.toString());
			}
		});
		Method.writeFile(m.args.getOutPrefix() + "_allele_struct.bed", structList, null);
	}
	
	public static void ensureMutDes(Map<String, IntervalTree<Peak>> peak_table) {
		for (Entry<String, IntervalTree<Mutation>> entry : m.muts.entrySet()) {
			for (Node<Mutation> mut_node : entry.getValue()) {
				Mutation mut = mut_node.getValue();
				if (peak_table.containsKey(entry.getKey())) {
					Iterator<Node<Peak>> peak_nodes = peak_table.get(entry.getKey()).overlappers(
							mut_node.getEnd() - m.args.getReadLength(), mut_node.getEnd() + m.args.getReadLength());
					if (peak_nodes.hasNext()) {
						if (!"Peak".equals(mut.getDescription())) {
							mut.setDescription("Peak");
						}
					}
					else if ("Peak".equals(mut.getDescription())){
						mut.setDescription("None");
					}
				}
				else if ("Peak".equals(mut.getDescription())){
					mut.setDescription("None");
				}
			}
		}
	}
	
	private String reverseFliq(String bases, int k) {
		return (k & 1) == 0 ? bases : reverseFliq(bases);
	}
	
	public static String reverseFliq(String bases) {
		if (bases == null || bases.length() == 0) {
			return bases;
		}
		StringBuilder sb = new StringBuilder();
		for (int i = bases.length() - 1; i >= 0; i--) {
			char c = BASE_MAP.containsKey(bases.charAt(i)) ? BASE_MAP.get(bases.charAt(i)) : bases.charAt(i);
			sb.append(c);
		}
		return sb.toString();
	}
	
//	public static String getMutBases(String bases, ArrayList<Node<Mutation>> muts, int start) {
//		StringBuffer sb = new StringBuffer();
//		sb.append(bases);
//		for (int i = 0; i < muts.size(); i++) {
//			sb.setCharAt(muts.get(i).getStart() - 1 - start, muts.get(i).getValue().getSNP()[0]);
//		}
//		return sb.toString();
//	}
//	
//	public static String getMutBases(String bases, ArrayList<Node<Mutation>> muts, int start, ArrayList<Integer> mut_points, int point_start) {
//		if (bases == null || muts == null) {
//			return bases;
//		}
//		StringBuffer sb = new StringBuffer();
//		sb.append(bases);
//		int index = binarySearch(start, muts);
//		for (int i = index; i < muts.size() && muts.get(i).getStart() < start + bases.length(); i++) {
//			sb.setCharAt(muts.get(i).getStart() - start, muts.get(i).getSNP()[0]);
//			if (mut_points != null) {
//				mut_points.add(muts.get(i).getStart() - start + point_start);
//			}
//		}
//		return sb.toString();
//	}
	
	private String buildUniqID(String chr, String annote) {
		StringBuffer out = new StringBuffer("@");
		String[] cols = annote.split("[\t\\-]");
		out.append(cols[0]);
		out.append(':');
		out.append(cols[1]);
		out.append(':');
		out.append(cols[2]);
		out.append(':');
		out.append(cols[cols.length - 2]);
		out.append(':');
		out.append(annote.substring(annote.length() - 8, annote.length()));
		out.append(':');
		out.append(getUniqueID());
//		int end = Integer.parseInt(cols[cols.length - 2]);
//		if (end >= 32229920 && end <= 32230068) {
//			++counts[end - 32229920];
//		}
		return out.toString();
	}
	
	private synchronized int getUniqueID() {
		return ++id_count;
	}
	
	public static int writeWithAnnote(BufferedWriter writer1, BufferedWriter writer2, BufferedWriter writer, List<Fastq> list){
		try {
			synchronized (SimMethod.class) {
				if (m.args.isPairEnd()) {
					writeInTwo(writer1, writer2, list);
				}
				else {
					write(writer1, list);
				}
				if (writer != null) {
					writeAnnote(writer, list);
				}
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
		return m.args.isPairEnd() ? list.size() / 2 : list.size();
	}
	
	private static void writeAnnote(BufferedWriter writer, List<Fastq> list) throws IOException {
		for (int i = 0; i + 1 < list.size(); i += m.args.isPairEnd() ? 2 : 1){
			writer.write(list.get(i).getAnnotation());
			writer.newLine();
		}
		writer.flush();
	}
	
	public static <T> void writeInTwo(BufferedWriter writer1, BufferedWriter writer2, List<T> list) throws IOException {
		for (int i = 0; i + 1 < list.size(); ++i){
			writer1.write(list.get(i).toString());
			writer1.newLine();
			writer2.write(list.get(++i).toString());
			writer2.newLine();
		}
		writer1.flush();
		writer2.flush();
	}
	
	private static <T> void write(BufferedWriter writer, List<T> list) throws IOException {
		for (int i = 0; i + 1 < list.size(); ++i){
			writer.write(list.get(i).toString());
			writer.newLine();
		}
		writer.flush();
	}
	
	public static Map<String, IntervalTree<Mutation>> getMuts() {
		return m.muts;
	}
	
	public static HashSet<Integer> randomNoRepeat(int start, int end, int num) {
		if (num > end - start) {
			return null;
		}
		HashSet<Integer> out = new HashSet<>();
		UniformIntegerDistribution uid = new UniformIntegerDistribution(start, end);
		while (out.size() < num) {
			out.add(uid.sample());
		}
		return out;
	}
	
//	public static int binarySearch(int target, ArrayList<? extends IntRegion> list, int start, int end) {
//		if (start >= end) {
//			return start;
//		}
//		int mid = (start + end) >> 1;
//		return target > list.get(mid).getStart() ? binarySearch(target, list, mid + 1, end) : binarySearch(target, list, start, mid);
//	}

	public static int binarySearch(int target, ArrayList<? extends IntRegion> list) {
		int index = Collections.binarySearch(list, new IntRegion(target, target + 1));
		return index < 0 ? -index - 1 : index;
	}
	
	public static <T extends IntRegion> ArrayList<IntRegion> mergerRegion(ArrayList<T> list){
		if (list == null || list.size() == 0) {
			return null;
		}
		ArrayList<IntRegion> out = new ArrayList<>();
		Collections.sort(list);
		out.add(new IntRegion(list.get(0).getStart(), list.get(0).getEnd()));
		IntRegion last_region = out.get(0);
		for (int i = 1; i < list.size(); i++) {
			IntRegion ir = list.get(i);
			if (ir.getStart() <= last_region.getEnd()) {
				if (ir.getEnd() > last_region.getEnd()) {
					last_region.resetStartAndEnd(last_region.getStart(), ir.getEnd());
				}
			}
			else {
				out.add(last_region = new IntRegion(ir.getStart(), ir.getEnd()));
			}
		}
		return out;
	}
	
	public static <T extends IntRegion> ArrayList<Integer> getRegionSizes(ArrayList<T> list){
		if (list == null || list.size() == 0) {
			return null;
		}
		ArrayList<Integer> out = new ArrayList<>();
		Collections.sort(list);
		for (int i = 0; i < list.size(); i++) {
			out.add(list.get(i).getLength());
		}
		return out;
	}
	
	public static String buildQualiString(int length, String chars) {
		if (chars == null || chars.length() < 1) {
			return null;
		}
		if (chars.length() == 1) {
			if (m.FixMapQ != null && length <= m.FixMapQ.length()){
				return m.FixMapQ.substring(0, length);
			}
			StringBuilder sb = new StringBuilder();
			while (--length >= 0) {
				sb.append(chars);
			}
			return m.FixMapQ = sb.toString();
		}
		UniformIntegerDistribution uid = new UniformIntegerDistribution(0, chars.length() - 1);
		StringBuffer sb = new StringBuffer();
		while(--length >= 0) {
			sb.append(chars.charAt(uid.sample()));
		}
		return sb.toString();
	}
	
	public static int sampleReadCount(double p) {
		return p <= 0.0 ? 0 : new PoissonDistribution(p).sample();
	}
	
	public static <T extends AbstractIntegerDistribution> int sampleReadLength(T aid, int min_len, int max_len) {
		int out = aid.sample();
		while (out < min_len) {
			out = aid.sample();
		}
		out = out > max_len ? max_len : out;
		return out;
	}
	
	public static <T extends AbstractIntegerDistribution> int sampleCutInRange(T aid, int lower, int upper) {
		int s = aid.sample();
		while (s < lower || s > upper) {
			s = aid.sample();
		}
		return s;
	}
	
	public static <T extends AbstractIntegerDistribution>  int sampleInRange(T aid, int lower, int upper) {
		int length = upper - lower + 1;
		int multi = (aid.getSupportUpperBound() - aid.getSupportLowerBound() + 1) / length;	
		if (multi == 0) {
			System.err.println("Cannot sample!");
			return 0;
		}
		int s = aid.sample();
		while (s > length * multi) {
			s = aid.sample();
		}
		s = lower + (s - aid.getSupportLowerBound()) % length;
		return s;
	}
	
	public static double sampleMutPer(double mut_per) {
		double s = mut_per_dis.sample();
		while(s < - mut_per || s > 1.0 - mut_per) {
			s = mut_per_dis.sample();
		}
		return s + mut_per;
	}
	
	public static double sampleMutPerOR(double mut_per) {
		return 1.0 - 1.0 / ( 1.0 + sampleORValue() * mut_per / (1.0 - mut_per));
	}
	
	private static double sampleORValue() {
		return Math.exp(Math.floor((OR_dis.sample() + 0.15) / 0.3) * 0.3);
	}
	
	public static <T extends AbstractRealDistribution> double sampleCutInRange(T ard, double lower, double upper) {
		double s = ard.sample();
		while (s < lower || s > upper) {
			s = ard.sample();
		}
		return s;
	}
	
	public static <T extends IntRegion> T moreNear(IntRegion target, T t1, T t2) {
		return Math.abs(t2.getStart() - target.getStart()) + Math.abs(t2.getEnd() - target.getEnd()) < 
				Math.abs(t1.getStart() - target.getStart()) + Math.abs(t1.getEnd() - target.getEnd()) ?
						t2 : t1;
	}
	
	public static boolean withIn(int target, int lower, int upper) {
		if (lower > upper) {
			lower ^= upper;
			upper ^= lower;
			lower ^= upper;
		}
		return target >= lower && target <= upper;
	}
	
	public static String toString(Object obj) {
		if (obj == null) {
			return null;
		}
		StringBuilder sb = new StringBuilder();
		HashSet<Object> visited = new HashSet<>();
		visited.add(null);
		toString(sb, obj, visited);
		return sb.toString();
	}
	
	private static void toString(StringBuilder sb, Object obj, HashSet<Object> visited) {
		if (obj == null) {
			sb.append("null");
		}
		if (visited.add(obj)) {
			Class<?> class1 = obj.getClass();
			if (class1 == String.class) {
				sb.append(obj);
			}
			else if (class1.isArray()) {
				sb.append(class1.getComponentType());
				sb.append("[]{");
				if (Array.getLength(obj) > 0) {
					if (class1.getComponentType().isPrimitive()) {
						sb.append(Array.get(obj, 0));
					}
					else {
						toString(sb, Array.get(obj, 0), visited);
					}
				}
				for (int i = 1; i < Array.getLength(obj); i++) {
					sb.append(',');
					if (class1.getComponentType().isPrimitive()) {
						sb.append(Array.get(obj, i));
					}
					else {
						toString(sb, Array.get(obj, i), visited);
					}
				}
				sb.append("}\n");
			}
			else {
				sb.append(class1.getName());
				do {
					sb.append('[');
					Field[] fields = class1.getDeclaredFields();
					AccessibleObject.setAccessible(fields, true);
					for (Field field : fields) {
						if (!Modifier.isStatic(field.getModifiers())) {
							if (sb.charAt(sb.length() - 1) != '[') {
								sb.append(',');
							}
							sb.append(field.getName());
							sb.append('=');
							try {
								Object val = field.get(obj);
								if (field.getType().isPrimitive()) {
									sb.append(val);
								}
								else {
									toString(sb, val, visited);
								}
							} catch (Exception e) {
								e.printStackTrace();
							}
							
						}
					}
					sb.append(']');
					class1 = class1.getSuperclass();
				}while(class1 != null);
			}
		}
		else {
			sb.append("...");
		}
	}
}
