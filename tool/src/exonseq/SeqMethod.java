package exonseq;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import bed.*;
import genome.*;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import main.InParam;
import main.Method;
import sim.SimMethod;

public class SeqMethod {
	private Genome genome = null;
	private static SeqMethod m = new SeqMethod();
	private Map<String, IntervalTree<Bed12>> bed_table = null;
	private Map<String, IntervalTree<Mutation>> mut_table = null;
	private Map<Bed12, HashMap<Transcript, Integer>> script_ponits = null;
	
	public static void run() {
		try {
			m.genome = Method.loadGenomeInfo(InParam.getParams().getExonFile(), InParam.getParams().getGenomeFile());
			Method.printNow("Genome complete!");
			m.bed_table = InParam.getParams().getPeakFile() == null ? null : Method.getInstance().readBedFile(InParam.getParams().getPeakFile(), true, null);
			Method.printNow("Load complete!");
			if (InParam.getParams().getMutFile() != null) {
				m.mut_table = Method.readMutation(InParam.getParams().getMutFile(), null);
				m.getSeq();
				return;
			}
			m.script_ponits = Method.getInstance().seqBamScript(m.genome, InParam.getParams().getInputSams(), m.bed_table);
			
			m.annote();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void getSeq() throws IOException {
		List<String> outList = new ArrayList<>();
		mut_table.forEach((chr, bedTree) -> {
			if (genome.getChr(chr) == null || genome.getChr(chr).getSeq() == null) {
				return;
			}
			for (Iterator<Node<Mutation>> bedNodes = bedTree.overlappers(0, Integer.MAX_VALUE); bedNodes.hasNext();) {
				Node<Mutation> bedNode = bedNodes.next();
				int start = Math.max(0, bedNode.getStart() - InParam.getParams().getAlignmentLength() - 1);
				int end = Math.min(genome.getChr(chr).getSeq().length(), bedNode.getEnd() + InParam.getParams().getAlignmentLength());
				StringBuilder sb = new StringBuilder();
				sb.append(chr);
				sb.append('\t');
				sb.append(start);
				sb.append('\t');
				sb.append(end);
				sb.append('\t');
				sb.append(genome.getChr(chr).getSeq().substring(start, end).toUpperCase());
				outList.add(sb.toString());
			}
		});
		Method.writeFile(InParam.getParams().getOutPrefix() + "_seq.bed", outList, null);
	}

	private void annote() throws IOException {
		final List<Bed12> annoted_bed = new ArrayList<>();
		m.bed_table.forEach((chr_string, bed_tree) -> {
			Chromosome chr = genome.getChr(chr_string);
			if (chr != null && chr.getSeq() != null) {
				for (Node<Bed12> bed_node : bed_tree) {
					annoted_bed.addAll(getAnnotedBeds(chr, bed_node));
				}
			}
			else {
				System.out.println("Warning: No " + chr_string + " in gtf file or genome.fa file");
			}
		});
		Method.writeFile(InParam.getParams().getOutPrefix() + "_annoted.bed", annoted_bed, null);
	}
	
	private List<Bed12> getAnnotedBeds(Chromosome chr, Node<Bed12> bed_node) {
		List<Bed12> out = new ArrayList<>();
		HashMap<String, Bed12> front_map = new HashMap<>();
		HashMap<String, Bed12> behind_map = new HashMap<>();
		for (Transcript script : locateScript(bed_node, chr)) {
			if (ensureBases(script, bed_node.getValue().getStart(), bed_node.getValue().getExamine())) {
				int seq_start = Math.max(getSeqIndex(script, bed_node.getValue().getStart()) - InParam.getParams().getAlignmentLength(), 0);
				int bed_start = getSeqIndex(script, bed_node.getValue().getStart());
				int bed_end = getSeqIndex(script, bed_node.getValue().getEnd());
				int seq_end = Math.min(getSeqIndex(script, bed_node.getValue().getEnd()) + InParam.getParams().getAlignmentLength(), script.getExonLength());
				String front_seq = script.getExonBases().substring(seq_start, bed_start);
				String behind_seq = script.getExonBases().substring(bed_end, seq_end);
				if (script.getGene().getStrand() == '-') {
					front_seq = SimMethod.reverseFliq(front_seq);
					behind_seq = SimMethod.reverseFliq(behind_seq);
				}
				if (front_map.containsKey(front_seq)) {
					Bed12 bed = front_map.get(front_seq);
					if (bed != null) {
						bed.setName(script.getScript_id() + ":" + seq_start + "|" + bed.getName());
					}
				}
				else {
					List<List<Integer>> blocks = getBlocks(script, seq_start, bed_start);
					int start = getGenomeIndex(script, seq_start);
					int end = getGenomeIndex(script, bed_start - 1) + 1;
					if (overlapRegion(m.bed_table.get(chr.getChr()), blocks, start)) {
						front_map.put(front_seq, null);
					}
					else {
						front_map.put(front_seq, new Bed12(bed_node.getValue().getStart(), bed_node.getValue().getEnd(), script.getChr(), 
								script.getScript_id() + ":" + seq_start + "|UpStream," + front_seq, bed_start - seq_start,
								script.getGene().getStrand(), start, end, script.getGene().getStrand() == '-' ? -1 : 1, 
								blocks.get(0).size(), blocks.get(0), blocks.get(1)));
						out.add(front_map.get(front_seq));
						front_map.get(front_seq).addInfo(bed_node.getValue().getDescription());
					}
				}
				if (behind_map.containsKey(behind_seq)) {
					Bed12 bed = behind_map.get(behind_seq);
					if (bed != null) {
						bed.setName(script.getScript_id() + ":" + seq_start + "|" + bed.getName());
					}
				}
				else {
					List<List<Integer>> blocks = getBlocks(script, bed_end, seq_end);
					int start = getGenomeIndex(script, bed_end);
					int end = getGenomeIndex(script, seq_end - 1) + 1;
					if (overlapRegion(m.bed_table.get(chr.getChr()), blocks, start)){
						behind_map.put(behind_seq, null);
					}
					else {
						behind_map.put(behind_seq, new Bed12(bed_node.getValue().getStart(), bed_node.getValue().getEnd(), script.getChr(), 
								script.getScript_id() + ":" + seq_start + "|DownStream," + behind_seq, seq_end - bed_end,
								script.getGene().getStrand(), start, end, script.getGene().getStrand() == '-' ? -1 : 1,
										blocks.get(0).size(), blocks.get(0), blocks.get(1)));
						out.add(behind_map.get(behind_seq));
						behind_map.get(behind_seq).addInfo(bed_node.getValue().getDescription());
					}
				}
			}
		}
		return out;
	}
	
	private List<Transcript> locateScript(Node<Bed12> bed_node, Chromosome chr) {
		List<Transcript> scripts = new ArrayList<>();
		int max_point = -1;
		HashMap<Integer, HashSet<Transcript>> start_scripts = new HashMap<>();
		for (Iterator<Node<Exon>> exon_nodes = chr.getExonTree().overlappers(bed_node.getStart(), bed_node.getStart());
				exon_nodes.hasNext();) {
			Node<Exon> exon_node = exon_nodes.next();
			for (Exon exon = exon_node.getValue(); exon != null; exon = exon.getAnotherExon()) {
				int point = script_ponits == null ? 0 : 
					script_ponits.getOrDefault(bed_node.getValue(), new HashMap<>())
					.getOrDefault(exon.getScript(), 0);
				if (point > max_point) {
					max_point = point;
					start_scripts.put(max_point, new HashSet<>());
				}
				if (point == max_point) {
					start_scripts.get(max_point).add(exon.getScript());
				}
			}
		}
		if (max_point >= 0) {
			if (bed_node.getStart() == bed_node.getEnd()) {
				scripts.addAll(start_scripts.get(max_point));
			}
			else {
				for (Iterator<Node<Exon>> exon_nodes = chr.getExonTree().overlappers(bed_node.getEnd(), bed_node.getEnd());
						exon_nodes.hasNext();) {
					Node<Exon> exon_node = exon_nodes.next();
					for (Exon exon = exon_node.getValue(); exon != null; exon = exon.getAnotherExon()) {
						if (start_scripts.get(max_point).contains(exon.getScript())) {
							scripts.add(exon.getScript());
						}
					}
				}
			}
		}
		return scripts;
	}
	
//	private Bed12 getAnnotedBed(Transcript script, Node<Bed12> bed_node) {
//		int seq_start = Math.max(getSeqIndex(script, bed_node.getValue().getStart()) - InParam.getParams().getAlignmentLength(), 0);
//		int seq_end = Math.min(getSeqIndex(script, bed_node.getValue().getEnd()) + InParam.getParams().getAlignmentLength(), script.getExonLength());
//		List<List<Integer>> blocks = getBlocks(script, seq_start, seq_end);
//		int start = getGenomeIndex(script, seq_start);
//		int end = getGenomeIndex(script, seq_end - 1);
//		return new Bed12(bed_node.getValue().getStart(), bed_node.getValue().getEnd(), script.getChr(), 
//				script.getScript_id() + ":" + script.getExonBases(genome.getChr(script.getChr()).getSeq()).substring(seq_start, seq_end), seq_end - seq_start,
//				script.getGene().getStrand(), start, end, 0, blocks.get(0).size(), blocks.get(0), blocks.get(1));
//	}
	
	private boolean ensureBases(Transcript script, int genome_index, String bases) {
		if (bases == null || bases.length() == 0) {
			return true;
		}
//		if (script.getGene().getStrand() == '-') {
//			return genome.getChr(script.getChr()) != null && 
//					genome_index + bases.length() <= genome.getChr(script.getChr()).getSeq().length() &&
//					genome.getChr(script.getChr()).getSeq().substring(genome_index, genome_index + bases.length()).equalsIgnoreCase(
//						SimMethod.reverseFliq(bases));
//		}
//		else {
			return genome.getChr(script.getChr()) != null && 
					genome_index + bases.length() <= genome.getChr(script.getChr()).getSeq().length() &&
					genome.getChr(script.getChr()).getSeq().substring(genome_index, genome_index + bases.length()).equalsIgnoreCase(bases);
//		}
	}
	
	public static int getSeqIndex(Transcript script, int genome_index) {
		int seq_index = 0;
		for (Exon exon = script.getExons(); exon != null; exon = exon.getDownStream()) {
			if (genome_index >= exon.getStart() && genome_index <= exon.getEnd()) {
				return seq_index + genome_index - exon.getStart();
			}
			seq_index += exon.getLength();
		}
		return -1;
	}
	
	public static int getGenomeIndex(Transcript script, int seq_index) {
		if (seq_index < 0) {
			return -1;
		}
		for (Exon exon = script.getExons(); exon != null; exon = exon.getDownStream()) {
			if (seq_index < exon.getLength()) {
				return seq_index + exon.getStart();
			}
			seq_index -= exon.getLength();
		}
		return -1;
	}
	
	private List<List<Integer>> getBlocks(Transcript script, int seq_start, int seq_end) {
		if (script == null) {
			return null;
		}
		List<List<Integer>> blocks = new ArrayList<>();
		List<Integer> block_starts = new ArrayList<>();
		List<Integer> block_sizes = new ArrayList<>();
		seq_start = Math.max(0, seq_start);
		seq_end = Math.min(script.getExonLength(), seq_end);
		int start = 0;
		for (Exon exon = script.getExons(); exon != null; exon = exon.getDownStream()) {
			if (seq_start < exon.getLength()) {
				if (seq_start >= 0) {
					block_starts.add(0);
					block_sizes.add(Math.min(seq_end, exon.getEnd() - exon.getStart()) - seq_start);
					start = exon.getStart() + seq_start;
				}
				else {
					block_starts.add(exon.getStart() - start);
					block_sizes.add(Math.min(seq_end, exon.getEnd() - exon.getStart()));
				}
			}
			seq_start -= exon.getLength();
			seq_end -= exon.getLength();
			if (seq_end <= 0) {
				break;
			}
		}
		blocks.add(block_sizes);
		blocks.add(block_starts);
		return blocks;
	}
	
	private boolean overlapRegion(IntervalTree<?> tree, List<List<Integer>> blocks, int start) {
		if (tree == null || blocks == null || blocks.size() < 2 || blocks.get(0).size() != blocks.get(1).size()) {
			return false;
		}
		for (int i = 0; i < blocks.get(0).size(); ++i) {
			if (tree.overlappers(start + blocks.get(1).get(i) + 1, start + blocks.get(0).get(i)).hasNext()) {
				return true;
			}
		}
		return false;
	}
}
