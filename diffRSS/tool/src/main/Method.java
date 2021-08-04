package main;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Logger;
import java.util.Set;
import java.util.spi.LocaleNameProvider;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.distribution.UniformIntegerDistribution;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import bed.*;
import genome.*;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import smile.stat.hypothesis.KSTest;
import mapping.MappingMethod;
import mapping.MappingStat;
import sim.Scaffolds;

public class Method {

	private static Method m = new Method();
	private boolean ip = false;
	private double ref_fre = 0.0;
	private Map<String, IntervalTree<Mutation>> mut_table = null;

	
	private Method() {
	}
	
	public static Method getInstance() {
		return m;
	}
	
	public static void stat() {
		if (InParam.getParams().getInputSams().size() < 0) {
			return;
		}
		String ref_file = InParam.getParams().getControlFile();
		StringBuffer header = new StringBuffer();
		try {
			HashMap<String, IntervalTree<Bed3Extend>> ref_tree = m.readBed3(ref_file, header);
			for (int i = 1; i < InParam.getParams().getTreatFiles().size(); ++i) {
				m.statSup(InParam.getParams().getTreatFiles().get(i), ref_tree);
				String tool = InParam.getParams().getTreatFiles().get(i).substring(
						InParam.getParams().getTreatFiles().get(i).lastIndexOf('/') + 1);
				tool = tool.substring(0, tool.indexOf('_') < 0 ? tool.length() : tool.indexOf('_'));
				header.append('\t');
				header.append(tool);
			}
			writeFile(InParam.getParams().getOutPrefix() + "_sup.txt", ref_tree, false, false, header.toString());
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void statSup(String file, HashMap<String, IntervalTree<Bed3Extend>> ref_table) throws IOException {
		HashMap<String, IntervalTree<Bed3Extend>> f_table = readBed3(file, null);
		final int dev = 3;
		ref_table.forEach((chr, ref_tree) -> {
			IntervalTree<Bed3Extend> f_tree = f_table.getOrDefault(chr, new IntervalTree<>());
			for (Node<Bed3Extend> ref_node : ref_tree) {
				boolean match = false;
				for (Iterator<Node<Bed3Extend>> f_nodes = f_tree.overlappers(ref_node.getStart() - dev, ref_node.getEnd() + dev);
						f_nodes.hasNext();) {
					Node<Bed3Extend> f_node = f_nodes.next();
					if (Math.abs(f_node.getStart() - ref_node.getStart()) <= dev &&
							Math.abs(f_node.getEnd() - ref_node.getEnd()) <= dev) {
						match = true;
						ref_node.getValue().appDes("\t" + f_node.getValue().getDescription());
						break;
					}
				}
				if (!match) {
					ref_node.getValue().appDes("\t0");
				}
			}
		});
	}

	public static void cover() {
		try {
//			CountMethod.countSNP();
			Genome genome = loadGenomeInfo(InParam.getParams().getExonFile(), null);
			m.mut_table = readMutation(InParam.getParams().getMutFile(), null);
			m.ip = false;
			HashMap<String, IntervalTree<MappingStat>> input_table = m.readBamFile(InParam.getParams().getControlFile());
			m.ip = true;
			HashMap<String, IntervalTree<MappingStat>> ip_table = m.readBamFile(InParam.getParams().getTreatFile());
			m.setPeakStrand(genome, InParam.getParams().getPeakFile());
			m.calESES(genome, input_table, ip_table);
			m.writeGenes(genome, InParam.getParams().getOutPrefix() + "_gene_score.bed");
			if (m.mut_table.size() > 0) {
				writeFile(InParam.getParams().getOutPrefix() + "_snp.bed", m.mut_table, true, true, Mutation.getHeader());
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void countSNP(SAMRecord record) {
		if (mut_table != null && mut_table.size() > 0) {
			IntervalTree<Mutation> mut_tree = mut_table.get(record.getReferenceName());
			if (mut_tree != null) {
				IntervalTree<Segment> seg_tree = BamMethod.cigarToSeq(record);
				for (Node<Segment> seg_node : seg_tree) {
					if (seg_node.getValue().getSeq() == null) {
						continue;
					}
					for (Iterator<Node<Mutation>> mut_nodes = mut_tree.overlappers(seg_node.getStart(), seg_node.getEnd());
							mut_nodes.hasNext();) {
						Node<Mutation> mut_node = mut_nodes.next();
						mut_node.getValue().incCount(seg_node.getValue().getSeq().charAt(mut_node.getStart() - seg_node.getStart()), ip);
					}
				}
			}
		}
	}
	
	private void calESES(Genome genome, HashMap<String, IntervalTree<MappingStat>> input_table,
			HashMap<String, IntervalTree<MappingStat>> ip_table) {
		final int window_size = 25;
		ArrayList<Double> input_counts = new ArrayList<>();
		ArrayList<Double> ip_counts = new ArrayList<>();
		final int[] count_index = {0};
		genome.getChrMap().entrySet().forEach(entry -> {
			IntervalTree<Gene> gene_tree = entry.getValue().getGeneTree();
			IntervalTree<MappingStat> input_tree = input_table.getOrDefault(entry.getKey(), new IntervalTree<>());
			IntervalTree<MappingStat> ip_tree = ip_table.getOrDefault(entry.getKey(), new IntervalTree<>());
			int index = 1;
			for (Gene gene = nextESESGene(gene_tree, input_tree, ip_tree, index); gene != null;
					gene = nextESESGene(gene_tree, input_tree, ip_tree, index)) {
				int len = 0;
				for (int i = gene.getStart() + 1; i <= gene.getEnd(); i += window_size) {
					input_counts.add((double) countCover(input_tree, i, i + window_size - 1));
					ip_counts.add((double) countCover(ip_tree, i, i + window_size - 1));
					++len;
				}
				int count = 0;
				for (Iterator<Node<MappingStat>> map_nodes = input_tree.overlappers(gene.getStart() + 1, gene.getEnd());
						map_nodes.hasNext();) {
					Node<MappingStat> map_node = map_nodes.next();
					count += map_node.getValue().getReads();
				}
				gene.setInputReadCount(count);
				count = 0;
				for (Iterator<Node<MappingStat>> map_nodes = ip_tree.overlappers(gene.getStart() + 1, gene.getEnd());
						map_nodes.hasNext();) {
					Node<MappingStat> map_node = map_nodes.next();
					count += map_node.getValue().getReads();
				}
				gene.setIpReadCount(count);
				gene.setIp_enrich(calESES(input_counts.subList(count_index[0], count_index[0] + len),
						ip_counts.subList(count_index[0], count_index[0] + len)));
				index = gene.getEnd() + 1;
				count_index[0] += len;
			}
		});
		System.out.println("Genome ESES: " + calESES(input_counts, ip_counts));
	}
	
	private Gene nextESESGene(IntervalTree<Gene> gene_tree, IntervalTree<MappingStat> input_tree,
			IntervalTree<MappingStat> ip_tree, int index) {
		Gene gene = null;
		while (gene == null && index > 0) {
			Node<MappingStat> input_node = input_tree.min(index, index);
			Node<MappingStat> ip_node = ip_tree.min(index, index);
			if (input_node == null) {
				if (ip_node == null) {
					index = -1;
				}
				else {
					index = ip_node.getStart() + 1;
					gene = nextESESGene(gene_tree, ip_node);
				}
			}
			else if (ip_node == null) {
				index = input_node.getStart() + 1;
				gene = nextESESGene(gene_tree, input_node);
			}
			else {
				if (ip_node.getStart() < input_node.getStart()) {
					index = ip_node.getStart() + 1;
					gene = nextESESGene(gene_tree, ip_node);
				}
				else {
					index = input_node.getStart() + 1;
					gene = nextESESGene(gene_tree, input_node);
				}
			}
		}
		return gene;
	}
	
	private Gene nextESESGene(IntervalTree<Gene> gene_tree, Node<MappingStat> node) {
		Gene gene = null;
		for (Iterator<Node<Gene>> gene_nodes = gene_tree.overlappers(node.getStart(), node.getEnd());
				gene_nodes.hasNext();) {
			Gene tmp = gene_nodes.next().getValue();
			if (gene == null || tmp.getLength() > gene.getLength()) {
				gene = tmp;
			}
		}
		return gene;
	}
	
	private int countCover(IntervalTree<MappingStat> map_tree, int start, int end) {
		int count = 0;
		for (Iterator<Node<MappingStat>> map_nodes = map_tree.overlappers(start, end); map_nodes.hasNext();) {
			Node<MappingStat> map_node = map_nodes.next();
			count += map_node.getValue().getReads();
		}
		return count;
	}
	
	private double calESES(List<Double> input_list, List<Double> ip_list) {
		if (ip_list == null || input_list == null || input_list.size() != ip_list.size()) {
			return -1.0;
		}
		double input_sum = 0.0;
		double ip_sum = 0.0;
		double eses = 0.0;
		for (int i = 0; i < input_list.size(); i++) {
			input_sum += input_list.get(i);
			ip_sum += ip_list.get(i);
		}
		if (input_sum == 0.0) {
			return 1.0;
		}
		if (ip_sum == 0.0) {
			return 0.0;
		}
		for (int i = 0; i < input_list.size(); i++) {
			input_list.set(i, input_list.get(i) / input_sum);
			ip_list.set(i, ip_list.get(i) / ip_sum);
		}
		Collections.sort(input_list);
		Collections.sort(ip_list);
		for (int i = 1; i < input_list.size(); i++) {
			input_list.set(i, input_list.get(i) + input_list.get(i - 1));
			ip_list.set(i, ip_list.get(i) + ip_list.get(i - 1));
		}
		for (int i = 0; i < input_list.size(); i++) {
			eses = Math.max(eses, input_list.get(i) - ip_list.get(i));
		}
		return eses;
	}
	
	private void setPeakStrand(Genome genome, String peakFile) throws IOException {
		if (peakFile == null) {
			return;
		}
		Map<String, IntervalTree<Peak>> peak_map = readPeak(peakFile, null);
		genome.getChrMap().forEach((chr_string, chr) -> {
			if (!peak_map.containsKey(chr_string)) {
				for (Node<Gene> gene_node : chr.getGeneTree()) {
					gene_node.getValue().setStrand('F');
				}
			}
			else {
				for (Node<Gene> gene_node : chr.getGeneTree()) {
					if (peak_map.get(chr_string).overlappers(gene_node.getStart(), gene_node.getEnd()).hasNext()) {
						gene_node.getValue().setStrand('T');
					}
					else {
						gene_node.getValue().setStrand('F');
					}
				}
			}
		});
	}
	
	private void writeGenes(Genome genome, String out_file) throws IOException {
		try (BufferedWriter bw = new BufferedWriter(new FileWriter(new File(out_file)))){
			for (Entry<String, Chromosome> entry : genome.getChrMap().entrySet()) {
				for (Node<Gene> node : entry.getValue().getGeneTree()) {
					if (node.getValue().getInputReadCount() + node.getValue().getIpReadCount() > 0) {
						bw.write(node.getValue().toString());
						bw.newLine();
					}
				}
			}
		}
	}

	public static void test(String[] args) {
		try {
			m.findUnused(args[0], args[1], args[2]);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void findUnused(String bed, String exon, String out) throws IOException {
		HashMap<String, IntervalTree<Bed12>> bed_table = readBedFile(bed, true, null);
		Genome genome = new Genome();
		genome.readAnnoteFile(exon);
		HashMap<String, IntervalTree<Bed12>> out_table = new HashMap<>();
		genome.getChrMap().forEach((chr_string, chr) -> {
			IntervalTree<Bed12> bed_tree = bed_table.getOrDefault(chr_string, new IntervalTree<>());
			IntervalTree<Bed12> out_tree = new IntervalTree<>();
			out_table.put(chr_string, out_tree);
			int index = 0;
			for (Node<Bed12> bed_node : bed_tree) {
				int start = -1;
				if (bed_node.getValue().getStrand() == '+') {
					start = bed_node.getStart() - 2;
				}
				else if (bed_node.getValue().getStrand() == '-') {
					start = bed_node.getEnd();
				}
				if (start < 0) {
					continue;
				}
				for (Iterator<Node<Exon>> exon_nodes = chr.getExonTree().overlappers(start, start + 1); 
						exon_nodes.hasNext();) {
					Node<Exon> exon_node = exon_nodes.next();
					if (bed_node.getValue().getStrand() == '+') {
						if (exon_node.getValue().getEnd() == exon_node.getValue().getScript().getEnd()
								&& (exon_node.getValue().getEnd() == start || exon_node.getValue().getEnd() == start + 1)) {
							Bed12 record = bed_node.getValue().deepClone();
							out_tree.put(bed_node.getStart(), bed_node.getStart() + ++index, record);
							record.addInfo("\t" + getExonInfo(exon_node.getValue()));
							record.setName(exon_node.getValue().getGene().getGene_id());
						}
					}
					else {
						if (exon_node.getValue().getStart() == exon_node.getValue().getScript().getStart()
								&& (exon_node.getValue().getStart() == start || exon_node.getValue().getStart() == start + 1)) {
							Bed12 record = bed_node.getValue().deepClone();
							out_tree.put(bed_node.getStart(), bed_node.getStart() + ++index, record);
							record.addInfo("\t" + getExonInfo(exon_node.getValue()));
							record.setName(exon_node.getValue().getGene().getGene_id());
						}
					}
				}
			}
		});
		writeFile(out, out_table, false, false, null);
	}

	private String getExonInfo(Exon exon) {
		StringBuilder sb = new StringBuilder();
		sb.append(exon.getChr());
		sb.append('-');
		sb.append(exon.getStart() + 1);
		sb.append('-');
		sb.append(exon.getEnd());
		sb.append('|');
		sb.append(exon.getScript());
		return sb.toString();
	}
	
	public static void run() {
		try {
			Map<String, IntervalTree<Mutation>> mut_map = readMutation(InParam.getParams().getMutFile(), null);
//			if (InParam.getParams().isSamFormat()) {
//				m.splitSamFile(InParam.getParams().getInputSams().get(0), InParam.getParams().getOutPrefix(), mut_map);
//			}
//			else {
//				m.splitBamFile(InParam.getParams().getInputSams().get(0), InParam.getParams().getOutPrefix(), mut_map);
//			}
			Genome genome = loadGenomeInfo(InParam.getParams().getExonFile(), InParam.getParams().getGenomeFile());
			if (InParam.getParams().getControlFile() != null) {
				m.ip = false;
				m.countMutInSam(InParam.getParams().getControlFile(), mut_map);
//				m.ref_fre = m.calRefFreq(mut_map);
				m.calMAF(mut_map, genome);
//				HashMap<String, IntervalTree<Gene>> out_file = new HashMap<>();
//				for (Entry<String, Chromosome> entry : genome.getChrMap().entrySet()) {
//					out_file.put(entry.getKey(), entry.getValue().getGeneTree());
//				}
//				writeFile(InParam.getParams().getOutPrefix() + "_gene_maf.txt", out_file, true, false, null);
				for (Entry<String, IntervalTree<Mutation>> entry : mut_map.entrySet()) {
					for (Node<Mutation> node : entry.getValue()) {
						Mutation mut = node.getValue();
						if (mut.getInputCount(mut.getRef()) == 0 || mut.getMutCount(false) == 0
								|| node.getValue().getTotalCount(false) < 10) {
							entry.getValue().remove(node.getStart(), node.getEnd());
						}
					}
				}
			}
			
			m.ip = true;
//			if (InParam.getParams().isSamFormat()) {
//				m.splitSamFile(InParam.getParams().getInputSams().get(1), InParam.getParams().getOutPrefix(), mut_map);
//			}
//			else {
//				m.splitBamFile(InParam.getParams().getInputSams().get(1), InParam.getParams().getOutPrefix(), mut_map);
//			}
			m.countMutInSam(InParam.getParams().getTreatFile(), mut_map);
//			m.compareMAF(mut_map, genome);
			writeFile(InParam.getParams().getOutPrefix() + "_mut.txt", m.mutOutput(genome, mut_map), Mutation.getSimpleHeader());
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public List<String> mutOutput(Genome genome, Map<String, IntervalTree<Mutation>> mut_table) {
		List<String> out = new ArrayList<>();
		char sep = '\t';
		List<Double> p_values = new ArrayList<>();
		for (Entry<String, IntervalTree<Mutation>> entry : mut_table.entrySet()) {
			for (Iterator<Node<Mutation>> nodes = entry.getValue().overlappers(Integer.MIN_VALUE, Integer.MAX_VALUE);
					nodes.hasNext();) {
				IntervalTree.Node<Mutation> node = nodes.next();
				Mutation mut = node.getValue();
				mut.remakeMostAlt(false);
				p_values.add(InParam.getParams().getControlFile() == null ? 
					CommonMethod.calPvalue(0.5, mut.getIpCount(mut.getSNP()[0]), mut.getIpCount(mut.getRef())) :
					CommonMethod.calPvalue(mut.getInputCount(mut.getSNP()[0]), mut.getInputCount(mut.getRef()),
						mut.getIpCount(mut.getSNP()[0]), mut.getIpCount(mut.getRef())));
			}
		}
		List<Double> fdr_values = new ArrayList<>(p_values);
		CommonMethod.adjustPValue(fdr_values, "bh");
		int index = 0;
		for (Entry<String, IntervalTree<Mutation>> entry : mut_table.entrySet()) {
			String seq = genome != null && genome.getChr(entry.getKey()) != null && 
					genome.getChr(entry.getKey()).getSeq() != null ? genome.getChr(entry.getKey()).getSeq() : null;
			for (Iterator<Node<Mutation>> nodes = entry.getValue().overlappers(Integer.MIN_VALUE, Integer.MAX_VALUE);
				nodes.hasNext(); ++index) {
				IntervalTree.Node<Mutation> node = nodes.next();
				StringBuffer sb = new StringBuffer();
				sb.append(entry.getKey());
				sb.append(sep);
				sb.append(node.getStart() - 1);
				sb.append(sep);
				sb.append(node.getEnd());
				sb.append(sep);
				sb.append(node.getValue().toSimpleString(p_values.get(index), fdr_values.get(index)));
				Transcript script = node.getValue().getScript();
				if (script != null) {
					sb.append(sep);
					sb.append(node.getEnd() - node.getValue().getScript().getStart());
					sb.append(sep);
					sb.append(script.getRegionFeature(node.getEnd(), true, true));
					sb.append(sep);
					sb.append(script.getUpStreamString(seq, node.getEnd(), 50));
					sb.append(sep);
					sb.append(script.getDownStreamString(seq, node.getEnd(), 50));
					sb.append(sep);
				}
				else {
					sb.append("\t0\tintergenic\t");
					if (seq != null && node.getEnd() >= 51 && node.getEnd() + 50 <= seq.length()) {
						sb.append(seq.substring(node.getEnd() - 51, node.getEnd() - 1));
						sb.append(sep);
						sb.append(seq.substring(node.getEnd(),node.getEnd() + 50));
						sb.append(sep);
					}
					else {
						sb.append(".\t.\t");
					}
				}
				if (seq != null) {
					sb.append(String.format("%.3f", BamMethod.getGCContent(seq, node.getEnd(), 1000)));
				}
				else {
					sb.append("0.000");
				}
				out.add(sb.toString());
			}
		}
		return out;
	}
	
	private double calRefFreq(HashMap<String, IntervalTree<Mutation>> mut_map) {
		double ref_num = 0.0;
		double total_num = 0.0;
		for (Entry<String, IntervalTree<Mutation>> entry : mut_map.entrySet()) {
			for (Node<Mutation> mut_node : entry.getValue()) {
				ref_num += mut_node.getValue().getCount(mut_node.getValue().getRef(), false);
				total_num += mut_node.getValue().getTotalCount(false);
			}
		}
		return ref_num / total_num;
	}

	public void runCount(String bed, String sam) {
		try {
			HashMap<String, IntervalTree<Bed12>> bed_map = readBedFile(bed, false, null);
			HashMap<String, IntervalTree<MappingStat>> sam_map = readBamFile(sam);
			HashMap<String, IntervalTree<Pair<Integer>>> out = countRegionInSam(bed_map, sam_map);
			writeFile(bed.substring(bed.lastIndexOf("/") + 1) + "_in_" + sam.substring(sam.lastIndexOf("/") + 1) + ".txt", out, true, true, null);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	private void countMutInSam(String sam_file, Map<String, IntervalTree<Mutation>> mut_map) throws IOException {
		try (SamReader reader = SamReaderFactory.makeDefault().open(new File(sam_file))) {
			for (SAMRecord record : reader) {
				countMutInRecord(record, mut_map);
			}
		}
	}
	
	private void countMutInRecord(SAMRecord record, Map<String, IntervalTree<Mutation>> mut_map) {
		if (record.getReadUnmappedFlag() || 
			!mut_map.containsKey(record.getReferenceName()) ||
			!mut_map.get(record.getReferenceName()).overlappers(record.getAlignmentStart(), record.getAlignmentEnd()).hasNext()) {
			return;
		}
		
		IntervalTree<Segment> cigar_tree = BamMethod.cigarToSeq(record);
		countMutInSegs(mut_map.get(record.getReferenceName()), cigar_tree);
	}

	private void countMutInSegs(IntervalTree<Mutation> mut_tree, IntervalTree<Segment> cigar_tree) {
		for (Node<Segment> cigar_node : cigar_tree) {
			for (Iterator<Node<Mutation>> nodes = mut_tree.overlappers(cigar_node.getStart(), cigar_node.getEnd());
					nodes.hasNext();) {
				Node<Mutation> node = nodes.next();
				if (cigar_node.getValue().getLength() == 0) {
					node.getValue().incCount('D', ip);
				}
				else {
					node.getValue().incCount(cigar_node.getValue().getSeq().charAt(node.getEnd() - cigar_node.getStart()), ip);
				}
			}
		}
	}

	public void feedSamFile(String file_name, String file_out) {
		try (SamReader reader = SamReaderFactory.makeDefault().open(new File(file_name));
				SAMFileWriter writer = new SAMFileWriterFactory().makeSAMWriter(reader.getFileHeader(), true, new File(file_out))) {
			int[] num = {1000, 500, 200};
			int[] count = {0, 0};
			String id = null;
			HashSet<Integer> quali = new HashSet<>();
			ArrayList<SAMRecord> records = new ArrayList<>();
			for (SAMRecord record : reader) {
				quali.add(record.getMappingQuality());
				if (record.getMappingQuality() > 1) {
					writer.addAlignment(record);
				}
				int pair_index = record.getReadPairedFlag() && record.getSecondOfPairFlag() ? 1 : 0;
				if (record.getReadName().equals(id)) {
					records.add(record);
					++count[pair_index];
				}
				else {
					if (id != null) {
						int max = Math.max(count[0], count[1]);
						if (max <= 3 && --num[max - 1] >= 0) {
							for (SAMRecord r : records) {
								writer.addAlignment(r);
							}
						}
					}
					count[0] = 0;
					count[1] = 0;
					id = record.getReadName();
					records.clear();
					records.add(record);
					++count[pair_index];
				}
				if (num[0] < 0 && num[1] < 0 && num[2] < 0) {
//					break;
				}
			}
			for (int i : quali) {
				System.out.println(i);
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void splitSamFile(String file_name, String file_out, HashMap<String, IntervalTree<Mutation>> bed_map) throws IOException {
		try (SamReader reader = SamReaderFactory.makeDefault().open(new File(file_name));
			SAMFileWriter mut_writer = new SAMFileWriterFactory().makeSAMWriter(setUnsort(reader.getFileHeader()), true, new File(file_out + "_alt.sam"));
			SAMFileWriter ref_writer = new SAMFileWriterFactory().makeSAMWriter(setUnsort(reader.getFileHeader()), true, new File(file_out + "_ref.sam"))){
			for (SAMRecord record : reader) {
				splitSamRecord(record, mut_writer, ref_writer, bed_map);
			}
		}
	}
	
	public static SAMFileHeader setUnsort(SAMFileHeader header) {
		SAMFileHeader out = header.clone();
		out.setSortOrder(SortOrder.unsorted);
		return out;
	}
	
	public static boolean isCoordinate(SAMFileHeader header) {
		return SortOrder.coordinate.equals(header.getSortOrder());
	}
	
	public static List<String> getChrs(SAMFileHeader header) {
		List<String> chrs = new ArrayList<>();
		for (SAMSequenceRecord seq : header.getSequenceDictionary().getSequences()) {
			chrs.add(seq.getSequenceName().intern());
		}
		return chrs;
	}
	
	public void splitBamFile(String file_name, String file_out, HashMap<String, IntervalTree<Mutation>> bed_map) throws IOException {
		try (SamReader reader = SamReaderFactory.makeDefault().open(new File(file_name));
			SAMFileWriter mut_writer = new SAMFileWriterFactory().makeBAMWriter(setUnsort(reader.getFileHeader()), true, new File(file_out + "_alt.bam"));
			SAMFileWriter ref_writer = new SAMFileWriterFactory().makeBAMWriter(setUnsort(reader.getFileHeader()), true, new File(file_out + "_ref.bam"))){
			for (SAMRecord record : reader) {
				splitSamRecord(record, mut_writer, ref_writer, bed_map);
			}
		}
	}
	
	public void splitSamRecord(SAMRecord record, SAMFileWriter mut_writer, SAMFileWriter ref_writer, HashMap<String, IntervalTree<Mutation>> bed_map) {
		if (record.getReadUnmappedFlag()) {
			return;
		}
		if (!bed_map.containsKey(record.getReferenceName())) {
			writeWithoutSNP(ref_writer, mut_writer, record, -1.0);
			return;
		}
		Iterator<Node<Mutation>> nodes = bed_map.get(record.getReferenceName()).overlappers(record.getAlignmentStart(), record.getAlignmentEnd());
		if (!nodes.hasNext()) {
			writeWithoutSNP(ref_writer, mut_writer, record, -1.0);
			return;
		}
		
		IntervalTree<Segment> cigar_tree = BamMethod.cigarToSeq(record);
		splitWrite(record, mut_writer, ref_writer, bed_map.get(record.getReferenceName()), cigar_tree);
	}
	
	private void splitWrite(SAMRecord record, SAMFileWriter mut_writer, SAMFileWriter ref_writer, IntervalTree<Mutation> mut_tree, IntervalTree<Segment> cigar_tree) {
		int align_start = record.getAlignmentStart();
		int cigar_start = 0;
		int cigar_end = 0;
		
		int last_mut = 0;
		int last_ref = 0;
		boolean last_cigar_mut = false;
		boolean last_cigar_ref = false;
		Iterator<Node<Segment>> cigar_nodes = cigar_tree.overlappers(Integer.MIN_VALUE, Integer.MAX_VALUE);
		while (cigar_nodes.hasNext()) {
			Node<Segment> cigar_node = cigar_nodes.next();
			Iterator<Node<Mutation>> nodes = mut_tree.overlappers(cigar_node.getStart(), cigar_node.getEnd());
			while (nodes.hasNext()) {
				Node<Mutation> node = nodes.next();
				node.getValue().incCount(cigar_node.getValue().getSeq().charAt(node.getEnd() - cigar_node.getStart()), ip);
				if (node.getValue().getRef() != cigar_node.getValue().getSeq().charAt(node.getEnd() - cigar_node.getStart())) {
					last_mut = node.getEnd();
					if (last_ref > 0) {
						cigar_end = (last_ref + last_mut) / 2 - cigar_node.getStart() + 1 + cigar_node.getValue().getStart();
						writeTwice(ref_writer, splitRecord(record, align_start, cigar_start, cigar_end));
						align_start = (last_ref + last_mut) / 2 + 1;
						cigar_start = cigar_end;
					}
					else if (last_cigar_ref) {
						writeTwice(ref_writer, splitRecord(record, align_start, cigar_start, cigar_end));
						align_start = cigar_node.getStart();
						cigar_start = cigar_node.getValue().getStart();
					}
					last_ref = 0;
					last_cigar_ref = false;
				}
				else {
					last_ref = node.getEnd();
					if (last_mut > 0) {
						cigar_end = (last_ref + last_mut) / 2 - cigar_node.getStart() + 1 + cigar_node.getValue().getStart();
						writeTwice(mut_writer, splitRecord(record, align_start, cigar_start, cigar_end));
						align_start = (last_ref + last_mut) / 2 + 1;
						cigar_start = cigar_end;
					}
					else if (last_cigar_mut) {
						writeTwice(mut_writer, splitRecord(record, align_start, cigar_start, cigar_end));
						align_start = cigar_node.getStart();
						cigar_start = cigar_node.getValue().getStart();
					}
					last_mut = 0;
					last_cigar_mut = false;
				}
			}
			last_cigar_mut = last_mut > 0;
			last_cigar_ref = last_ref > 0;
			last_mut = 0;
			last_ref = 0;
			cigar_end = cigar_node.getValue().getEnd();
		}
		if (!last_cigar_mut) {
			if (!last_cigar_ref) {
				writeWithoutSNP(ref_writer, mut_writer, record, -1.0);
			}
			else {
				writeTwice(ref_writer, cigar_start == 0 ? record : splitRecord(record, align_start, cigar_start, cigar_end));
			}
		}
		else if (!last_cigar_ref) {
			writeTwice(mut_writer, cigar_start == 0 ? record : splitRecord(record, align_start, cigar_start, cigar_end));
		}
	}
	
	private void writeWithoutSNP(SAMFileWriter ref_writer, SAMFileWriter mut_writer, SAMRecord record, double ref_chance) {
		if (0.0 > ref_chance || 1.0 < ref_chance) {
//			mut_writer.addAlignment(record);
//			ref_writer.addAlignment(record);
		}
		else {
			if (Math.random() < ref_chance) {
				ref_writer.addAlignment(record);
			}
			else {
				mut_writer.addAlignment(record);
			}
		}
	}
	
	private void writeTwice(SAMFileWriter writer, SAMRecord record) {
//		writer.addAlignment(record);
//		record.setReadName(record.getReadName() + ":X2");
//		writer.addAlignment(record);
	}
	
	private SAMRecord splitRecord(SAMRecord record, int align_start, int cigar_start, int cigar_end) {
		if (cigar_start == 0 && cigar_end == record.getReadLength()) {
			return record;
		}
		SAMRecord out = record.deepCopy();
		out.setAlignmentStart(align_start);
		out.setBaseQualityString(out.getBaseQualityString().substring(cigar_start, cigar_end));
		out.setReadString(out.getReadString().substring(cigar_start, cigar_end));
		out.setCigarString(buildSplitCigar(record.getCigarString(), cigar_start, cigar_end));
		return out;
	}
	
	public String buildSplitCigar(String cigar, int start, int end) {
		StringBuilder sb = new StringBuilder();
		int total_len = 0;
		int hard_clip = 0;
		int len = 0;
		int index_start = 0;
		for (int i = 0; i < cigar.length(); i++) {
			switch (cigar.charAt(i)) {
			case '0':
			case '1':
			case '2':
			case '3':
			case '4':
			case '5':
			case '6':
			case '7':
			case '8':
			case '9':
				len = len * 10 + cigar.charAt(i) - '0';
				break;

			case 'M':
			case '=':
			case 'X':
			case 'I':
			case 'S':
				if (start >= total_len && start < total_len + len) {
					if (start + hard_clip > 0) {
						sb.append(start + hard_clip);
						sb.append('H');
					}
					index_start = i + 1;
					sb.append(Math.min(total_len + len, end) - start);
					sb.append(cigar.charAt(i));
				}
				if (end > total_len && total_len + len >= end) {
					if (index_start < i) {
						sb.append(cigar.substring(index_start, i - (int) Math.log10(len) - 1));
						sb.append(end - total_len);
						sb.append(cigar.charAt(i));
					}
				}
				total_len += len;
				len = 0;
				hard_clip = 0;
				break;
			
			case 'H':
				hard_clip = len;
			case 'D':
			case 'N':
			case 'P':
				len = 0;
			default:
				break;
			}
		}
		if (total_len + hard_clip > end) {
			sb.append(total_len + hard_clip - end);
			sb.append('H');
		}
		return sb.toString();
	}
	
	public void compareMAF(HashMap<String, IntervalTree<Mutation>> mut_map, Genome genome) {
		for (Entry<String, Chromosome> entry : genome.getChrMap().entrySet()) {
			compareMAF(mut_map.get(entry.getKey()), entry.getValue());
		}
	}
	
	private void compareMAF(IntervalTree<Mutation> mut_tree, Chromosome chr) {
		if (mut_tree == null) {
			return;
		}
		Iterator<Node<Mutation>> nodes = mut_tree.iterator();
		while (nodes.hasNext()) {
			Node<Mutation> node = nodes.next();
			int gene_len = 0;
			Iterator<Node<Gene>> gene_nodes = chr.getGeneTree().overlappers(node.getStart(), node.getEnd());
			if (!gene_nodes.hasNext()) {
				node.getValue().setMaf(chr.getMaf());
			}
			while (gene_nodes.hasNext()) {
				Node<Gene> gene_node = gene_nodes.next();
				if (gene_node.getLength() <= gene_len) {
					continue;
				}
				node.getValue().setMaf(gene_node.getValue().getMAF());
			}
		}
//		chr.setMaf(calMAF(mut_tree, Integer.MIN_VALUE, Integer.MAX_VALUE));
	}
	
	public void calMAF(Map<String, IntervalTree<Mutation>> mut_map, Genome genome) {
		for (Entry<String, Chromosome> entry : genome.getChrMap().entrySet()) {
			calMAF(mut_map.get(entry.getKey()), entry.getValue());
		}
	}
	
	private void calMAF(IntervalTree<Mutation> mut_tree, Chromosome chr) {
		if (mut_tree == null) {
			return;
		}
		chr.setMaf(calMAF(mut_tree, Integer.MIN_VALUE, Integer.MAX_VALUE));
		Iterator<Node<Gene>> nodes = chr.getGeneTree().iterator();
		while (nodes.hasNext()) {
			Node<Gene> node = nodes.next();
			node.getValue().setMAF(calMAF(mut_tree, node.getStart(), node.getEnd()));
		}
	}
	
	public double calMAFWithRef(IntervalTree<Mutation> mut_tree, int start, int end) {
		double allele = 1.0;
		double total = 1.0;
		Iterator<Node<Mutation>> muts = mut_tree.overlappers(start, end);
		while (muts.hasNext()) {
			Node<Mutation> mut = muts.next();
			if (mut.getValue().isMutMajor()) {
				allele *= (1.0 - m.ref_fre);
				total *= m.ref_fre;
			}
			else {
				allele *= m.ref_fre;
				total *= (1.0 - m.ref_fre);
			}
		}
		return allele / allele + total;
	}
	
	public static Genome loadGenomeInfo(String exon_file, String seq_file) throws IOException {
		Genome genome = new Genome();
		genome.readAnnoteFile(exon_file);
		genome.readSeqFile(seq_file);
		return genome;
	}
	
//	public static boolean isNotSNP(String seq, int index, boolean positive, char ref) {
//		return positive ? seq.charAt(index) == ref : seq.charAt(seq.length() - 1 - index) == baseReverse(ref);
//	}
//	
//	private static boolean isSNP(String seq, int index, boolean positive, char mut) {
//		return positive ? seq.charAt(index) == mut : seq.charAt(seq.length() - 1 - index) == baseReverse(mut);
//	}
//	
//	private static char baseReverse(char c) {
//		switch (c) {
//		case 'A':
//			return 'T';
//		case 'T':
//			return 'A';
//		case 'C':
//			return 'G';
//		case 'G':
//			return 'C';
//		default:
//			return c;
//		}
//	}
	
//	public List<List<Gene>> readGeneFile(String file_name) throws IOException {
//		if (file_name == null) {
//			return null;
//		}
//		List<List<Gene>> out = new ArrayList<>(2);
//		out.add(new ArrayList<>());
//		out.add(new ArrayList<>());
//		try(BufferedReader reader = new BufferedReader(new FileReader(new File(file_name)))) {
//			String line = null;
//			while((line = reader.readLine()) != null) {
//				if (line.length() == 0 || line.charAt(0) == '#') {
//					continue;
//				}
//				String[] cols = line.split("\t");
//				if (cols.length >= 6 && cols[4].length() > 0) {
//					Gene gene = new Gene(null, null, null, null, '.', 0, 0);
//					switch (cols[4].charAt(0)) {
//					case 'f':
//					case 'F':
//						out.get(1).add(gene);
//						break;
//
//					case 't':
//					case 'T':
//						out.get(0).add(gene);
//						break;
//						
//					default:
//						continue;
//					}
//					gene.setIp_enrich(Double.parseDouble(cols[5]));
//					if (cols.length >= 8) {
//						gene.setInputReadCount(Integer.parseInt(cols[6]));
//						gene.setIpReadCount(Integer.parseInt(cols[7]));
//					}
//				}
//			}
//		}
//		return out;
//	}
	
	public static Map<String, Gene> readGeneFile(String file_name) throws IOException {
		if (file_name == null) {
			return null;
		}
		HashMap<String, Gene> out = new HashMap<>();
		try(BufferedReader reader = new BufferedReader(new FileReader(new File(file_name)))) {
			String line = null;
			while((line = reader.readLine()) != null) {
				if (line.length() == 0 || line.charAt(0) == '#') {
					continue;
				}
				String[] cols = line.split("\t");
				if (cols.length >= 4) {
					Gene gene = null;
					if (cols.length >= 6) {
						gene = new Gene(null, null, null, null, '.', 0, 0);
						gene.setIp_enrich(Double.parseDouble(cols[5]));
						if (cols.length >= 8) {
							gene.setInputReadCount(Integer.parseInt(cols[6]));
							gene.setIpReadCount(Integer.parseInt(cols[7]));
						}
					}
					out.put(cols[3].intern(), gene);
				}
			}
		}
		return out;
	}
	
	public HashMap<String, IntervalTree<Bed12>> readBedFile(String file_name, 
			HashMap<String, IntervalTree<Bed12>> out) throws IOException {
		out = out == null ? new HashMap<>() : out;
		try(BufferedReader reader = new BufferedReader(new FileReader(new File(file_name)))) {
			String line = null;
			while((line = reader.readLine()) != null) {
				if (line.length() == 0 || line.charAt(0) == '#') {
					continue;
				}
				String[] cols = line.split("\t");
				String chr = cols[0];
				int start = Integer.parseInt(cols[1]);
				int end = Integer.parseInt(cols[2]);
				Bed12 bed = null;
				if (!out.containsKey(chr)) {
					out.put(chr, new IntervalTree<>());
				}
				if (cols.length >= 12) {
					bed = new Bed12(start, end, chr, cols[3], Double.parseDouble(cols[4]), cols[5].charAt(0), 
						Integer.parseInt(cols[6]), Integer.parseInt(cols[7]), Integer.parseInt(cols[8]), Integer.parseInt(cols[9]),
						CommonMethod.lineToList(cols[10], ","), CommonMethod.lineToList(cols[11], ","));
				}
				else if (cols.length >= 6) {
					bed = new Bed12(start, end, chr, cols[3], Double.parseDouble(cols[4]), cols[5].charAt(0));
				}
				else {
					bed = new Bed12(chr, start, end);
				}
				if (start >= end) {
//					System.err.println(line);
					continue;
				}
				out.get(chr).put(start + 1, end, bed);
			}
		}
		return out;
	}
	
	public HashMap<String, IntervalTree<Bed12>> readBedFile(String file_name, boolean read_info, 
			HashMap<String, IntervalTree<Bed12>> out) throws IOException {
		out = out == null ? new HashMap<>() : out;
		try(BufferedReader reader = new BufferedReader(new FileReader(new File(file_name)))) {
			String line = null;
			while((line = reader.readLine()) != null) {
				if (line.length() == 0 || line.charAt(0) == '#') {
					continue;
				}
				String[] cols = line.split("\t");
				if (cols.length < 3) {
					continue;
				}
				String chr = cols[0];
				int start = Integer.parseInt(cols[1]);
				int end = start + 1;
				if (cols[2].matches("[0-9]*")) {
					end = Integer.parseInt(cols[2]);
				}
				Bed12 bed = null;
				if (!out.containsKey(chr)) {
					out.put(chr, new IntervalTree<>());
				}
				if (cols.length >= 12) {
					bed = new Bed12(start, end, chr, cols[3], Double.parseDouble(cols[4]), cols[5].charAt(0), 
						Integer.parseInt(cols[6]), Integer.parseInt(cols[7]), Integer.parseInt(cols[8]), Integer.parseInt(cols[9]),
						CommonMethod.lineToList(cols[10], ","), CommonMethod.lineToList(cols[11], ","));
				}
				else if (cols.length >= 6) {
					bed = new Bed12(start, end, chr, cols[3], Double.parseDouble(cols[4]), cols[5].charAt(0));
				}
				else {
					bed = new Bed12(chr, start, end);
				}
				if (read_info) {
					int index = -1;
					if (cols.length >= 12) {
						for (int i = 0; i < 12; i++) {
							if ((index = line.indexOf('\t', index + 1)) < 0) {
								break;
							}
						}
					}
					else if (cols.length >= 6) {
						for (int i = 0; i < 6; i++) {
							if ((index = line.indexOf('\t', index + 1)) < 0) {
								break;
							}
						}
					}
					else {
						for (int i = 0; i < 3; i++) {
							if ((index = line.indexOf('\t', index + 1)) < 0) {
								break;
							}
						}
					}
					bed.addInfo(line.substring(index + 1));
				}
				if (cols.length >= InParam.getParams().getExam() && InParam.getParams().getExam() > 0) {
					bed.setExamine(cols[InParam.getParams().getExam() - 1]);
				}
				if (start >= end) {
//					System.err.println(line);
					continue;
				}
				out.get(chr).put(start + 1, end, bed);
			}
		}
		return out;
	}
	
	public Map<String, IntervalTree<CircleClip>> readCircFile(String file_name, Map<String, 
			IntervalTree<CircleClip>> out, Genome genome) throws IOException {
		out = out == null ? new HashMap<>() : out;
		try(BufferedReader reader = new BufferedReader(new FileReader(new File(file_name)))) {
			String line = null;
			Map<String, Transcript> script_map = genome != null ? genome.getScriptMap() : new HashMap<>();
			while((line = reader.readLine()) != null) {
				if (line.length() == 0 || line.charAt(0) == '#') {
					continue;
				}
				String[] cols = line.split("\t");
				if (cols.length >= 3) {
					String chr = cols[0];
					int start = Integer.parseInt(cols[1]);
					int end = Integer.parseInt(cols[2]);
					Transcript script = cols.length > 3 ? script_map.get(cols[3]) : null;
					if (!out.containsKey(chr)) { 
						out.put(chr, new IntervalTree<>());
					}
					out.get(chr).put(start + 1, end, new CircleClip(start, end, script, cols.length > 4 ? isTrue(cols[4]) : false));
				}
			}
		}
		return out;
	}
	
	public static Map<String, IntervalTree<Mutation>> readMutation(String file_name, Map<String, IntervalTree<Mutation>> out) throws IOException {
		out = out == null ? new HashMap<>() : out;
		boolean vcf_flag = file_name.endsWith("vcf");
		try (BufferedReader reader = new BufferedReader(new FileReader(new File(file_name)))) {
			String line = null;
			while((line = reader.readLine()) != null) {
				if (line.length() == 0 || line.charAt(0) == '#') {
					continue;
				}
				readMutLine(out, line, vcf_flag);
			}
		}
		return out;
	}
	
	public static String readMutation(BufferedReader reader, Map<String, IntervalTree<Mutation>> out, boolean vcf_flag,
			String chr, String line) throws IOException {
		if (line == null || (line.length() > 0 && chr != null && !line.startsWith(chr))) {
			return line;
		}
		readMutLine(out, line, vcf_flag);
		while((line = reader.readLine()) != null) {
			if (line.length() == 0 || line.charAt(0) == '#') {
				continue;
			}
			if (!line.startsWith(chr)) {
				return line;
			}
			readMutLine(out, line, vcf_flag);
		}
		return line;
	}
	
	public static Map<String, IntervalTree<Mutation>> readMutWithLine(String file_name, int line_num) throws IOException {
		HashMap<String, IntervalTree<Mutation>> out = new HashMap<>();
		boolean vcf_flag = file_name.endsWith("vcf");
		int line_index = 0;
		try (BufferedReader reader = new BufferedReader(new FileReader(new File(file_name)))) {
			String line = null;
			while((line = reader.readLine()) != null) {
				if (line.length() == 0 || line.charAt(0) == '#') {
					continue;
				}
				if (++line_index >= line_num) {
					readMutLine(out, line, vcf_flag);
					break;
				}
			}
		}
		return out;
	}
	
	private static void readMutLine(Map<String, IntervalTree<Mutation>> out, String line, boolean vcf_flag) {
		if (line == null) {
			return;
		}
		String[] cols = line.split("\t");
		if (cols.length >= 5) {
			String chr = cols[0].intern();
			String id = null;
			int start = Integer.parseInt(cols[1]);
			int end = 0;
			if (vcf_flag) {
				id = cols[2];
				end = start--;
			}
			else {
				id = cols.length >= 6 ? cols[5] : id;
				end = Integer.parseInt(cols[2]);
			}
			char ref = cols[3].toUpperCase().charAt(0);
			String[] muts = cols[4].toUpperCase().split(",");
			char[] mut = new char[muts.length];
			for (int i = 0; i < muts.length; i++) {
				mut[i] = muts[i].charAt(0);
			}
			if (!out.containsKey(chr)) {
				out.put(chr, new IntervalTree<>());
			}
			Mutation mutation = new Mutation(ref, mut, id); 
			out.get(chr).put(start + 1, end, mutation);
		}
	}
	
	public static Map<String, IntervalTree<Scaffolds>> readScaffolds(String file, Map<String, IntervalTree<Scaffolds>> out) throws IOException {
		out = out == null ? new HashMap<>() : out;
		try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
			String line = null;
			while ((line = reader.readLine()) != null) {
				String[] cols = line.split("\t");
				String chr = cols[0].intern();
				out.computeIfAbsent(chr, k -> new IntervalTree<>());
				int start = Integer.parseInt(cols[1]);
				int end = Integer.parseInt(cols[2]);
				Scaffolds scaf = out.get(chr).put(start + 1, end, null);
				scaf = scaf == null ? new Scaffolds(start, end) : scaf;
				if (cols.length < 5) {
					scaf.putScaf(cols[3]);
				}
				else {
					scaf.putScaf(cols[3], Double.parseDouble(cols[4]));
				}
				out.get(chr).put(start + 1, end, scaf);
			}
		}
		return out;
	}
	
	public static void setStrictScript(Map<String, IntervalTree<Mutation>> mutTable, Genome genome, String chr) {
		if (mutTable == null || mutTable.get(chr) == null) {
			return;
		}
		if (genome == null || genome.getChr(chr) == null) {
			mutTable.remove(chr);
			return;
		}
		for (Node<Mutation> mut_node : mutTable.get(chr)) {
			Mutation mut = mut_node.getValue();
			for (Iterator<Node<Transcript>> script_nodes = genome.getChr(chr).getScrpitTree()
					.overlappers(mut_node.getStart(), mut_node.getEnd()); script_nodes.hasNext();) {
				Node<Transcript> script_node = script_nodes.next();
				for (Transcript script = script_node.getValue(); script != null; script = script.getAnotherScript()) {
					if (mut.getScript() == null || script.getExonLength() > mut.getScript().getExonLength()) {
						mut.setScript(script);
					}
				}
			}
			if (mut.getScript() == null) {
				mutTable.get(chr).remove(mut_node.getStart(), mut_node.getEnd());
			}
		}
	}
	
	public Map<String, IntervalTree<Peak>> readPeak(String file_name, Map<String, IntervalTree<Peak>> out) throws IOException {
		out = out == null ? new HashMap<>() : out;
		try (BufferedReader reader = new BufferedReader(new FileReader(new File(file_name)))) {
			String line = null;
			while ((line = reader.readLine()) != null) {
				if (line.length() == 0 || line.charAt(0) == '#') {
					continue;
				}
				String[] cols = line.split("\t");
				if (cols.length >= 3) {
					String chr = cols[0];
					int start = Integer.parseInt(cols[1]);
					int peak = Integer.parseInt(cols[2]);
					if (!out.containsKey(chr)) {
						out.put(chr, new IntervalTree<>());
					}
					out.get(chr).put(start + 1, peak, new Peak(peak, cols.length >= 4 ? cols[3] : null));
				}
			}
		}
		return out;
	}
	
	public HashMap<String, IntervalTree<Bed3Extend>> readBed3(String file_name, StringBuffer header) throws IOException {
		HashMap<String, IntervalTree<Bed3Extend>> out = new HashMap<>();
		try (BufferedReader reader = new BufferedReader(new FileReader(new File(file_name)))) {
			String line = null;
			while ((line = reader.readLine()) != null) {
				if (line.length() == 0 || line.charAt(0) == '#') {
					if (header != null) {
						header.append(line);
					}
					continue;
				}
				String[] cols = line.split("\t");
				if (cols.length >= 3) {
					String chr = cols[0];
					int start = Integer.parseInt(cols[1]);
					int end = Integer.parseInt(cols[2]);
					if (!out.containsKey(chr)) {
						out.put(chr, new IntervalTree<>());
					}
					int index = -1;
					for (int i = 0; i < 3; i++) {
						if ((index = line.indexOf('\t', index + 1)) < 0) {
							break;
						}
					}
					out.get(chr).put(start + 1, end, index >= 0 ? new Bed3Extend(chr, start, end, line.substring(index + 1)) : new Bed3Extend(chr, start, end));
				}
			}
		}
		return out;
	}
	
	public static double[] getCol(double[][] mat, int col) {
		double[] out = new double[mat.length];
		for (int i = 0; i < mat.length; ++i) {
			out[i] = mat[i][col];
		}
		return out;
	}
	
	public static double[][] readMatrix(String file) throws IOException {
		return readMatrix(file, "[ ,\t]+", 0);
	}
	
	public static double[][] readMatrix(String file, String sep, int skip) throws IOException {
		List<double[]> outL = new ArrayList<>();
		String line = null;
		try (BufferedReader reader = new BufferedReader(new FileReader(file))){
			while ((line = reader.readLine()) != null) {
				if (--skip >= 0) {
					continue;
				}
				String[] cols = line.split(sep);
				double[] nums = new double[cols.length];
				for (int i = 0; i < cols.length; ++i) {
					nums[i] = Double.parseDouble(cols[i]);
				}
				outL.add(nums);
			}
		}
		double[][] outM = new double[outL.size()][];
		for (int i = 0; i < outL.size(); ++i) {
			outM[i] = outL.get(i);
		}
		return outM;
	}
	
	public static List<SamReader> openBamFiles(List<String> bamFiles) {
		List<SamReader> readers = new ArrayList<>(bamFiles.size());
		for (int i = 0; i < bamFiles.size(); ++i) {
			SamReader reader = SamReaderFactory.makeDefault().open(new File(bamFiles.get(i)));
			if (!Method.isCoordinate(reader.getFileHeader())) {
				return readers;
			}
			readers.add(reader);
		}
		return readers;
	}
	
	public HashMap<String, IntervalTree<MappingStat>> readBamFile(String file_name) throws IOException{
		HashMap<String, IntervalTree<MappingStat>> out = new HashMap<>();
		try(SamReader reader = SamReaderFactory.makeDefault().open(new File(file_name))) {
			for (SAMRecord record : reader) {
				if (record.getReadUnmappedFlag()) {
					MappingStat.incReadNum(false);
					continue;
				}
				String chr = record.getReferenceName();
				if (!out.containsKey(chr)) {
					out.put(chr, new IntervalTree<>());
				}
				MappingMethod.putMappingStat(out.get(chr), record.getAlignmentStart(), record.getCigarString());
				countSNP(record);
			}
		}
		return out;
	}
	
	public <T extends IntRegion> HashMap<T, HashMap<Transcript, Integer>> seqBamScript(Genome genome, List<String> file_names, HashMap<String, IntervalTree<T>> table) throws IOException {
		if (file_names == null) {
			return null;
		}
		HashMap<T, HashMap<Transcript, Integer>> out = new HashMap<>();
		for (String file_name : file_names) {
			try(SamReader reader = SamReaderFactory.makeDefault().open(new File(file_name))) {
				String id = null;
				String seq = "";
				String chr = null;
				HashMap<T, Integer> points = new HashMap<>();
				for (SAMRecord record : reader) {
					if (!record.getReadName().equals(id)) {
						if (chr != null && !"".equals(chr)) {
							for (Entry<T, Integer> entry : points.entrySet()) {
								if (!out.containsKey(entry.getKey())) {
									out.put(entry.getKey(), new HashMap<>());
								}
								seqBamScript(out.get(entry.getKey()), genome.getChr(chr), seq, entry.getKey().getStart(), entry.getValue());
							}
						}
						seq = "";
						points = new HashMap<>();
						chr = null;
					}
					id = record.getReadName();
					if (record.getReadUnmappedFlag()) {
						continue;
					}
					if (chr != null && !chr.equals(record.getReferenceName())) {
						chr = "";
						continue;
					}
					chr = record.getReferenceName();
					getSeqIndex(points, record, table.get(chr));
					seq = record.getReadString().length() > seq.length() ? record.getReadString() : seq;
				}
			}
		}
		return out;
	}
	
	private void seqBamScript(HashMap<Transcript, Integer> script_map, Chromosome chr, String seq, int pos, int index) {
		for (Iterator<Node<Exon>> exon_nodes = chr.getExonTree().overlappers(pos + 1, pos + 1);
				exon_nodes.hasNext();) {
			Node<Exon> exon_node = exon_nodes.next();
			Transcript script = exon_node.getValue().getScript();
			int dis = script.disHammingExonSeq(seq, pos, index, chr.getSeq());
			if (index < 0 || index >= seq.length()) {
				System.out.println("debug");
			}
			dis -= seq.charAt(index) == chr.getSeq().charAt(pos) ? 0 : 1;
			if (dis == 0) {
				script_map.put(script, script_map.getOrDefault(script, 0) + 1);
			}
			else if (dis < 0) {
				System.err.println("Error Hamming Distance");
			}
		}
	}
	
	public <T extends IntRegion> void getSeqIndex(HashMap<T, Integer> map, SAMRecord record, IntervalTree<T> tree) {
		IntervalTree<Segment> seg_tree = BamMethod.cigarToSeq(null, record.getCigarString(), record.getAlignmentStart(), false);
		for (Node<Segment> seg_node : seg_tree) {
			for (Iterator<Node<T>> nodes = tree.overlappers(seg_node.getStart(), seg_node.getEnd()); 
					nodes.hasNext();) {
				Node<T> node = nodes.next();
				if (node.getStart() >= seg_node.getStart() && node.getEnd() <= seg_node.getEnd()) {
					int hard_clip = 'H' == CigarOperator.enumToBinary(record.getCigar().getCigarElement(0).getOperator()) ?
							record.getCigar().getCigarElement(0).getLength() : 0;
					map.put(node.getValue(), node.getStart() - seg_node.getStart() + seg_node.getValue().getStart() + hard_clip);
				}
			}
		}
	}
	
	public boolean isTrue(String s) {
		return s != null && s.length() > 0 && (s.charAt(0) == 'T' || s.charAt(0) == 't');
	}
	
	public static String toString(int[] array, String format, String sep) {
		StringBuilder sb = new StringBuilder();
		for (int i : array) {
			sb.append(String.format(format, i));
			sb.append(sep);
		}
		if (sb.length() > 0) {
			sb.setLength(sb.length() - sep.length());
		}
		return sb.toString();
	}
	
	public static String toString(double[] array, String format, String sep) {
		StringBuilder sb = new StringBuilder();
		for (double d : array) {
			sb.append(String.format(format, d));
			sb.append(sep);
		}
		if (sb.length() > 0) {
			sb.setLength(sb.length() - sep.length());
		}
		return sb.toString();
	}
	
//	public static double calORvalue(int in_a1, int in_a2, int ip_a1, int ip_a2) {
//		return (double) (ip_a1 * in_a2) / (double) (ip_a2 * in_a1);
//	}
	
//	public static double calPvalue(int in_a1, int in_a2, int ip_a1, int ip_a2) {
//		if (ip_a1 + ip_a2 + in_a1 + in_a2 == 0) {
//			return 1.0;
//		}
//		HypergeometricDistribution hd = new HypergeometricDistribution(ip_a1 + ip_a2 + in_a1 + in_a2, ip_a1 + in_a1, ip_a1 + ip_a2);
//		return hd.cumulativeProbability(ip_a1 - 1, ip_a1 + in_a1);
////		return (double) ip_a1 / (double) ip_a2 > (double) in_a1 / (double) in_a2 ? hd.cumulativeProbability(ip_a1 - 1, ip_a1 + in_a1)
////				: hd.cumulativeProbability(-1, ip_a1);
//	}
	
//	public static double calPvalue(double p, int ip_a1, int ip_a2) {
//		BinomialDistribution bd = new BinomialDistribution(ip_a1 + ip_a2, p);
//		boolean gain = (double) ip_a1 / (double) (ip_a1 + ip_a2) > p;
//		if (gain) {
//			return 1.0 - bd.cumulativeProbability(ip_a1 - 1);
//		}
//		else {
//			return bd.cumulativeProbability(ip_a1);
//		}
//	}
	
//	private void adjustP_Value(List<Double> p_value, String method){
//		if ("bon".equals(method)) {
//			double n = (double) p_value.size();
//			for (int i = 0; i < p_value.size() ; ++i) {
//				p_value.set(i, p_value.get(i) * n / (i + 1));
//			}
//		}
//		else if ("bh".equals(method)) {
//			List<double[]> sort = new ArrayList<>(p_value.size());
//			for (int i = 0; i < p_value.size(); ++i) {
//				sort.add(new double[] {p_value.get(i), i});
//			}
//			sort.sort((p1, p2) -> {
//				Double d = p1[0];
//				return d.compareTo(p2[0]);
//			});
//			double value = Double.MAX_VALUE;
//			for (int i = sort.size(); i > 0 ; --i) {
//				value = Math.min(value, sort.get(i - 1)[0] * (double) sort.size() / (double) i);
//				p_value.set((int) sort.get(i - 1)[1], value);
//			}
//		}
//		else {
//			System.out.println("Warning: unkown method of adjusting p-value, adjusting is disabled");
//		}
//	}
	
//	public static IntervalTree<Segment> cigarToSeq(SAMRecord record) {
//		return cigarToSeq(record, false);
//	}
//	
//	public static IntervalTree<Segment> cigarToSeq(SAMRecord record, boolean seqTreeFlag) {
//		if (record.getReadUnmappedFlag()) {
//			return new IntervalTree<>();
//		}
//		return cigarToSeq(record.getReadString(), record.getCigarString(), record.getAlignmentStart(), seqTreeFlag);
//	}
//	
//	public static IntervalTree<Segment> cigarToSeq(String seq, String cigar, int align, boolean seqTreeFlag) {
//		IntervalTree<Segment> out = new IntervalTree<>();
//		int start = 0;
//		int len = 0;
//		for (int i = 0; i < cigar.length(); i++) {
//			switch (cigar.charAt(i)) {
//			case '0':
//			case '1':
//			case '2':
//			case '3':
//			case '4':
//			case '5':
//			case '6':
//			case '7':
//			case '8':
//			case '9':
//				len = len * 10 + cigar.charAt(i) - '0';
//				break;
//
//			case 'M':
//			case '=':
//			case 'X':
//				if (seqTreeFlag) {
//					out.put(start + 1, start + len, new Segment(align - 1, align + len - 1, 
//							seq == null ? null : seq.substring(start, start + len)));
//				}
//				else {
//					out.put(align, align + len - 1, new Segment(start, start + len, 
//							seq == null ? null : seq.substring(start, start + len)));
//				}
//				align += len;
//			
//			case 'I':
//			case 'S':
//				start += len;
//				len = 0;
//				break;
//				
//			case 'D':
//				if (!seqTreeFlag) {
//					out.put(align, align + len - 1, new Segment(start, start, null));
//				}
//			case 'N':
//				align += len;
//				len = 0;
//				break;
//				
//			case 'H':
//			case 'P':
//				len = 0;
//			default:
//				break;
//			}
//		}
//		return out;
//	}
	
//	public static int indexOfKth(String source, String needle, int k) {
//		for (int index = source.indexOf(needle); index != -1; index = source.indexOf(needle, index + 1)) {
//			if (--k == 0) {
//				return index;
//			}
//		}
//		return -1;
//	}
	
//	public static <K, V> K randKeyWithValue(Map<K, V> map, V value) {
//		final List<K> keys = new ArrayList<>();
//		map.forEach((key, v) -> {
//			if (value != null && value.equals(v)) {
//				keys.add(key);
//			}
//			else if (v == null && value == null) {
//				keys.add(key);
//			}
//		});
//		return randElement(keys);
//	}
	
//	public static <T> T randElement(Collection<T> c) {
//		if (c == null || c.size() == 0) {
//			return null;
//		}
//		int index = (int) (Math.random() * (c.size())) + 1;
//		for (T t : c) {
//			if (--index == 0) {
//				return t;
//			}
//		}
//		return null;
//	}
	
//	public static <T> Collection<T> randNElement(Collection<T> c, int k) {
//		if (c == null || c.size() <= k) {
//			return c;
//		}
//		Set<Integer> indexs = randInSize(c.size(), k);
//		Set<T> out = new HashSet<>();
//		int index = 0;
//		for (T t : c) {
//			if (indexs.contains(index)) {
//				out.add(t);
//			}
//			++index;
//		}
//		return out;
//	}
	
//	public static Set<Integer> randInSize(int l, int k) {
//		Set<Integer> out = new HashSet<>();
//		if (l <= k) {
//			for (int i = 0; i < l; i++) {
//				out.add(i);
//			}
//			return out;
//		}
//		UniformIntegerDistribution uid = new UniformIntegerDistribution(0, l - 1);
//		boolean rm_flag = l - k > k;
//		k = rm_flag ? l - k : k;
//		if (rm_flag) {
//			for (int i = 0; i < l; i++) {
//				out.add(i);
//			}
//		}
//		while (k > 0) {
//			if (rm_flag) {
//				if (out.remove(uid.sample())) {
//					--k;
//				}
//			}
//			else {
//				if (out.add(uid.sample())) {
//					--k;
//				}
//			}
//		}
//		return out;
//	}
	
	public static HashMap<String, IntervalTree<Pair<Integer>>> countRegionInSam(HashMap<String, IntervalTree<Bed12>> regions,
			HashMap<String, IntervalTree<MappingStat>> mapping){
		HashMap<String, IntervalTree<Pair<Integer>>> out = new HashMap<>();
		for (Entry<String, IntervalTree<Bed12>> entry : regions.entrySet()) {
			if (mapping.containsKey(entry.getKey())) {
				out.put(entry.getKey(), new IntervalTree<>());
				Iterator<Node<Bed12>> nodes = entry.getValue().iterator();
				while (nodes.hasNext()) {
					IntervalTree.Node<Bed12> node = nodes.next();
					int strict = 0;
					int loose = 0;
					Iterator<Node<MappingStat>> map_nodes = mapping.get(entry.getKey()).overlappers(node.getStart(), node.getEnd());
					while (map_nodes.hasNext()) {
						IntervalTree.Node<MappingStat> map_node = map_nodes.next();
						loose += map_node.getValue().getReads();
						strict += map_node.getStart() <= node.getStart() && map_node.getEnd() >= node.getEnd() ? map_node.getValue().getReads() : 0;
					}
					out.get(entry.getKey()).put(node.getStart(), node.getEnd(), new Pair<>(strict, loose));
				}
			}
		}
		return out;
	}
	
	public static double calMAF(IntervalTree<Mutation> mut_tree, int start, int end) {
		double allele = 0.0;
		double total = 0.0;
		Iterator<Node<Mutation>> muts = mut_tree.overlappers(start, end);
		while (muts.hasNext()) {
			Node<Mutation> mut = muts.next();
			allele += mut.getValue().getMajorAlleleCount(false);
			total += mut.getValue().getTotalCount(false);
		}
		return allele / total;
	}
	
	public static void printNow(String prefix) {
		Date time = new Date();
		System.out.printf("%s at %tF %tT\n", prefix, time, time);
	}
	
	public static void logGlobal(String prefix) {
		Logger.getGlobal().info(prefix);
	}
	
	public static void writeMutation(String file_name, Genome genome, 
			HashMap<String, IntervalTree<Mutation>> out, String title) throws IOException{
		try(BufferedWriter writer = new BufferedWriter(new FileWriter(new File(file_name)))) {
			char sep = '\t';
			if (title != null) {
				writer.write(title);
				writer.newLine();
			}
			for (Entry<String, IntervalTree<Mutation>> entry : out.entrySet()) {
				Iterator<Node<Mutation>> nodes = entry.getValue().overlappers(Integer.MIN_VALUE, Integer.MAX_VALUE);
				String seq = genome != null && genome.getChr(entry.getKey()) != null && 
						genome.getChr(entry.getKey()).getSeq() != null ? genome.getChr(entry.getKey()).getSeq() : null;
				while (nodes.hasNext()) {
					IntervalTree.Node<Mutation> node = nodes.next();
					Mutation mut = node.getValue();
					if (mut.getCount(mut.getRef(), false) == 0 || mut.getMutCount(false) == 0
							|| node.getValue().getTotalCount(false) < 10) {
						continue;
					}
					writer.write(entry.getKey());
					writer.write(sep);
					writer.write(String.valueOf(node.getStart() - 1));
					writer.write(sep);
					writer.write(String.valueOf(node.getEnd()));
					writer.write(sep);
					writer.write(node.getValue().toSimpleString());
					Transcript script = node.getValue().getScript();
					if (script != null) {
						writer.write(sep);
						writer.write(String.valueOf(script.getGene().getStrand() == '+' ? script.inExons(node.getStart() - 1) + 1 :
							script.getExonLength() - script.inExons(node.getStart() - 1)));
						writer.write(sep);
						writer.write(script.getRegionFeature(node.getEnd(), true, true));
						writer.write(sep);
						writer.write(script.getUpStreamString(seq, node.getEnd(), 50));
						writer.write(sep);
						writer.write(script.getDownStreamString(seq, node.getEnd(), 50));
						writer.write(sep);
					}
					else {
						writer.write("\t0\tintergenic\t");
						if (seq != null && node.getEnd() >= 51 && node.getEnd() + 50 <= seq.length()) {
							writer.write(seq.substring(node.getEnd() - 51, node.getEnd() - 1));
							writer.write(sep);
							writer.write(seq.substring(node.getEnd(),node.getEnd() + 50));
							writer.write(sep);
						}
						else {
							writer.write(".\t.\t");
						}
					}
					if (seq != null) {
						writer.write(String.format("%.3f", BamMethod.getGCContent(seq, node.getEnd(), 1000)));
					}
					else {
						writer.write("0.000");
					}
					writer.newLine();
				}
				writer.flush();
			}
		}
	}
	
	public static <T> void writeFile(String file_name, Map<String, IntervalTree<T>> out, boolean with_key, boolean with_pos,
			String title) throws IOException{
		try(BufferedWriter writer = new BufferedWriter(new FileWriter(new File(file_name)))) {
			char sep = '\t';
			if (title != null) {
				writer.write(title);
				writer.newLine();
			}
			for (Entry<String, IntervalTree<T>> entry : out.entrySet()) {
				Iterator<Node<T>> nodes = entry.getValue().overlappers(Integer.MIN_VALUE, Integer.MAX_VALUE);
				while (nodes.hasNext()) {
					IntervalTree.Node<T> node = nodes.next();
					if (with_key) {
						writer.write(entry.getKey());
						writer.write(sep);
					}
					if (with_pos) {
						writer.write(String.valueOf(node.getStart()));
						writer.write(sep);
						writer.write(String.valueOf(node.getEnd()));
						writer.write(sep);
					}
					writer.write(node.getValue().toString());
					writer.newLine();
				}
				writer.flush();
			}
		}
	}
	
	public static <T> void writeFile(String file_name, Collection<T> set, String title) throws IOException {
		try(BufferedWriter writer = new BufferedWriter(new FileWriter(new File(file_name)))) {
			if (title != null) {
				writer.write(title);
				writer.newLine();
			}
			for (T t : set) {
				writer.write(t.toString());
				writer.newLine();
			}
			writer.flush();
		}
	}
	
	public static <T> void appendFile(String file_name, Collection<T> set, String title) throws IOException {
		try(BufferedWriter writer = new BufferedWriter(new FileWriter(new File(file_name), true))) {
			if (title != null) {
				writer.write(title);
				writer.newLine();
			}
			for (T t : set) {
				writer.write(t.toString());
				writer.newLine();
			}
			writer.flush();
		}
	}
}
