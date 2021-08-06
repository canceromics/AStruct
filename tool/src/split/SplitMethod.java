package split;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;
import java.util.logging.Logger;

import bed.Mutation;
import bed.Segment;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import main.BamMethod;
import main.InParam;
import main.Method;

public class SplitMethod {

	public static void run() {
		try {
			Map<String, IntervalTree<Mutation>> mut_table = Method.readMutWithLine(InParam.getParams().getMutFile(),
					InParam.getParams().getMinBackRead());
			if (mut_table.size() < 1) {
				Logger.getLogger("error").severe("None SNP loaded");
				return;
			}
			SplitMethod spm = new SplitMethod();
			spm.splitBam(InParam.getParams().getTreatFile(), mut_table, InParam.getParams().getOutPrefix(), false);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void splitBam(String file_name, Map<String, IntervalTree<Mutation>> mut_table, String out_prefix, boolean bam) {
		if (file_name == null || out_prefix == null) {
			return;
		}
		try (SamReader reader = SamReaderFactory.makeDefault().open(new File(file_name));
			SAMFileWriter ref_writer = bam ? new SAMFileWriterFactory().makeBAMWriter(Method.setUnsort(reader.getFileHeader()),
				true, new File(out_prefix + "_ref.bam")) : new SAMFileWriterFactory().makeSAMWriter(
					Method.setUnsort(reader.getFileHeader()), true, new File(out_prefix + "_ref.sam"));
			SAMFileWriter alt_writer = bam ? new SAMFileWriterFactory().makeBAMWriter(Method.setUnsort(reader.getFileHeader()),
				true, new File(out_prefix + "_alt.bam")) : new SAMFileWriterFactory().makeSAMWriter(
					Method.setUnsort(reader.getFileHeader()), true, new File(out_prefix + "_alt.sam"))) {
			for (SAMRecord record : reader) {
				switch (aroundSNPs(record, mut_table)) {
				case 1:
					ref_writer.addAlignment(record);
					break;

				case 2:
					alt_writer.addAlignment(record);
					break;
					
				default:
					break;
				}
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void combineSam(String f1, String f2, String fo) {
		try (SamReader reader1 = SamReaderFactory.makeDefault().open(new File(f1));
			SamReader reader2 = SamReaderFactory.makeDefault().open(new File(f2));
			SAMFileWriter writer = new SAMFileWriterFactory().makeSAMWriter(Method.setUnsort(reader1.getFileHeader()),
				true, new File(fo))) {
			Iterator<SAMRecord> it1 = reader1.iterator();
			Iterator<SAMRecord> it2 = reader2.iterator();
			while (it1.hasNext() || it2.hasNext()) {
				if (it1.hasNext()) {
					writer.addAlignment(it1.next());
				}
				if (it2.hasNext()) {
					writer.addAlignment(it2.next());
				}
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private int aroundSNPs(SAMRecord record, Map<String, IntervalTree<Mutation>> mut_table) {
		if (mut_table == null || !mut_table.containsKey(record.getReferenceName())) {
			return -1;
		}
		IntervalTree<Mutation> mut_tree = mut_table.get(record.getReferenceName());
		IntervalTree<Segment> seg_tree = BamMethod.cigarToSeq(record);
		for (Node<Segment> seg_node : seg_tree) {
			Segment seg = seg_node.getValue();
			Node<Mutation> mut_node = mut_tree.minOverlapper(seg_node.getStart(), seg_node.getEnd());
			if (mut_node == null) {
				continue;
			}
			if (record.getCigarLength() > 1) {
				Logger.getLogger("debug").info("debug");
			}
			if (seg.getStart() == seg.getEnd()) {
				return 2;
			}
			char record_char = seg.getSeq().charAt(mut_node.getEnd() - seg_node.getStart());
			if (record_char == mut_node.getValue().getRef()) {
				return 1;
			}
			else if (mut_node.getValue().isMut(record_char)) {
				return 2;
			}
		}
		return -1;
	}
	
	public static void main(String[] args) {
		SplitMethod spm = new SplitMethod();
		spm.combineSam("D:/work/workspace/data/Circ_sim_alt.sam", "D:/work/workspace/data/Circ_sim_ref.sam",
				"D:/work/workspace/data/combine.sam");
	}
}
