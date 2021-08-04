package main;

import bed.Segment;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.IntervalTree;

public class BamMethod {

	private static final char[] BASES = {'A', 'T', 'C', 'G'};
	
	public static IntervalTree<Segment> cigarToSeq(SAMRecord record) {
		return cigarToSeq(record, false);
	}
	
	public static IntervalTree<Segment> cigarToSeq(SAMRecord record, boolean seqTreeFlag) {
		if (record.getReadUnmappedFlag()) {
			return new IntervalTree<>();
		}
		return cigarToSeq(record.getReadString(), record.getCigarString(), record.getAlignmentStart(), seqTreeFlag);
	}
	
	public static IntervalTree<Segment> cigarToSeq(String seq, String cigar, int align, boolean seqTreeFlag) {
		IntervalTree<Segment> out = new IntervalTree<>();
		int start = 0;
		int len = 0;
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
				if (seqTreeFlag) {
					out.put(start + 1, start + len, new Segment(align - 1, align + len - 1, 
							seq == null ? null : seq.substring(start, start + len)));
				}
				else {
					out.put(align, align + len - 1, new Segment(start, start + len, 
							seq == null ? null : seq.substring(start, start + len)));
				}
				align += len;
			
			case 'I':
			case 'S':
				start += len;
				len = 0;
				break;
				
			case 'D':
				if (!seqTreeFlag) {
					out.put(align, align + len - 1, new Segment(start, start, null));
				}
			case 'N':
				align += len;
				len = 0;
				break;
				
			case 'H':
			case 'P':
				len = 0;
			default:
				break;
			}
		}
		return out;
	}
	
	public static double getGCContent(String chr_seq, int pos, int len) {
		int sum = 0;
		int start = Math.max(0, pos - len - 1);
		int end = Math.min(chr_seq.length(), pos + len);
		for (int i = start; i < end; ++i) {
			char c = chr_seq.charAt(i);
			sum += c == 'C' || c == 'G' ? 1 : 0;
		}
		return (double) sum / (double) (end - start);
	}
	
	public static char getMutBase(char ref) {
		char alt = BASES[CommonMethod.randInt(BASES.length)];
		while (alt == ref) {
			alt = BASES[CommonMethod.randInt(BASES.length)];
		}
		return alt;
	}
	
}
