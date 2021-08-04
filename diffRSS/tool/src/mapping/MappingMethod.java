package mapping;

import htsjdk.samtools.util.IntervalTree;

public class MappingMethod {

	public static void putMappingStat(IntervalTree<MappingStat> mapping_tree, int start, String cigar) {
		int[] cigar_value = cigarToRegion(cigar);
		MappingStat ms = mapping_tree.put(start, start + cigar_value[2] - 1, null);
		if (ms == null) {
			ms = new MappingStat();
		}
		ms.incReads(false, false);
		mapping_tree.put(start, start + cigar_value[2] - 1, ms);
		MappingStat.incReadNum(true);
	}
	
	public static void putMappingStat(IntervalTree<MappingStat> mapping_tree, int start, int end, boolean front, boolean behind) {
		MappingStat ms = mapping_tree.put(start, end, null);
		if (ms == null) {
			ms = new MappingStat();
		}
		ms.incReads(behind, front);
		mapping_tree.put(start, end, ms);
		MappingStat.incReadNum(true);
	}
	
	public static int[] cigarToRegion(String cigar) {
		boolean start_flag = false;
		int start = 1;
		int len = 0;
		int end = 0;
		int seg = 0;
		for (int i = 0; i < cigar.length(); ++i) {
			char c = cigar.charAt(i);
			switch (c) {
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
				seg = seg * 10 + c - '0';
				break;
				
			case 'M':
				len += seg;
				
			case 'I':
				start_flag = true;
				end += seg;
				seg = 0;
				break;
				
			case 'N':
			case 'D':
				len += seg;
				seg = 0;
				break;
				
			case 'S':
			case 'H':
				start += start_flag ? 0 : seg;
				end += start_flag ? 0 : seg;
				seg = 0;
				break;
				
			default:
				System.err.printf("Warning: %s has undefined char: %c\n", cigar, c);
				seg = 0;
				break;
			}
		}
		return new int[] {start, end, len};
	}
	
}
