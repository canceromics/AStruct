package mapping;

public class MappingStat {
	
	private static int input_reads = 0;
	private static int mapping_input_reads = 0;
	
	private int reads = 0;
	private int both = 0;
	private int no_behind = 0;
	private int no_front = 0;
	
	public void incReads(boolean behind_flag, boolean front_flag) {
		++reads;
		no_behind += behind_flag ? 0 : 1;
		no_front += front_flag ? 0 : 1;
		both += behind_flag || front_flag ? 0 : 1;
	}
	
	public void resetReads() {
		reads = 0;
		both = 0;
		no_behind = 0;
		no_front = 0;
	}
	
	public int getReads() {
		return reads;
	}
	
	public int getNoBehindReads() {
		return no_behind;
	}
	
	public int getNoFrontReads() {
		return no_front;
	}
	
	public int getBothReads() {
		return both;
	}
	
	public static int getInputReads() {
		return input_reads;
	}
	
	public static int getMappingInputReads() {
		return mapping_input_reads;
	}
	
	public static void incReadNum(boolean mapping_flag) {
		if (mapping_flag) {
			++mapping_input_reads;
		}
		++input_reads;
	}
	
	public static void resetReadNum() {
		mapping_input_reads = 0;
		input_reads = 0;
	}
	
}
