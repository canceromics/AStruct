package bed;

public class Junction {

	private String chr;
	private int start;
	private int end;
	
	public Junction(String chr, int start, int end) {
		this.chr = chr;
		this.start = start;
		this.end = end;
	}

	public String getChr() {
		return chr;
	}
	
	public int getStart() {
		return start;
	}
	
	public int getEnd() {
		return end;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(chr);
		sb.append('\t');
		sb.append(start);
		sb.append('\t');
		sb.append(end);
		return sb.toString();
	}
}
