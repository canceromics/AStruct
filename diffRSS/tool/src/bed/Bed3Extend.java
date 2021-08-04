package bed;

import genome.IntRegion;

public class Bed3Extend extends IntRegion{
	
	private String chr = null;
	private StringBuffer description = null;
	
	public Bed3Extend(String chr, int start, int end) {
		super(start, end);
		this.chr = chr;
		this.description = new StringBuffer();
	}
	
	public Bed3Extend(String chr, int start, int end, String description) {
		super(start, end);
		this.chr = chr;
		this.description = new StringBuffer();
		this.description.append(description);
	}
	
	public StringBuffer getDescription() {
		return description;
	}
	
	public void appDes(String des) {
		this.description.append(des);
	}
	
	@Override
	public String toString() {
		char sep = '\t';
		StringBuilder sb = new StringBuilder();
		sb.append(chr);
		sb.append(sep);
		sb.append(getStart());
		sb.append(sep);
		sb.append(getEnd());
		sb.append(sep);
		sb.append(description);
		return sb.toString();
	}
}
