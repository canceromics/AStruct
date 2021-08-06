package bed;

import java.util.ArrayList;
import java.util.List;

import genome.IntRegion;

public class Bed12 extends IntRegion{

	private String chr = null;
	private String name = null;
	private double score = 0.0;
	private char strand = '.';
	private int thick_start = 0;
	private int thick_end = 0;
	private int item_RGB = 0;
	private int block_count = 0;
	private List<Integer> block_sizes = null;
	private List<Integer> block_starts = null;
	private StringBuffer description = null;
	private String examine = null;
	
	public Bed12(String chr, int start, int end) {
		super(start, end);
		this.chr = chr;
	}
	
	public Bed12(int start, int end, String chr, String name, double score, char strand) {
		super(start, end);
		this.chr = chr;
		this.name = name;
		this.score = score;
		this.strand = strand;
	}

	public Bed12(int start, int end, String chr, String name, double score, char strand, int thick_start, int thick_end,
			int item_RGB, int block_count, List<Integer> block_sizes, List<Integer> block_starts) {
		super(start, end);
		this.chr = chr;
		this.name = name;
		this.score = score;
		this.strand = strand;
		this.thick_start = thick_start;
		this.thick_end = thick_end;
		this.item_RGB = item_RGB;
		this.block_count = block_count;
		this.block_sizes = block_sizes;
		this.block_starts = block_starts;
	}
	
	public char getStrand() {
		return strand;
	}

	public String getName() {
		return name;
	}
	
	public void setName(String name) {
		this.name = name;
	}
	
	public StringBuffer getDescription() {
		return description;
	}
	
	public void addInfo(String info) {
		if (description == null) {
			description = new StringBuffer();
		}
		description.append(info);
	}
	
	public void addInfo(StringBuffer info) {
		if (description == null) {
			description = new StringBuffer();
		}
		description.append(info);
	}
	
	public String getExamine() {
		return examine;
	}
	
	public void setExamine(String examine) {
		this.examine = examine;
	}
	
	public Bed12 deepClone() {
		Bed12 bed = block_count > 0 ? new Bed12(getStart(), getEnd(), chr, name, score, strand, 
				thick_start, thick_end, item_RGB, block_count, new ArrayList<>(block_sizes), new ArrayList<>(block_starts))
				: new Bed12(getStart(), getEnd(), chr, name, score, strand);
		bed.description = new StringBuffer(description);
		bed.examine = examine;
		return bed;
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
		if (name != null) {
			sb.append(sep);
			sb.append(name);
			sb.append(sep);
			sb.append(score);
			sb.append(sep);
			sb.append(strand);
			if (block_count > 0) {
				sb.append(sep);
				sb.append(thick_start);
				sb.append(sep);
				sb.append(thick_end);
				sb.append(sep);
				sb.append(item_RGB);
				sb.append(sep);
				sb.append(block_count);
				sb.append(sep);
				sb.append(block_sizes.get(0));
				for (int i = 1; i < block_count; i++) {
					sb.append(',');
					sb.append(block_sizes.get(i));
				}
				sb.append(sep);
				sb.append(block_starts.get(0));
				for (int i = 1; i < block_count; i++) {
					sb.append(',');
					sb.append(block_starts.get(i));
				}
			}
		}
		if (description != null) {
			sb.append(sep);
			sb.append(description);
		}
		return sb.toString();
	}
	
}
