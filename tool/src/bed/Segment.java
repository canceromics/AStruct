package bed;

import genome.IntRegion;

public class Segment extends IntRegion{

	private String seq = null;
	private boolean positive = true;
	
	public Segment(int start, int end, String seq) {
		super(start, end);
		this.seq = seq;
	}

	public String getSeq() {
		return seq;
	}

	public boolean isPositive() {
		return positive;
	}

	public void setPositive(boolean positive) {
		this.positive = positive;
	}
}
