package bed;

import genome.Exon;
import genome.IntRegion;
import genome.Transcript;
import main.InParam;

public class CircleClip extends IntRegion {

	private Peak peak = null;
	private int index = Integer.MAX_VALUE;
	private boolean exon_bounder = false;
	private boolean gt_ag = false;
	private boolean visited = false;
//	private int stat = 0;
	private int linear_ip = 0;
	private int linear_input = 0;
	private int circ_ip = 0;
	private int circ_input = 0;
	private Transcript script = null;
	
	public CircleClip(int start, int end, Transcript script, boolean gt_ag) {
		super(start, end);
		this.script = script;
		this.gt_ag = gt_ag;
	}

	public void setScript(Transcript t) {
		this.script = t;
	}
	
	public Transcript getScript() {
		return script;
	}
	
	public Peak getPeak() {
		return peak;
	}
	
	public int getPeakDis() {
		return peak == null ? -1 : Math.min(index, (script == null ? getEnd() : script.getExonLength()) - index);
	}
	
	public boolean isNearPeak() {
		return peak != null && (index < InParam.getParams().getReadLength() 
				|| index + InParam.getParams().getReadLength() > (script == null ? getEnd() : script.getExonLength()));
	}

	public boolean isFarPeak() {
		return peak != null && !isNearPeak();
	}
	
	public void setPeak(Peak peak) {
		this.peak = peak;
	}
	
	public void setPeak(Peak peak, int index) {
		this.peak = peak;
		this.index = index;
	}
	
	public void setPeakIndex() {
		int len = 0;
		for (Exon exon = script.getExons(); exon != null; exon = exon.getDownStream()) {
			if (peak.getStart() >= exon.getStart() && peak.getEnd() <= exon.getEnd()) {
				index = peak.getStart() - exon.getStart() + len;
				break;
			}
			len += exon.getLength();
		}
	}
	
	public void setGTAG(boolean b) {
		this.gt_ag = b;
	}
	
//	public void setStat(int stat) {
//		this.stat = stat;
//	}
	
	public boolean isGTAG() {
		return gt_ag;
	}
	
	public boolean isExon() {
		return exon_bounder;
	}

	public void setExon(boolean exon_bounder) {
		this.exon_bounder = exon_bounder;
	}

//	public int getLinear() {
//		return linear;
//	}
	
	public int incCircNum(boolean circ_flag, boolean ip_flag) {
		return circ_flag ? (ip_flag ? ++circ_ip : ++circ_input) : (ip_flag ? ++linear_ip : ++linear_input);
	}
	
//	public int getCirc() {
//		return circ;
//	}
	
//	public int getStat() {
//		return stat;
//	}
	
	public void visited() {
		visited = true;
	}
	
	public boolean isVisited() {
		return visited;
	}
	
	public static String getHeader() {
		return "#Chr\tStart\tEnd\tTranscript\tDis_m6A\tBSJ_T\tBSJ_C\tlinear_T\tlinear_C";
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(super.toString());
		sb.append('\t');
		sb.append(script);
//		sb.append('\t');
//		sb.append(exon_bounder);
//		sb.append('\t');
//		sb.append(gt_ag);
//		sb.append('\t');
//		sb.append(peak == null ? "None" : isNearPeak() ? "Near" : "Far");
		sb.append('\t');
		sb.append(getPeakDis());
//		sb.append('\t');
//		sb.append(stat > 1 ? (stat == 2 ? "Script" : "Gene") : (stat == 0 ? "None" : "Exon"));
		sb.append('\t');
		sb.append(circ_ip);
		sb.append('\t');
		sb.append(circ_input);
		sb.append('\t');
		sb.append(linear_ip);
		sb.append('\t');
		sb.append(linear_input);
		
		return sb.toString();
	}

}
