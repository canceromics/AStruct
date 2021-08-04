package genome;

import java.util.HashMap;
import java.util.Map;

import bed.CircleClip;
import bed.Mutation;
import htsjdk.samtools.util.IntervalTree;

public class Transcript extends IntRegion{

	private Gene gene = null;
	private String script_id = null;
	private String script_type = null;
	private String exon_seq = null;
	private Map<String, Exon> exons = null;
	private int exon_length = 0;
	private int cover = 0;
	private Transcript same_position = null;
	private Transcript same_gene = null;
	private CircleClip circ = null;
	private Mutation mut = null;
	private double expression = 1.0;
	private boolean circ_flag = false;
	private boolean peak_flag = false;
	
	public Transcript(Gene gene, String script_id, String script_type, int start, int end) {
		super(start, end);
		this.gene = gene;
		this.script_id = script_id;
		this.script_type = script_type;
		this.exons = new HashMap<>();
	}
	
	public Transcript(Gene gene, String script_id, int start, int end, Map<String, Exon> exons) {
		super(start, end);
		this.gene = gene;
		this.script_id = script_id;
		this.exons = exons == null ? new HashMap<>() : exons;
	}

	public String getChr() {
		return gene == null ? null : gene.getChr();
	}
	
	public Gene getGene() {
		return gene;
	}
	
	public void setGene(Gene gene) {
		this.gene = gene;
	}
	
	public int getExonLength() {
		return exon_length == 0 ? setExonLength() : exon_length;
	}
	
	public Transcript getGeneScript() {
		return same_gene;
	}
	
	void setGeneScript(Transcript script) {
		same_gene = script;
	}
	
	public void addGeneScript(Transcript script) {
		if (script != null) {
			if (same_gene != null) {
				Transcript s = script;
				for (; s.same_gene != null; s = s.same_gene) {}
				s.same_gene = this.same_gene;
			}
			this.same_gene = script;
			if (gene != null) {
				gene.getChromosome().setNewScriptTree();
			}
		}
	}
	
	public Transcript getAnotherScript() {
		return same_position;
	}
	
	public void addAnotherScript(Transcript script) {
		if (script != null && same_position != null) {
			Transcript s = script;
			for (; s.same_position != null; s = s.same_position) {}
			s.same_position = this.same_position;
		}
		this.same_position = script;
	}
	
	private int setExonLength() {
		Exon exon = getExons();
		exon_length = 0;
		while (exon != null) {
			exon_length += exon.getLength();
			exon = exon.getDownStream();
		}
		return exon_length;
	}
	
	public String getScript_id() {
		return script_id;
	}

	public String getScript_type() {
		return script_type;
	}
	
	public Exon getExonHead(String key) {
		return exons.get(key);
	}
	
	public void addFirstExon(String key, Exon exon) {
		if (exon != null) {
			exon.addLastExon(exons.put(key, exon));
		}
	}
	
	public void addLastExon(String key, Exon exon) {
		if (!exons.containsKey(key)) {
			exons.put(key, exon);
		}
		else {
			exons.get(key).addLastExon(exon);
		}
	}
	
	public Exon getExons() {
		return exons.get("exon");
	}
	
	public void setExons(Exon exon) {
		exons.put("exon", exon);
	}
	
	public String getExonBases() {
		String bases = gene.getChromosome().getSeq();
		StringBuilder sb = new StringBuilder();
		for (Exon e = getExons(); e != null; e = e.getDownStream()) {
			sb.append(bases.substring(e.getStart(), e.getEnd()));
		}
		return sb.toString();
	}
	
	public int disHammingExonSeq(String seq, int pos, int index, String bases) {
		int script_index = inExons(pos);
		if (script_index < 0) {
			return -1;
		}
		int seq_start = index - (index < script_index ? index : script_index);
		int script_start = script_index - (script_index < index ? script_index : index);
		int seq_len = Math.min(index, script_index) + Math.min(seq.length() - index, getExonLength() - script_index);
		int dis = 0;
		for (int i = 0; i < seq_len; ++i) {
			dis += seq.charAt(i + seq_start) == getExonSeq().charAt(i + script_start) ? 0 : 1;
		}
		return dis;
	}
	
	public int getExonSeqSite(int index) {
		for (Exon e = getExons(); e != null; e = e.getDownStream()) {
			if (index < e.getLength()) {
				return index + e.getStart() + 1;
			}
			index -= e.getLength();
		}
		return -1;
	}
	
	public int inExons(int pos) {
		int out = 0;
		for (Exon e = getExons(); e != null; e = e.getDownStream()) {
			if (pos < e.getEnd() && pos >= e.getStart()) {
				out += pos - e.getStart();
				return out;
			}
			out += e.getLength();
		}
		return -1;
	}
	
	public Exon reachExon(int pos) {
		for (Exon e = getExons(); e != null; e = e.getDownStream()) {
			if (pos < e.getEnd() && pos >= e.getStart()) {
				return e;
			}
		}
		return null;
	}
	
	public void addExonInTree(IntervalTree<Exon> exon_tree) {
		if (exon_tree != null) {
			int start = Integer.MAX_VALUE;
			int end = Integer.MIN_VALUE;
			if (exons.containsKey("exon")) {
				exons.put("exon", exons.get("exon").sort());
			}
			for (Exon e = this.exons.get("exon"); e != null; e = e.getDownStream()) {
				e.addAnotherExon(exon_tree.put(e.getStart() + 1, e.getEnd(), e));
				start = Math.min(start, e.getStart());
				end = Math.max(end, e.getEnd());
			}
			if (getStart() == getEnd()) {
				resetStartAndEnd(start, end);
				if (gene != null) {
					gene.getChromosome().setNewScriptTree();
				}
			}
		}
	}
	
	public boolean sort() {
		if (getExons() == null) {
			return false;
		}
		exons.put("exon", getExons().sort());
		resetStartAndEnd();
		resetUTR();
		setExonLength();
		return true;
	}
	
	public void resetStartAndEnd() {
		if (getStart() == getEnd()) {
			int start = Integer.MAX_VALUE;
			int end = Integer.MIN_VALUE;
			for (Exon e = getExons(); e != null; e = e.getDownStream()) {
				start = Math.min(start, e.getStart());
				end = Math.max(end, e.getEnd());
			}
			resetStartAndEnd(start, end);
			if (gene != null) {
				gene.getChromosome().setNewScriptTree();
			}
		}
	}

	public void resetUTR() {
		if (exons.containsKey("CDS")) {
			exons.put("CDS", exons.get("CDS").sort());
		}
		else {
			return;
		}
		exons.remove("UTR");
		boolean firstCDS = true;
		for (Exon cds = exons.get("CDS"), exon = exons.get("exon"); exon != null;exon = exon.getDownStream()) {
			if (firstCDS) {
				if (exon.getEnd() <= cds.getStart()) {
					addLastExon("UTR", new Exon(this, exon.getStart(), exon.getEnd()));
				}
				else {
					if (exon.getStart() < cds.getStart()) {
						addLastExon("UTR", new Exon(this, exon.getStart(), cds.getStart()));
					}
					while (cds.getDownStream() != null) {
						cds = cds.getDownStream();
					}
					firstCDS = false;
				}
			}
			if (cds.getEnd() < exon.getEnd()) {
				addLastExon("UTR", new Exon(this, Math.max(cds.getEnd(), exon.getStart()), exon.getEnd()));
			}
		}
	}
	
	public boolean isCirc() {
		return circ != null;
	}

	public boolean isMut() {
		return mut != null;
	}
	
	public boolean isCircTrans() {
		return circ_flag;
	}
	
	public void setMut(Mutation mut) {
		this.mut = mut;
	}
	
	public void setOnceMut(Mutation mut) {
		if (this.mut == null) {
			this.mut = mut;
		}
	}
	
	public Mutation getMut() {
		return mut;
	}
	
	public boolean isPeak() {
		return peak_flag;
	}
	
	public CircleClip getCirc() {
		return circ;
	}
	
	public String getExonSeq() {
		return exon_seq == null ? exon_seq = getExonBases() : exon_seq;
	}
	
	public void setCirc(CircleClip circ) {
		this.circ = circ;
	}
	
	public void setCircTrans(boolean circ_flag) {
		this.circ_flag = circ_flag;
	}
	
	public void setPeak(boolean peak_flag) {
		this.peak_flag = peak_flag;
	}
	
	public void setExpression(double expression) {
		this.expression = expression;
	}
	
	public double getExpression() {
		return expression;
	}
	
	public void setCover(int cover) {
		this.cover = cover;
	}
	
	public int getCover() {
		return cover;
	}
	
	public void incCover() {
		++cover;
	}
	
	public <T extends IntRegion> String getRegionFeature(T bed) {
		return getRegionFeature(bed, true, true);
	}
	
	public String getRegionFeature(int pos) {
		return getRegionFeature(pos, true, true);
	}
	
	public <T extends IntRegion> String getRegionFeature(T bed, boolean exon_flag, boolean intron_flag) {
		StringBuffer sb = new StringBuffer();
		if (this.getExonHead("UTR") != null) {
			boolean first_utr = false;
			boolean last_utr = false;
			for (Exon exon = getExonHead("UTR"); exon != null; exon = exon.getDownStream()) {
				first_utr |= bed.getStart() >= exon.getStart() && bed.getStart() < exon.getEnd();
				last_utr |= bed.getEnd() > exon.getStart() && bed.getEnd() <= exon.getEnd();
			}
			if (getGene().getStrand() == '-') {
				sb.append(last_utr ? "5'UTR," : "");
				sb.append(first_utr ? "3'UTR," : "");
			}
			else {
				sb.append(first_utr ? "5'UTR," : "");
				sb.append(last_utr ? "3'UTR," : "");
			}
		}
		if (this.getExonHead("CDS") != null) {
			for (Exon exon = getExonHead("CDS"); exon != null; exon = exon.getDownStream()) {
				if (bed.getEnd() > exon.getStart() && exon.getEnd() > bed.getStart()) {
					sb.append("CDS,");
					break;
				}
			}
		}
		if (sb.length() == 0 && exon_flag) {
			for (Exon exon = this.getExons(); exon != null; exon = exon.getDownStream()) {
				if (bed.getEnd() > exon.getStart() && exon.getEnd() > bed.getStart()) {
					sb.append("exon,");
					break;
				}
			}
		}
		if (sb.length() == 0) {
			if (intron_flag) {
				sb.append("intron");
			}
			else {
				sb.append('.');
			}
		}
		else {
			sb.setLength(sb.length() - 1);
		}
		return sb.toString();
	}
	
	public String getRegionFeature(int pos, boolean exon_flag, boolean intron_flag) {
		StringBuffer sb = new StringBuffer();
		if (this.getExonHead("UTR") != null) {
			Exon utr = null;
			for (Exon exon = this.getExonHead("UTR"); exon != null; exon = exon.getDownStream()) {
				utr = pos > exon.getStart() && pos <= exon.getEnd() ? exon : utr;
			}
			sb.append(getUTR(utr));
		}
		if (this.getExonHead("CDS") != null) {
			for (Exon exon = this.getExonHead("CDS"); exon != null; exon = exon.getDownStream()) {
				if (pos > exon.getStart() && exon.getEnd() >= pos) {
					sb.append("CDS,");
					break;
				}
			}
		}
		if (sb.length() == 0 && exon_flag) {
			for (Exon exon = this.getExons(); exon != null; exon = exon.getDownStream()) {
				if (pos > exon.getStart() && exon.getEnd() >= pos) {
					sb.append("exon,");
					break;
				}
			}
		}
		if (sb.length() == 0) {
			if (intron_flag) {
				sb.append("intron");
			}
			else {
				sb.append('.');
			}
		}
		else {
			sb.setLength(sb.length() - 1);
		}
		return sb.toString();
	}
	
	private String getUTR(Exon exon) {
		if (exon == null) {
			return "";
		}
		if (exon.getStart() - getStart() < getEnd() - exon.getEnd()) {
			return gene.getStrand() == '-' ? "3'UTR," : "5'UTR,";
		}
		else {
			return gene.getStrand() == '-' ? "5'UTR," : "3'UTR,";
		}
	}
	
	public String getUpStreamString(String chr_seq, int pos, int len) {
		if (chr_seq == null) {
			return ".";
		}
		int script_index = 0;
		for (Exon exon = getExons(); exon != null; exon = exon.getDownStream()) {
			if (pos > exon.getStart() && pos <= exon.getEnd()) {
				script_index += pos - exon.getStart() - 1;
				break;
			}
			script_index += exon.getEnd() - exon.getStart();
		}
		if (script_index == getExonLength()) {
			return pos < 51 ? "." : chr_seq.substring(pos - 51, pos - 1);
		}
		else {
			int start_index = Math.max(0, script_index - len);
			return getExonBases().substring(start_index, script_index);
		}
	}
	
	public String getDownStreamString(String chr_seq, int pos, int len) {
		if (chr_seq == null) {
			return ".";
		}
		int script_index = 0;
		for (Exon exon = getExons(); exon != null; exon = exon.getDownStream()) {
			if (pos > exon.getStart() && pos <= exon.getEnd()) {
				script_index += pos - exon.getStart();
				break;
			}
			script_index += exon.getEnd() - exon.getStart();
		}
		if (script_index == getExonLength()) {
			return pos + 50 > chr_seq.length() ? "." : chr_seq.substring(pos, pos + 50);
		}
		else {
			int end_index = Math.min(getExonLength(), script_index + len);
			return getExonBases().substring(script_index, end_index);
		}
	}
	
	public String getExonGCContent(String feature, int pos) {
		StringBuilder sb = new StringBuilder();
		int[] count = new int[4];
		int stat = 0;
		switch (feature) {
		case "exon":
			stat = 0;
			break;
			
		case "CDS":
			stat = 4;
			break;

		case "UTR":
			stat = 3;
			break;
			
		case "5'UTR":
			stat = 1;
			break;
			
		case "3'UTR":
			stat = 2;
			break;
			
		default:
			stat = -1;
			break;
		}
		for (Exon exon = exons.get("CDS"); exon != null; exon = exon.getDownStream()) {
			count[1] += exon.getLength();
			if (stat == 4) {
				count[3] += exon.getGCBases();
			}
		}
		for (Exon exon = exons.get("UTR"); exon != null; exon = exon.getDownStream()) {
			if (exon.getStart() - getStart() < getEnd() - exon.getEnd() ^ gene.getStrand() == '-') {
				count[0] += exon.getLength();
				if ((stat & 1) > 0) {
					count[3] += exon.getGCBases();
				}
			}
			else {
				count[2] += exon.getLength();
				if ((stat & 2) > 0) {
					count[3] += exon.getGCBases();
				}
			}
		}
		if (stat == 0) {
			for (Exon exon = exons.get("exon"); exon != null; exon = exon.getDownStream()) {
				count[3] += exon.getGCBases();
			}
		}
		if (stat < 0) {
			return getIntronGC(pos);
		}
		sb.append(count[0]);
		sb.append('|');
		sb.append(count[1]);
		sb.append('|');
		sb.append(count[2]);
		sb.append('|');
		int baseSum = stat == 0 ? getExonLength() : stat == 4 ? count[1] : 0;
		baseSum += (stat & 1) > 0 ? count[0] : 0;
		baseSum += (stat & 2) > 0 ? count[2] : 0;
		sb.append(String.format("%.4f", (double) count[3] / (double) baseSum));
		return sb.toString();
	}
	
	private String getIntronGC(int pos) {
		if (pos <= getStart() || pos > getEnd()) {
			return "0|0|0|-1.000";
		}
		else {
			int inStart = 0;
			int inEnd = 0;
			for (Exon exon = exons.get("exon"); exon != null; exon = exon.getDownStream()) {
				if (exon.getEnd() < pos) {
					inStart = exon.getEnd();
				}
				if (exon.getStart() >= pos) {
					inEnd = exon.getStart();
					break;
				}
			}
			int sum = 0;
			for (int i = inStart; i < inEnd; ++i) {
				char c = Character.toUpperCase(getGene().getChromosome().getSeq().charAt(i));
				sum += c == 'G' || c == 'C' ? 1 : 0;
			}
			return "0|0|0|" + String.format("%.4f", (double) sum / (double) (inEnd - inStart));
		}
	}
	
	@Override
	public String toString() {
		return script_id;
	}

}
