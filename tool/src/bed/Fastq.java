package bed;

import main.CommonMethod;
import main.InParam;
import sim.SimMethod;

public class Fastq {
	
	private String id = "@";
	private String seq = null;
	private String add_info = "+";
	private String quali = null;
	private String annote = null;
	
	public Fastq(String id, String seq, String add_info, String quali, String annote) {
		super();
		this.id = id;
		this.seq = seq;
		this.add_info = add_info;
		this.quali = quali;
		this.annote = annote;
	}

	public Fastq(String id, String seq, String quali, String annote) {
		super();
		this.id = id;
		this.seq = seq;
		this.quali = quali;
		this.annote = annote;
	}
	
	public Fastq(String id, String seq, String annote) {
		super();
		this.id = id;
		this.seq = seq;
		this.annote = annote;
		this.quali = SimMethod.buildQualiString(seq.length(), InParam.getParams().getQuality());
	}

	public String getID() {
		return id;
	}
	
	public void setAnnotation(String annote) {
		this.annote = annote;
	}
	
	public String getAnnotation() {
		return annote;
	}
	
	public int seqLen() {
		return this.seq.length();
	}
	
	public void merge(Fastq anotherFastq) {
		if (anotherFastq != null) {
			this.id = this.id.length() < anotherFastq.id.length() ? anotherFastq.id : this.id;
			this.add_info = this.add_info + anotherFastq.add_info.substring(1);
			this.seq = this.seq + anotherFastq.seq;
			this.quali = this.quali + anotherFastq.quali;
			if (anotherFastq.annote != null && CommonMethod.indexOfKth(anotherFastq.annote, "\t", 3) > 0) {
				this.annote = this.annote + anotherFastq.annote.substring(CommonMethod.indexOfKth(anotherFastq.annote, "\t", 3) + 1);
			}
		}
	}
	
	private boolean checkID() {
		return id.length() > 0 && id.charAt(0) == '@';
	}
	
	private boolean checkSeq() {
		for (int i = 0; i < seq.length(); i++) {
			switch(seq.charAt(i)) {
			case 'A':
			case 'T':
			case 'U':
			case 'G':
			case 'C':
			case 'N':
				break;
			default:
				return false;
			}
		}
		return true;
	}
	
	private boolean checkAdd_Info() {
		return add_info.length() > 0 && add_info.charAt(0) == '+';
	}
	
	private boolean checkQuality() {
		if (quali.length() != seq.length()) {
			return false;
		}
		for (int i = 0; i < quali.length(); i++) {
			if (quali.charAt(i) < 33 || quali.charAt(i) > 126) {
				return false;
			}
		}
		return true;
	}
	
	public boolean check() {
		return checkID() && checkSeq() && checkAdd_Info() && checkQuality();
	}

	@Override
	public String toString() {
		StringBuffer tmp = new StringBuffer();
		tmp.append(id);
		tmp.append('\n');
		tmp.append(seq);
		tmp.append('\n');
		tmp.append(add_info);
		tmp.append('\n');
		tmp.append(quali);
		return tmp.toString();
	}
	
	public String toStringWithCheck() {
		if (!check()) {
			System.err.println("Warning: Wrong read:\n" + toString());
			return null;
		}
		return toString();
	}
}
