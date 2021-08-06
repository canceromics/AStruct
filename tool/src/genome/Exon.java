package genome;

public class Exon extends IntRegion{

	private Transcript script = null;
	private Exon same_position = null;
	private Exon next_exon = null;
	private int score = 0;
	
	public Exon(Transcript script, int start, int end) {
		super(start, end);
		this.script = script;
	}
	
	public void addLastExon(Exon exon) {
		Exon e = this;
		while (e.next_exon != null) {
			e = e.next_exon;
		}
		e.next_exon = exon;
		if (script != null && script.getGene() != null) {
			script.getGene().getChromosome().setNewExonTree();
		}
	}
	
	public void setDownStream(Exon exon) {
		next_exon = exon;
	}
	
	public Exon sort() {
		Exon head = new Exon(null, 0, 0);	
		for (Exon p = this; p != null;) {
			Exon q = p.next_exon;
//			p.next_exon = null;
			for (Exon tail = head; ; tail = tail.next_exon) {
				if (tail.next_exon == null || tail.next_exon.getStart() > p.getStart()) {
					Exon e = tail.next_exon;
					tail.next_exon = p;
					p.next_exon = e;
					break;
				}
			}
			p = q;
		}
		return head.next_exon;
	}
	
	int getGCBases() {
		int sum = 0;
		for (int i = getStart(); i < getEnd(); ++i) {
			char c = Character.toUpperCase(getGene().getChromosome().getSeq().charAt(i));
			sum += c == 'G' || c == 'C' ? 1 : 0;
		}
		return sum;
	}
	
	public Exon getDownStream() {
		return next_exon;
	}
	
	public void addAnotherExon(Exon exon) {
		this.same_position = exon;
	}
	
	public Exon getAnotherExon() {
		return same_position;
	}
	
	public String getChr() {
		return script == null ? null : script.getChr();
	}
	
	public Gene getGene() {
		return script == null ? null : script.getGene();
	}
	
	public Transcript getScript() {
		return script;
	}
	
	public void setScript(Transcript script) {
		this.script = script;
	}
	
	public int getScore() {
		return score;
	}
	
	public void setScore(int score) {
		this.score = score;
	}
}
