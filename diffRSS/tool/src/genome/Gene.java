package genome;

import java.util.HashMap;

import htsjdk.samtools.util.IntervalTree;

public class Gene extends IntRegion{
	
	private Chromosome chr = null;
	private String gene_id = null;
	private String gene_symbol = null;
	private String gene_type = null;
	private char strand = '.';
	private Transcript script = null;
	private Gene same_position = null;
	private double sim_maf = 1.0;
	private double maf = 1.0;
	private double ip_enrich = 0.0;
	private int input_read_count = 0;
	private int ip_read_count = 0;
	private int peak_count = 0;
	
	public Gene(Chromosome chr, String gene_id, String gene_symbol, String gene_type, char strand, int start, int end) {
		super(start, end);
		this.chr = chr;
		this.gene_symbol = gene_symbol;
		this.gene_id = gene_id;
		this.gene_type = gene_type;
		this.strand = strand;
	}

	public String getChr() {
		return chr == null ? null : chr.getChr();
	}
	
	public Chromosome getChromosome() {
		return chr;
	}
	
	public String getGene_id() {
		return gene_id;
	}

	public String getGene_symbol() {
		return gene_symbol;
	}

	public String getGene_type() {
		return gene_type;
	}
	
	public char getStrand() {
		return strand;
	}
	
	public void setStrand(char strand) {
		this.strand = strand;
	}
	
	public Gene getAnotherGene() {
		return same_position;
	}
	
	public void setAnotherGene(Gene gene) {
		this.same_position = gene;
	}
	
	public double getMAF() {
		return maf;
	}
	
	public void setMAF(double maf) {
		this.maf = maf;
	}
	
	public double getSimMaf() {
		return sim_maf;
	}
	
	public void setSimMaf(double sim_maf) {
		this.sim_maf = sim_maf;
	}
	
	public int getInputReadCount() {
		return input_read_count;
	}

	public void setInputReadCount(int read_count) {
		this.input_read_count = read_count;
	}
	
	public int getIpReadCount() {
		return ip_read_count;
	}

	public void setIpReadCount(int ip_read_count) {
		this.ip_read_count = ip_read_count;
	}
	
	public double getIp_enrich() {
		return ip_enrich;
	}

	public void setIp_enrich(double ip_enrich) {
		this.ip_enrich = ip_enrich;
	}

	public int getPeakCount() {
		return peak_count;
	}
	
	public void setPeakCount(int peak_count) {
		this.peak_count = peak_count;
	}
	
	public int incPeakCount() {
		return ++peak_count;
	}
	
	public boolean hasPeak() {
		return peak_count > 0;
	}
	
	public void addAnotherGene(Gene gene) {
		if (this == gene) {
			return;
		}
		if (gene != null && this.same_position != null) {
			Gene g = gene;
			for (; g.same_position != null; g = g.same_position) {}
			g.same_position = this.same_position;
		}
		this.same_position = gene;
	}
	
	public Transcript getScript() {
		return script;
	}
	
	public void addScript(Transcript script) {
		if (script != null) {
			script.addGeneScript(this.script);
			this.script = script;
		}
	}
	
	public void mergeGene(Gene gene) {
		this.gene_id = this.gene_id + "," + gene.gene_id;
		this.gene_symbol = this.gene_symbol + ',' + gene_symbol;
		this.strand = this.strand == gene.strand ? this.strand : '.';
		this.script.addGeneScript(gene.script);
		resetStartAndEnd(Math.min(getStart(), gene.getStart()), Math.max(getEnd(), gene.getEnd()));
	}
	
	public HashMap<String, Transcript> getScriptMap() {
		HashMap<String, Transcript> out = new HashMap<>();
		for (Transcript script = this.script; script != null; script = script.getGeneScript()) {
			out.put(script.getScript_id(), script);
		}
		return out;
	}
	
	public void addScriptInTree(IntervalTree<Transcript> script_tree) {
		if (script_tree != null) {
			int start = Integer.MAX_VALUE;
			int end = Integer.MIN_VALUE;
			for (Transcript script = this.script; script != null; script = script.getGeneScript()) {
				script.resetStartAndEnd();
				script.addAnotherScript(script_tree.put(script.getStart() + 1, script.getEnd(), script));
				start = Math.min(start, script.getStart());
				end = Math.max(end, script.getEnd());
			}
			if (getStart() == getEnd()) {
				resetStartAndEnd(start, end);
			}
		}
	}
	
	public void addExonInTree(IntervalTree<Exon> exon_tree) {
		if (exon_tree != null) {
			int start = Integer.MAX_VALUE;
			int end = Integer.MIN_VALUE;
			for (Transcript script = this.script; script != null; script = script.getGeneScript()) {
				script.addExonInTree(exon_tree);
				start = Math.min(start, script.getStart());
				end = Math.max(end, script.getEnd());
			}
			if (getStart() == getEnd()) {
				resetStartAndEnd(start, end);
			}
		}
	}
	
	public boolean sort() {
		Transcript t = new Transcript(null, null, 0, 0, null);
		t.setGeneScript(script);
		for (Transcript s = t; s.getGeneScript() != null;) {
			if (!s.getGeneScript().sort()) {
				s.setGeneScript(s.getGeneScript().getGeneScript());
			}
			else {
				s = s.getGeneScript();
			}
		}
		script = t.getGeneScript();
		resetStartAndEnd();
		return script != null;
	}
	
	public void resetStartAndEnd() {
		if (getStart() == getEnd()) {
			int start = Integer.MAX_VALUE;
			int end = Integer.MIN_VALUE;
			for (Transcript s = this.script; s != null; s = s.getGeneScript()) {
				s.resetStartAndEnd();
				start = Math.min(start, s.getStart());
				end = Math.max(end, s.getEnd());
			}
			resetStartAndEnd(start, end);
		}
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(getChr());
		sb.append('\t');
		sb.append(getStart());
		sb.append('\t');
		sb.append(getEnd());
		sb.append('\t');
		sb.append(gene_id);
		sb.append('\t');
		sb.append(strand);
		sb.append('\t');
		sb.append(String.format("%.3f", ip_enrich));
		sb.append('\t');
		sb.append(input_read_count);
		sb.append('\t');
		sb.append(ip_read_count);
		return sb.toString();
	}
}
