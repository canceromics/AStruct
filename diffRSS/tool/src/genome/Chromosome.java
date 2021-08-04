package genome;

import java.util.HashMap;
import java.util.Iterator;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

public class Chromosome extends IntRegion {

	public static final String[] CLIPPER = {"GT", "AG", "CT", "AC"}; 
	
	private String chr_symbol = null;
	private IntervalTree<Gene> gene_tree = null;
	private IntervalTree<Transcript> script_tree = null;
	private boolean need_new_script_tree = true;
	private IntervalTree<Exon> exon_tree = null;
	private boolean need_new_exon_tree = true;
	private String seq = null;
	private double maf = 0.0;
	
	public Chromosome(int length, String chr) {
		super(0, length);
		this.chr_symbol = chr;
	}

	public String getChr() {
		return chr_symbol;
	}
	
	public String getSeq() {
		return seq;
	}
	
	public void setSeq(String seq) {
		this.seq = seq;
		resetStartAndEnd(0, seq.length());
	}
	
	public double getMaf() {
		return maf;
	}
	
	public void setMaf(double maf) {
		this.maf = maf;
	}
	
	public IntervalTree<Gene> getGeneTree() {
		return gene_tree == null ? gene_tree = new IntervalTree<>() : gene_tree;
	}
	
	public IntervalTree<Exon> getExonTree() {
		return need_new_exon_tree ? buildExonTree() : exon_tree;
	}
	
	void setNewExonTree() {
		this.need_new_exon_tree = true;
	}
	
	public IntervalTree<Transcript> getScrpitTree() {
		return need_new_script_tree ? buildScriptTree() : script_tree;
	}
	
	void setNewScriptTree() {
		this.need_new_script_tree = true;
	}
	
	public HashMap<String, Gene> getGeneMap() {
		HashMap<String, Gene> out = new HashMap<>();
		for (Node<Gene> gene_node : getGeneTree()) {
			for (Gene gene = gene_node.getValue(); gene != null; gene = gene.getAnotherGene()) {
				out.put(gene.getGene_id(), gene);
			}
		}
		return out;
	}
	
	public HashMap<String, Transcript> getScriptMap() {
		HashMap<String, Transcript> out = new HashMap<>();
		for (Node<Gene> gene_node : getGeneTree()) {
			for (Gene gene = gene_node.getValue(); gene != null; gene = gene.getAnotherGene()) {
				out.putAll(gene.getScriptMap());
			}
		}
		return out;
	}
	
	private IntervalTree<Exon> buildExonTree() {
		exon_tree = new IntervalTree<>();
		Iterator<Node<Gene>> nodes = getGeneTree().iterator();
		while (nodes.hasNext()) {
			Node<Gene> node = nodes.next();
			node.getValue().addExonInTree(exon_tree);
		}
		need_new_exon_tree = false;
		return exon_tree;
	}
	
	private IntervalTree<Transcript> buildScriptTree() {
		if (gene_tree == null) {
			return null;
		}
		script_tree = new IntervalTree<>();
		Iterator<Node<Gene>> nodes = getGeneTree().iterator();
		while (nodes.hasNext()) {
			Node<Gene> node = nodes.next();
			node.getValue().addScriptInTree(script_tree);
		}
		need_new_script_tree = false;
		return script_tree;
	}
	
	public void sort() {
		for (Node<Gene> gene_node : getGeneTree()) {
			Gene g = new Gene(null, null, null, null, '.', 0, 0);
			g.setAnotherGene(gene_node.getValue());
			for (Gene gene = g; gene.getAnotherGene() != null; ) {
				if (!gene.getAnotherGene().sort()) {
					gene.setAnotherGene(gene.getAnotherGene().getAnotherGene());
				}
				else {
					gene = gene.getAnotherGene();
				}
			}
			if (g.getAnotherGene() == null) {
				getGeneTree().remove(gene_node.getStart(), gene_node.getEnd());
			}
			else {
				gene_node.setValue(g.getAnotherGene());
			}
		}
		Node<Gene> node = getGeneTree().find(0, 0);
		if (node != null) {
			for (Gene gene = node.getValue(); gene != null;) {
				Gene next = gene.getAnotherGene();
				gene.setAnotherGene(getGeneTree().put(gene.getStart() + 1, gene.getEnd(), gene));
				gene = next;
			}
			getGeneTree().remove(0, 0);
		}
	}
	
	public boolean hasClipSignal(int start, int end, int dev) {
		return hasClipSignal(start, end, dev, CLIPPER);
	}
	
	public boolean hasClipSignal(int start, int end, int dev, String[] clippers) {
		if (seq == null) {
			return false;
		}
		for (int i = 1; i < clippers.length; i += 2) {
			String upClip = clippers[i - 1];
			String downClip = clippers[i];
			String upSeq = seq.substring(Math.max(0, start - dev), Math.min(seq.length(), start + 2 + dev));
			String downSeq = seq.substring(Math.max(0, end - 3 - dev), Math.min(seq.length(), end - 1 + dev));
			if (upSeq.contains(upClip) && downSeq.contains(downClip)) {
				return true;
			}
		}
		return false;
	}
	
	public boolean withinGene(int start, int end) {
		for (Iterator<Node<Gene>> geneNodes = getGeneTree().overlappers(start, start); geneNodes.hasNext();) {
			Node<Gene> geneNode = geneNodes.next();
			if (end <= geneNode.getEnd() && end >= geneNode.getStart()) {
				return true;
			}
		}
		return false;
	}
	
	public boolean isExonBoundry(int start, int end) {
		boolean exonBound = false;
		for (Iterator<Node<Exon>> exonNodes = getExonTree().overlappers(start, start); exonNodes.hasNext();) {
			Node<Exon> exonNode = exonNodes.next();
			if (exonNode.getEnd() == start) {
				exonBound = true;
				break;
			}
		}
		if (!exonBound) {
			return false;
		}
		for (Iterator<Node<Exon>> exonNodes = getExonTree().overlappers(end, end); exonNodes.hasNext();) {
			Node<Exon> exonNode = exonNodes.next();
			if (exonNode.getStart() == end) {
				return true;
			}
		}
		return false;
	}
	
	@Override
	public String toString() {
		return chr_symbol;
	}
}
