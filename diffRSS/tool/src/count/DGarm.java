package count;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

public class DGarm {

	private int cover = 0;
	private IntervalTree<Integer> dgTree = null;

	public DGarm() {
		dgTree = new IntervalTree<>();
		dgTree.setSentinel(0);
	}
	
	public void put(int start, int end, int c) {
		int old = dgTree.put(start, end, 0);
		old += c;
		cover += c;
		dgTree.put(start, end, old);
	}
	
	public int getGapCover(int start, int end) {
		Node<Integer> n = dgTree.find(start, end);
		return n == null ? 0 : n.getValue();
	}
	
	public int getCover() {
		return cover;
	}
}
