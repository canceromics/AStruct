package mapping;

import java.util.Iterator;

import bed.Segment;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import main.BamMethod;

public class Alignment {
	
	private Alignment next = null;
	private Segment seg = null;
	
	public Alignment(Segment seg) {
		this.seg = seg;
	}
	
	public Alignment(Segment seg, Alignment front) {
		this.seg = seg;
		if (front != null) {
			this.next = front.next;
			front.next = this;
		}
	}
	
	public Alignment getNext() {
		return next;
	}

	public void setNext(Alignment next) {
		this.next = next;
	}

	public Segment getSeg() {
		return seg;
	}

	public void setSeg(Segment seg) {
		this.seg = seg;
	}

	public void addAlignment(SAMRecord record) {
		IntervalTree<Segment> segs = BamMethod.cigarToSeq(record);
		int last_end = -1;
		int last_start = -1;
		StringBuilder sb = new StringBuilder();
		for (Iterator<Node<Segment>> seg_nodes = segs.overlappers(0, Integer.MAX_VALUE); seg_nodes.hasNext();) {
			Node<Segment> seg_node = seg_nodes.next();
			if (last_end < seg_node.getStart()){
				if (last_end > 0) {
					new Alignment(new Segment(last_start, last_end - 1, sb.toString()), this);
				}
				sb = new StringBuilder();
				last_start = seg_node.getStart() - 1;
			}
			if (seg_node.getValue().getSeq() == null) {
				for (int i = seg_node.getStart(); i <= seg_node.getEnd(); ++i) {
					sb.append('D');
				}
			}
			else {
				sb.append(seg_node.getValue().getSeq());
			}
			last_end = seg_node.getEnd() + 1;
		}
		if (last_end > 0) {
			new Alignment(new Segment(last_start, last_end - 1, sb.toString()), this);
		}
	}
	
	public void clear() {
		seg = null;
		for (Alignment align = next; align != null; ) {
			Alignment tmp = align.next;
			align.next = null;
			align.seg = null;
			align = tmp;
		}
		next = null;
	}
}
