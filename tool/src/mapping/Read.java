package mapping;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import bed.Segment;
import genome.Genome;
import genome.Chromosome;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import main.BamMethod;

public class Read {

	private String id = null;
	private String seq = null;
	private Map<String, IntervalTree<Segment>> segs = null;
	private Read pair = null;
	
	public Read(String id) {
		this.id = id;
		segs = new HashMap<>();
	}
	
	public Read(String id, Read pair) {
		this.id = id;
		this.pair = pair;
		segs = new HashMap<>();
	}

	public Read getPair() {
		return pair;
	}

	public void setPair(Read pair) {
		this.pair = pair;
	}

	public String getID() {
		return id;
	}

	public String getSeq() {
		return seq;
	}
	
	public Map<String, IntervalTree<Segment>> getSegs() {
		return segs;
	}
	
	public boolean addSegs(SAMRecord record) {
		return addSegs(record, 0, false, null);
	}
	
	public boolean addSegs(SAMRecord record, int mapQ, boolean seqTreeFlag, Map<String, IntervalTree<MappingStat>> mapTable) {
		if (record.getReadUnmappedFlag() || record.getMappingQuality() < mapQ) {
			return false;
		}
		String chr = record.getReferenceName().intern();
		if (seq == null || seq.length() < record.getReadLength()) {
			seq = record.getReadString();
		}
		IntervalTree<Segment> segTree = BamMethod.cigarToSeq(record, seqTreeFlag);
		int leftClip = getLeftClipBaseNum(record);
		segs.computeIfAbsent(chr, k -> new IntervalTree<>());
		for (Node<Segment> segNode : segTree) {
			segNode.getValue().setPositive(!record.getReadNegativeStrandFlag());
			segs.get(chr).put(leftClip + segNode.getStart(), leftClip + segNode.getEnd(), segNode.getValue());
			if (mapTable != null) {
				mapTable.computeIfAbsent(chr, k -> new IntervalTree<>());
				MappingMethod.putMappingStat(mapTable.get(chr), segNode.getValue().getStart() + 1,
						segNode.getValue().getEnd(), false, false);
			}
		}
		return true;
	}

	public void getJunction(Map<String, IntervalTree<IntervalTree<String>>> junctions, Map<String, Map<Integer, Double>> baseCount, int dev, Genome genome) {
		segs.forEach((chr, segTree) -> {
			if (segTree.size() < 2) {
				return;
			}
			for (Node<Segment> segNode : segTree) {
				for (Iterator<Node<Segment>> downNodes = segTree.overlappers(segNode.getEnd() + 1 - dev, segNode.getEnd() + 1 + dev);
						downNodes.hasNext();) {
					Node<Segment> downNode = downNodes.next();
					if (downNode != segNode && Math.abs(downNode.getStart() - segNode.getEnd() - 1) <= dev
							&& segNode.getValue().isPositive() == downNode.getValue().isPositive()) {
						int start = segNode.getValue().getEnd();
						int end = downNode.getValue().getStart() + 1;
						if (isValidJunction(genome == null ? null : genome.getChr(chr), start, end, dev)) {
							junctions.computeIfAbsent(chr, k -> new IntervalTree<>());
//							StringBuilder count = new StringBuilder();
//							count.append(id);
////							count.append('|');
////							count.append(segNode.getValue().getSeq());
////							count.append('|');
////							count.append(SimMethod.reverseFliq(downNode.getValue().getSeq()));
//							count.append('|');
//							String stat = Read.getMatchScore(segNode.getValue().getSeq(), 
//									SimMethod.reverseFliq(downNode.getValue().getSeq()), 0);
//							count.append(Read.getMatchScore(segNode.getValue().getSeq(), 
//									SimMethod.reverseFliq(downNode.getValue().getSeq()), 0));
//							count.append('|');
//							count.append(Math.min(segNode.getValue().getLength(), downNode.getValue().getLength()));
//							count.append(',');
//							if (start > end) {
//								count.insert(0, '-');
//								start ^= end;
//								end ^= start;
//								start ^= end;
//							}
//							count.append(junctions.get(chr).put(start, end, null));
							if (start <= end) {
								IntervalTree<String> oldList = junctions.get(chr).put(segNode.getValue().getStart() + 1, segNode.getValue().getEnd(), null);
								oldList = oldList == null ? new IntervalTree<>() : oldList;
								String old = oldList.put(downNode.getValue().getStart() + 1, downNode.getValue().getEnd(), null);
								old = old == null ? id + "," : old + id + ",";
								oldList.put(downNode.getValue().getStart() + 1, downNode.getValue().getEnd(), old);
								junctions.get(chr).put(segNode.getValue().getStart() + 1, segNode.getValue().getEnd(), oldList);
							}
							else {
								IntervalTree<String> oldList = junctions.get(chr).put(downNode.getValue().getStart() + 1, downNode.getValue().getEnd(), null);
								oldList = oldList == null ? new IntervalTree<>() : oldList;
								String old = oldList.put(segNode.getValue().getStart() + 1, segNode.getValue().getEnd(), null);
								old = old == null ? id + "," : old + id + ",";
								oldList.put(segNode.getValue().getStart() + 1, segNode.getValue().getEnd(), old);
								junctions.get(chr).put(downNode.getValue().getStart() + 1, downNode.getValue().getEnd(), oldList);
							}
							
//							if ("NT_167214.1".equals(chr)) {
//								baseCount.computeIfAbsent(chr, k -> new TreeMap<>());
//								incBaseCount(baseCount.get(chr), segNode.getValue(), downNode.getValue(), stat);
//							}
						}
					}
				}
			}
		});
	}
	
	private boolean isValidJunction(Chromosome chr, int start, int end, int dev) {
		return chr == null || (!chr.hasClipSignal(start, end, dev) && chr.withinGene(start, end) && !chr.isExonBoundry(start, end));
	}

	private int getLeftClipBaseNum(SAMRecord record) {
		if (!record.getCigar().isLeftClipped()) {
			return 0;
		}
		if (record.getCigar().getCigarElement(0).getOperator().consumesReadBases()) {
			return 0;
		}
		return record.getCigar().getCigarElement(0).getLength();
	}
	
	private void incBaseCount(Map<Integer, Double> baseCount, Segment seq, Segment ref, String stat) {
		int seqIndex = 1;
		int refIndex = 0;
		double score = 1.0;
//		for (int i = 0; i < stat.length(); ++i) {
//			score += stat.charAt(i) == 'M' ? 1.0 : 0.0;
//		}
//		score /= Math.min(seq.getLength(), ref.getLength());
		for (int i = 0; i < stat.length(); ++i) {
			switch (stat.charAt(i)) {
			case 'D':
				--refIndex;
				break;

			case 'I':
				++seqIndex;
				break;
				
			case 'M':
				baseCount.put(seqIndex + seq.getStart(), baseCount.getOrDefault(seqIndex + seq.getStart(), 0.0) + score);
				baseCount.put(refIndex + ref.getEnd(), baseCount.getOrDefault(refIndex + ref.getEnd(), 0.0) + score);
			case 'X':
				--refIndex;
				++seqIndex;
				break;
				
			default:
				break;
			}
		}
	}
	
	/**
	 * 
	 * @param seq sequence for match
	 * @param ref reference sequence
	 * @param method 0 means local; 1 means semi-global; 2 means global
	 * @return match stat sequence: 'D' means only in ref; 'I' means only in seq; 'X' means mismatch; 'M' means match
	 */
	public static String getMatchScore(String seq, String ref, int method, int[] statScore) {
		statScore = statScore == null || statScore.length < 4 ? new int[] {1, -1, -2, -2} : statScore;
		char[] stat = {'M', 'X', 'D', 'I'};
		int[][][] scoreMatrix = new int[ref.length() + 1][seq.length() + 1][2];
		method = Math.max(0, Math.min(2, method));
		for (int i = 1; i < scoreMatrix.length; ++i) {
			scoreMatrix[i][0][1] = 2;
			scoreMatrix[i][0][0] = method > 1 ? i * statScore[2] : 0;
		}
		int maxRef = 0;
		int maxSeq = 0;
		for (int j = 1; j <= seq.length(); ++j) {
			scoreMatrix[0][j][1] = 3;
			scoreMatrix[0][j][0] = method > 0 ? j * statScore[3] : 0;
			for (int i = 1; i < scoreMatrix.length; ++i) {
				int[] score = {scoreMatrix[i - 1][j - 1][0] + statScore[0], scoreMatrix[i - 1][j - 1][0] + statScore[1], 
						scoreMatrix[i - 1][j][0] + statScore[2], scoreMatrix[i][j - 1][0] + statScore[3], 0};
				int maxPath = seq.charAt(j - 1) == ref.charAt(i - 1) ? 0 : 1;
				maxPath = score[maxPath] < score[3] ? 3 : maxPath;
				maxPath = score[maxPath] < score[2] ? 2 : maxPath;
				maxPath = method == 0 && score[maxPath] < score[4] ? 4 : maxPath;
				scoreMatrix[i][j][1] = maxPath;
				scoreMatrix[i][j][0] = score[maxPath];
				if (scoreMatrix[maxRef][maxSeq][0] <= score[maxPath]) {
					maxRef = i;
					maxSeq = j;
				}
			}
		}
		
		maxRef = method == 2 ? ref.length() : maxRef;
		maxSeq = method > 0 ? seq.length() : maxSeq;
		if (method == 1) {
			for (int i = 0; i <= ref.length(); ++i) {
				maxRef = scoreMatrix[maxRef][seq.length()][0] < scoreMatrix[i][seq.length()][0] ? i : maxRef;
			}
		}
		
		char[] matchStat = new char[ref.length() + seq.length()];
		for (int i = maxSeq; i < seq.length(); ++i) {
			matchStat[i + maxRef] = stat[3];
		}
		for (int i = maxRef; i < ref.length(); ++i) {
			matchStat[seq.length() + i] = stat[2];
		}
		int maxScore = scoreMatrix[maxRef][maxSeq][0];
		while (maxRef > 0 || maxSeq > 0) {
			int statCode = scoreMatrix[maxRef][maxSeq][1];
			switch (statCode) {
			case 0:
			case 1:
				matchStat[--maxRef + --maxSeq] = stat[statCode];
				break;
				
			case 2:
				matchStat[--maxRef + maxSeq] = stat[statCode];
				break;
				
			case 3:
				matchStat[maxRef + --maxSeq] = stat[statCode];
				break;

			case 4:
				while (maxSeq > 0) {
					matchStat[maxRef + --maxSeq] = stat[3];
				}
				while (maxRef > 0) {
					matchStat[--maxRef] = stat[2];
				}
				
				break;
				
			default:
				break;
			}
		}
		int maxMatch = 0;
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < matchStat.length; ++i) {
			if (matchStat[i] > 0) {
				sb.append(matchStat[i]);
				maxMatch += matchStat[i] == stat[0] ? 1 : 0;
			}
		}
		sb.append('|');
		sb.append(maxScore);
		sb.append('|');
		sb.append(maxMatch);
		return sb.toString();
	}
}
