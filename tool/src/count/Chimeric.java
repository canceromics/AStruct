package count;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import genome.Genome;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import main.CommonMethod;
import main.InParam;
import main.Method;
import mapping.MappingStat;
import mapping.Read;

public class Chimeric {
	
	private int minOver = 5;
	private int minGap = 10;
	private int maxArm = 30;
	
	private static Chimeric chim = new Chimeric();
	private InParam args = null;
	
	private Chimeric() {
	}
	
	public static void run() {
		chim.args = InParam.getParams();
		try {
			Genome genome = Method.getInstance().loadGenomeInfo(chim.args.getExonFile(), chim.args.getGenomeFile());
			Map<String, Map<Integer, Double>> baseCount = new HashMap<>();
			Map<String, IntervalTree<MappingStat>> mapTable = new HashMap<>();
			Map<String, IntervalTree<IntervalTree<String>>> junctions = chim.getJunctionInBam(chim.args.getTreatFile(), genome, baseCount, mapTable);
			junctions = chim.mergeJunc(junctions);
			junctions = chim.mergeDG(junctions);
			chim.dgScore(junctions);
			chim.ngTag(junctions);
			
//			String chr = "NT_167214.1";
//			int start = 109078;
//			int end = 110946;
//			int[] juncCounts = chim.getJunctionCounts(start, end, junctions.get(chr));
//			int[] baseCounts = chim.getBaseCounts(start, end, mapTable.get(chr));
//			int[] matchCounts = chim.getMatchCounts(start, end, baseCount.get(chr));
//			double[] juncScore = RiboSNP.getInstance().getWinShapeScores(new int[][] {juncCounts, baseCounts}, start, end, null);
//			double[] matchScore = RiboSNP.getInstance().getWinShapeScores(new int[][] {matchCounts, baseCounts}, start, end, null);
//			System.out.println(Arrays.toString(juncScore));
//			System.out.println(Arrays.toString(matchScore));
			
			Method.writeFile(chim.args.getOutPrefix() + "_dg.txt", chim.dgOutput(junctions), null);
//			Method.writeFile(chim.args.getOutPrefix() + "_junc.txt", junctions, true, true, null);
//			Method.writeFile(chim.args.getOutPrefix() + "_basematch.txt", baseCount.get(chr).entrySet(), null);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private Map<String, IntervalTree<IntervalTree<String>>> getJunctionInBam(String file, Genome genome, 
			Map<String, Map<Integer, Double>> baseCount, Map<String, IntervalTree<MappingStat>> mapTable) throws IOException {
		Map<String, IntervalTree<IntervalTree<String>>> out = new HashMap<>();
		try (SamReader reader = SamReaderFactory.makeDefault().open(new File(file))) {
			Iterator<SAMRecord> rit = reader.iterator();
			Read read = null;
			while (rit.hasNext()) {
				SAMRecord record = rit.next();
				String id = record.getReadName();
				if (read == null || !id.equals(read.getID())) {
					if (read != null) {
						read.getJunction(out, baseCount, 0, genome);
						if (read.getPair() != null) {
							read.getPair().getJunction(out, baseCount, 0, genome);
						}
					}
					read = new Read(id);
				}
				if (record.isSecondaryOrSupplementary()) {
					continue;
				}
				if (record.getReadPairedFlag() && record.getSecondOfPairFlag()) {
					if (read.getPair() == null) {
						read.setPair(new Read(id, read));
					}
					read.getPair().addSegs(record, 0, true, mapTable);
				}
				else {
					read.addSegs(record, 0, true, mapTable);
				}
			}
			if (read != null) {
				read.getJunction(out, baseCount, 0, genome);
				if (read.getPair() != null) {
					read.getPair().getJunction(out, baseCount, 0, genome);
				}
			}
		}
		return out;
	}
	
	private Map<String, IntervalTree<IntervalTree<String>>> mergeJunc(Map<String, IntervalTree<IntervalTree<String>>> junctions){
		Map<String, IntervalTree<IntervalTree<String>>> dgroup = new HashMap<>();
		junctions.forEach((chr, juncTree) -> {
			for (Iterator<Node<IntervalTree<String>>> juncNodes = juncTree.overlappers(0, Integer.MAX_VALUE);
					juncNodes.hasNext();) {
				Node<IntervalTree<String>> juncNode = juncNodes.next();
				dgroup.computeIfAbsent(chr, k -> new IntervalTree<>());
				mergeJuncToDG(dgroup.get(chr), juncNode);
			}
		});
		return dgroup;
	}
	
	private void mergeJuncToDG(IntervalTree<IntervalTree<String>> dgTree, Node<IntervalTree<String>> juncTreeNode) {
		for (Node<String> juncNode : juncTreeNode.getValue()) {
			mergeJuncToDG(dgTree, juncTreeNode, juncNode);
		}
	}
	
	private void mergeJuncToDG(IntervalTree<IntervalTree<String>> dgTree, Node<IntervalTree<String>> juncTreeNode, Node<String> juncNode) {
		Node<String> dgNode = null;
		Node<IntervalTree<String>> dgTreeNode = null;
		for (Iterator<Node<IntervalTree<String>>> dgTNodes = dgTree.overlappers(juncTreeNode.getStart(), juncTreeNode.getEnd());
				dgTNodes.hasNext();) {
			Node<IntervalTree<String>> dgTNode = dgTNodes.next();
			if (isMinOverlap(dgTNode, juncTreeNode)) {
				for (Iterator<Node<String>> dgNodes = dgTNode.getValue().overlappers(juncNode.getStart(), juncNode.getEnd());
						dgNodes.hasNext();) {
					Node<String> tmpNode = dgNodes.next();
					if (isMinOverlap(tmpNode, juncNode)) {
						dgNode = tmpNode;
						dgTreeNode = dgTNode;
						break;
					}
				}
			}
		}
		if (dgNode != null) {
			boolean upSame = true;
			boolean downSame = true;
			int upStart = Math.max(dgTreeNode.getStart(), juncTreeNode.getStart());
			int upEnd = Math.min(dgTreeNode.getEnd(), juncTreeNode.getEnd());
			int downStart = Math.max(dgNode.getStart(), juncNode.getStart());
			int downEnd = Math.min(dgNode.getEnd(), juncNode.getEnd());
			upSame = upStart == dgTreeNode.getStart() && upEnd == dgTreeNode.getEnd();
			downSame = downStart == dgNode.getStart() && downEnd == dgNode.getEnd();
			
			if (upSame && downSame) {
				dgNode.setValue(dgNode.getValue() + juncNode.getValue());
			}
			else if (upSame) {
				dgTreeNode.getValue().remove(dgNode.getStart(), dgNode.getEnd());
				String old = dgTreeNode.getValue().put(downStart, downEnd, null);
				old = old == null ? dgNode.getValue() + juncNode.getValue() : old + dgNode.getValue() + juncNode.getValue();
				dgTreeNode.getValue().put(downStart, downEnd, old);
			}
			else {
				dgTreeNode.getValue().remove(dgNode.getStart(), dgNode.getEnd());
				if (dgTreeNode.getValue().size() <= 0) {
					dgTree.remove(dgTreeNode.getStart(), dgTreeNode.getEnd());
				}
				IntervalTree<String> oldTree = dgTree.put(upStart, upEnd, null);
				oldTree = oldTree == null ? new IntervalTree<>() : oldTree;
				String old = oldTree.put(downStart, downEnd, null);
				old = old == null ? dgNode.getValue() + juncNode.getValue() : old + dgNode.getValue() + juncNode.getValue();
				oldTree.put(downStart, downEnd, old);
				dgTree.put(upStart, upEnd, oldTree);
			}
		}
		else {
			IntervalTree<String> tmpTree = new IntervalTree<>();
			tmpTree.put(juncNode.getStart(), juncNode.getEnd(), juncNode.getValue());
			dgTree.put(juncTreeNode.getStart(), juncTreeNode.getEnd(), tmpTree);
		}
	}
	
	private <T> boolean isMinOverlap(Node<T> node1, Node<T> node2) {
		return Math.min(node1.getEnd(), node2.getEnd()) - Math.max(node1.getStart(), node2.getStart()) + 1 >= minOver;
	}
	
	private Map<String, IntervalTree<IntervalTree<String>>> mergeDG(Map<String, IntervalTree<IntervalTree<String>>> rawDG){
		Map<String, IntervalTree<IntervalTree<String>>> dgroup = new HashMap<>();
		rawDG.forEach((chr, juncTree) -> {
			for (Iterator<Node<IntervalTree<String>>> juncNodes = juncTree.overlappers(0, Integer.MAX_VALUE);
					juncNodes.hasNext();) {
				Node<IntervalTree<String>> juncNode = juncNodes.next();
				dgroup.computeIfAbsent(chr, k -> new IntervalTree<>());
				mergeDG(dgroup.get(chr), juncNode);
			}
		});
		return dgroup;
	}
	
	private void mergeDG(IntervalTree<IntervalTree<String>> dgTree, Node<IntervalTree<String>> juncTreeNode) {
		for (Node<String> juncNode : juncTreeNode.getValue()) {
			mergeDG(dgTree, juncTreeNode, juncNode);
		}
	}
	
	private void mergeDG(IntervalTree<IntervalTree<String>> dgTree, Node<IntervalTree<String>> juncTreeNode, Node<String> juncNode) {
		Node<String> dgNode = null;
		Node<IntervalTree<String>> dgTreeNode = null;
		for (Iterator<Node<IntervalTree<String>>> dgTNodes = dgTree.overlappers(juncTreeNode.getStart(), juncTreeNode.getEnd());
				dgTNodes.hasNext();) {
			Node<IntervalTree<String>> dgTNode = dgTNodes.next();
			if (canMerge(dgTNode, juncTreeNode)) {
				for (Iterator<Node<String>> dgNodes = dgTNode.getValue().overlappers(juncNode.getStart(), juncNode.getEnd());
						dgNodes.hasNext();) {
					Node<String> tmpNode = dgNodes.next();
					if (canMerge(tmpNode, juncNode)) {
						dgNode = tmpNode;
						dgTreeNode = dgTNode;
						break;
					}
				}
			}
		}
		if (dgNode != null) {
			boolean upSame = true;
			boolean downSame = true;
			int upStart = Math.min(dgTreeNode.getStart(), juncTreeNode.getStart());
			int upEnd = Math.max(dgTreeNode.getEnd(), juncTreeNode.getEnd());
			int downStart = Math.min(dgNode.getStart(), juncNode.getStart());
			int downEnd = Math.max(dgNode.getEnd(), juncNode.getEnd());
			upSame = upStart == dgTreeNode.getStart() && upEnd == dgTreeNode.getEnd();
			downSame = downStart == dgNode.getStart() && downEnd == dgNode.getEnd();
			
			if (upSame && downSame) {
				dgNode.setValue(dgNode.getValue() + juncNode.getValue());
			}
			else if (upSame) {
				dgTreeNode.getValue().remove(dgNode.getStart(), dgNode.getEnd());
				String old = dgTreeNode.getValue().put(downStart, downEnd, null);
				old = old == null ? dgNode.getValue() + juncNode.getValue() : old + dgNode.getValue() + juncNode.getValue();
				dgTreeNode.getValue().put(downStart, downEnd, old);
			}
			else {
				dgTreeNode.getValue().remove(dgNode.getStart(), dgNode.getEnd());
				if (dgTreeNode.getValue().size() <= 0) {
					dgTree.remove(dgTreeNode.getStart(), dgTreeNode.getEnd());
				}
				IntervalTree<String> oldTree = dgTree.put(upStart, upEnd, null);
				oldTree = oldTree == null ? new IntervalTree<>() : oldTree;
				String old = oldTree.put(downStart, downEnd, null);
				old = old == null ? dgNode.getValue() + juncNode.getValue() : old + dgNode.getValue() + juncNode.getValue();
				oldTree.put(downStart, downEnd, old);
				dgTree.put(upStart, upEnd, oldTree);
			}
		}
		else {
			IntervalTree<String> tmpTree = new IntervalTree<>();
			tmpTree.put(juncNode.getStart(), juncNode.getEnd(), juncNode.getValue());
			dgTree.put(juncTreeNode.getStart(), juncTreeNode.getEnd(), tmpTree);
		}
	}
	
	private void dgScore(Map<String, IntervalTree<IntervalTree<String>>> dgTable) {
		IntervalTree<DGarm> upTree = new IntervalTree<>();
		IntervalTree<DGarm> downTree = new IntervalTree<>();
		dgTable.forEach((chr, dgTree) -> dgScore(dgTree, upTree, downTree));
	}

	private void dgScore(IntervalTree<IntervalTree<String>> dgTree, IntervalTree<DGarm> upTree,
			IntervalTree<DGarm> downTree) {
		dgTree.forEach(dgTreeNode -> putIn(upTree, downTree, dgTreeNode));
		dgTree.forEach(dgTreeNode -> {
			dgTreeNode.getValue().forEach(dgNode -> {
				int up_down = upTree.find(dgTreeNode.getStart(), dgTreeNode.getEnd()).getValue()
						.getGapCover(dgNode.getStart(), dgNode.getEnd());
				int upCover = 0;
				for (Iterator<Node<DGarm>> upNodes = upTree.overlappers(dgTreeNode.getStart(), dgTreeNode.getEnd());
						upNodes.hasNext();) {
					Node<DGarm> upNode = upNodes.next();
					upCover += upNode.getValue().getCover();
				}
				int downCover = 0;
				for (Iterator<Node<DGarm>> downNodes = downTree.overlappers(dgNode.getStart(), dgNode.getEnd());
						downNodes.hasNext();) {
					Node<DGarm> downNode = downNodes.next();
					downCover += downNode.getValue().getCover();
				}
				double score = (double) up_down / Math.sqrt((double) upCover * downCover);
				if (score <= 0.01) {
					dgTreeNode.getValue().remove(dgNode.getStart(), dgNode.getEnd());
					if (dgTreeNode.getValue().size() <= 0) {
						dgTree.remove(dgTreeNode.getStart(), dgTreeNode.getEnd());
					}
				}
				else {
					dgNode.setValue(dgNode.getValue() + String.format("\t%.4f", score));
				}
			});
		});
	}
	
	private void putIn(IntervalTree<DGarm> upTree, IntervalTree<DGarm> downTree, Node<IntervalTree<String>> dgTreeNode) {
		DGarm upArm = new DGarm();
		upTree.put(dgTreeNode.getStart(), dgTreeNode.getEnd(), upArm);
		
		dgTreeNode.getValue().forEach(dgNode -> {
			int cover = getSepNum(dgNode.getValue(), ",");
			upArm.put(dgNode.getStart(), dgNode.getEnd(), cover);
			
			DGarm downArm = downTree.put(dgNode.getStart(), dgNode.getEnd(), null);
			downArm = downArm == null ? new DGarm() : downArm;
			downTree.put(dgNode.getStart(), dgNode.getEnd(), downArm);
			downArm.put(dgTreeNode.getStart(), dgTreeNode.getEnd(), cover);
		});
	}
	
	private int getSepNum(String s, String sep) {
		int num = 0;
		int index = -sep.length();
		while ((index = s.indexOf(sep, index + sep.length())) != -1) {
			++num;
		}
		return num;
	}

	private <T> boolean canMerge(Node<T> node1, Node<T> node2) {
		boolean out = node1.getEnd() < node2.getEnd() ? 
				node2.getStart() - node1.getEnd() < minGap : node1.getStart() - node2.getEnd() < minGap;
		out &= Math.max(node1.getEnd(), node2.getEnd()) - Math.min(node1.getStart(), node2.getStart()) < maxArm;
		return out;
	}

	private void ngTag(Map<String, IntervalTree<IntervalTree<String>>> dgTable) {
		dgTable.forEach((chr, dgTree) -> {
			List<IntervalTree<Integer>> ngList = new ArrayList<>();
			dgTree.forEach(dgTreeNode -> {
				dgTreeNode.getValue().forEach(dgNode -> {
					IntervalTree<Integer> ngTree = null;
					for (int i = 0; i < ngList.size(); ++i) {
						if (!ngList.get(i).overlappers(dgTreeNode.getStart(), dgTreeNode.getEnd()).hasNext() 
								&& !ngList.get(i).overlappers(dgNode.getStart(), dgNode.getEnd()).hasNext()) {
							ngTree = ngList.get(i);
							break;
						}
					}
					if (ngTree == null) {
						ngTree = new IntervalTree<>();
						ngList.add(ngTree);
						ngTree.setSentinel(ngList.size());
					}
					ngTree.put(dgTreeNode.getStart(), dgTreeNode.getEnd(), 0);
					ngTree.put(dgNode.getStart(), dgNode.getEnd(), 0);
					dgNode.setValue(dgNode.getValue() + String.format("\tNG%d", ngTree.getSentinel()));
				});
			});
			if (!"NT_167214.1".equals(chr)) {
				return;
			}
			try {
				double[] stand18s = Method.getCol(Method.readMatrix("D:/work/office/chim_18s.csv", ",", 1), 1);
				double[] stand28s = Method.getCol(Method.readMatrix("D:/work/office/chim_28s.csv", ",", 1), 1);
				int start18s = 109078;
				int start28s = 113348;
				for (int i = 0; i < ngList.size(); ++i) {
					double[] score18s = new double[stand18s.length];
					double[] score28s = new double[stand28s.length];
					for (int j = 0; j < stand18s.length; ++j) {
						score18s[j] = ngList.get(i).overlappers(start18s + j, start18s + j).hasNext() ? 1.0 : 0.0;
					}
					for (int j = 0; j < stand28s.length; ++j) {
						score28s[j] = ngList.get(i).overlappers(start28s + j, start28s + j).hasNext() ? 1.0 : 0.0;
					}
					System.out.printf("NG%d:\t18s:%.3f\t28s:%.3f\n", i + 1, 
						CommonMethod.calAUC(CommonMethod.calROC(score18s, stand18s, 0.5)), 
						CommonMethod.calAUC(CommonMethod.calROC(score28s, stand28s, 0.5)));
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		});
	}
	
	private List<String> dgOutput(Map<String, IntervalTree<IntervalTree<String>>> dgTable) {
		List<String> out = new ArrayList<>();
		dgTable.forEach((chr, dgTree) ->{
			for (Iterator<Node<IntervalTree<String>>> dgTreeNodes = dgTree.overlappers(0, Integer.MAX_VALUE);
					dgTreeNodes.hasNext();) {
				Node<IntervalTree<String>> dgTreeNode = dgTreeNodes.next();
				for (Iterator<Node<String>> dgNodes = dgTreeNode.getValue().overlappers(0, Integer.MAX_VALUE);
						dgNodes.hasNext();) {
					Node<String> dgNode = dgNodes.next();
					StringBuilder sb = new StringBuilder();
					sb.append(chr);
					sb.append('\t');
					sb.append(dgTreeNode.getStart());
					sb.append('\t');
					sb.append(dgTreeNode.getEnd());
					sb.append('\t');
					sb.append(dgNode.getStart());
					sb.append('\t');
					sb.append(dgNode.getEnd());
					sb.append('\t');
					sb.append(dgNode.getValue());
					out.add(sb.toString());
				}
			}
		});
		return out;
	}
}
