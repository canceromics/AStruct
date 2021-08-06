package count;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

import bed.Mutation;
import bed.Peak;
import genome.Gene;
import genome.Genome;
import genome.IntRegion;
import genome.Transcript;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import main.CommonMethod;
import main.InParam;
import main.Method;
import mapping.Alignment;

public class RiboSNP {

	private static RiboSNP cm = new RiboSNP();
	private double sub_fac = 0.25;
	private double div_fac = 10.0;
	private double win_fac = 0.05;
	private boolean ip_flag = false;
	private boolean permutate = false;
	private boolean replicate = false;
	private int alignLen = 50;
	private int window = 50;
	private int winStep = 3;
	private int treatReads = 0;
	private int controlReads = 0;
	private List<String> chrs = null;
	private Map<String, IntervalTree<Mutation>> mutTable = new HashMap<>();
	private IntervalTree<Mutation> emptyMutTree = new IntervalTree<>();
	private List<Iterator<SAMRecord>> itsT = null;
	private List<Iterator<SAMRecord>> itsC = null;
	private List<SAMRecord> recordsT = null;
	private List<SAMRecord> recordsC = null;
	private IntervalTree<List<String>> treeT = null;
	private IntervalTree<List<String>> treeC = null;
	private IntervalTree<List<String>> treeTR1 = null;
	private IntervalTree<List<String>> treeTR2 = null;
	private IntervalTree<List<String>> treeTH = null;
	private IntervalTree<List<String>> treeCH = null;
	private InParam args = null;
	
	private RiboSNP() {
	}
	
	public static void run() {
		try {
			cm.args = InParam.getParams();
			Genome genome = Method.loadGenomeInfo(cm.args.getExonFile(), cm.args.getGenomeFile());
//			Genome genome = new Genome();
			if (!cm.args.getTreatFiles2().isEmpty()) {
				if (cm.args.getMutFile() != null) {
					cm.getSampleSnpinSortBam(genome);
				}
				else {
					cm.getSamplesinSortBam(genome);
				}
				return;
			}
//			Method.readMutation(cm.args.getMutFile(), cm.mutTable);
//			for (String chr : cm.mutTable.keySet()) {
//				Method.setStrictScript(cm.mutTable.get(chr), genome.getChr(chr));
//			}
			cm.getSNPinSortBam(cm.args.getTreatFiles(), cm.args.getControlFiles(), genome);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static RiboSNP getInstance() {
		return cm;
	}
	
	public static void main(String[] args) {
		try(BufferedReader testR = new BufferedReader(new FileReader("D:/work/workspace/data/Test_esdc.txt"))) {
//			SamReader readerT = SamReaderFactory.makeDefault().open(new File(args[0]));
//			SamReader readerC = SamReaderFactory.makeDefault().open(new File(args[1]));
//			BufferedReader reader = new BufferedReader(new FileReader(new File(args[2])))) {
//
			int[][] counts = {
					{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 2, 2, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 2, 0, 1, 2, 2, 2, 2, 1, 1, 1, 0, 0, 2, 0, 3, 0, 1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 3, 2, 0, 2, 3, 1, 0, 3, 2, 1, 1, 2, 1, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 0, 2, 1, 4, 0, 0, 0, 0, 2, 0, 0, 0, 0, 4, 0, 0, 0, 3, 0, 0, 0, 0, 1, 1, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 5, 0, 0, 0, 0, 0, 0, 2, 0, 0, 3, 3, 5, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, 0, 0, 0},
					{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 7, 9, 10, 11, 12, 12, 12, 13, 14, 15, 16, 16, 16, 16, 16, 18, 18, 19, 21, 23, 25, 27, 28, 29, 30, 30, 30, 32, 32, 35, 35, 36, 37, 37, 39, 39, 39, 39, 39, 39, 39, 39, 42, 44, 44, 46, 49, 50, 50, 53, 55, 56, 57, 59, 60, 60, 63, 65, 65, 65, 65, 65, 65, 65, 65, 67, 68, 69, 69, 71, 72, 76, 76, 76, 76, 76, 78, 78, 78, 78, 78, 82, 82, 82, 82, 85, 85, 85, 85, 85, 86, 87, 88, 90, 92, 93, 93, 93, 93, 93, 93, 93, 93, 93, 94, 99, 99, 99, 99, 99, 99, 99, 101, 101, 101, 104, 107, 112, 112, 112, 112, 112, 112, 112, 112, 112, 116, 117, 117, 117},
					{0, 1, 0, 0, 0, 4, 0, 0, 0, 3, 0, 0, 1, 3, 1, 0, 0, 2, 3, 3, 1, 1, 3, 1, 2, 2, 5, 2, 1, 0, 2, 1, 0, 2, 1, 1, 3, 0, 2, 1, 1, 3, 0, 1, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 3, 1, 1, 2, 0, 5, 2, 0, 2, 1, 1, 0, 1, 3, 1, 1, 0, 3, 0, 2, 2, 3, 0, 1, 0, 2, 0, 1, 1, 3, 3, 1, 1, 0, 2, 3, 3, 1, 1, 3, 1, 1, 0, 4, 2, 0, 2, 1, 1, 1, 2, 1, 3, 2, 0, 1, 2, 2, 1, 1, 1, 1, 0, 1, 0, 2, 0, 1, 1, 0, 5, 1, 3, 1, 5, 0, 0, 4, 1, 1, 0, 4, 2, 1, 0, 1, 0, 3, 1, 2, 2, 1, 0, 0, 0, 2, 2, 2, 5, 3, 1, 2, 0, 1, 0},
					{1, 1, 2, 2, 2, 2, 6, 6, 6, 6, 9, 9, 9, 10, 13, 14, 14, 14, 16, 19, 22, 23, 24, 27, 28, 30, 32, 37, 39, 40, 40, 42, 43, 43, 45, 46, 47, 50, 50, 52, 53, 54, 57, 57, 58, 61, 63, 64, 65, 66, 67, 69, 70, 71, 72, 75, 76, 77, 79, 79, 84, 86, 86, 88, 89, 90, 90, 91, 94, 95, 96, 96, 99, 99, 101, 103, 106, 106, 107, 107, 109, 109, 110, 111, 114, 117, 118, 119, 119, 121, 124, 127, 128, 129, 132, 133, 134, 134, 138, 140, 140, 142, 143, 144, 145, 147, 148, 151, 153, 153, 154, 156, 158, 159, 160, 161, 162, 162, 163, 163, 165, 165, 166, 167, 167, 172, 173, 176, 177, 182, 182, 182, 186, 187, 188, 188, 192, 194, 195, 195, 196, 196, 199, 200, 202, 204, 205, 205, 205, 205, 207, 209, 211, 216, 219, 220, 222, 222, 223}
			};
			System.out.println(Arrays.toString(cm.getWinShapeScores(counts, 1, counts[0].length, null)));
			
			double[][] dataMat = Method.readMatrix("D:/work/workspace/data/beforeNorm_check1.txt");
			List<double[]> rocPoints = CommonMethod.calROC(dataMat, 0.5);
			for (double[] point :rocPoints) {
				System.out.println(Arrays.toString(point));
			}
			Mutation mut = new Mutation('T', new char[] {'A', 'G'}, "rs1");
			for (int i = 0; i < 45; ++i) {
				if (i < 35) {
					mut.incCount('T', true);
				}
				if (i < 20) {
					mut.incCount('A', false);
				}
				if (i < 15) {
					mut.incCount('A', true);
				}
				if (i < 10) {
					mut.incCount('G', false);
				}
				if (i < 8) {
					mut.incCount('G', true);
				}
				mut.incCount('T', false);
			}
			System.out.println(mut.isValidSNP(false, 2));
			cm.args = InParam.getParams();
			IntervalTree<Mutation> snpTree = new IntervalTree<>();
			snpTree.put(5, 5, mut);
			if (cm.removeSNP(snpTree, snpTree.find(5, 5))) {
				System.out.println("debug");
			}
			mut.remakeMostAlt(false);
			
			String line = null;
			List<double[]> readL = new ArrayList<>();
			KolmogorovSmirnovTest kst = new KolmogorovSmirnovTest();
			String name = null;
			double[] scores1 = null;
			double[] scores2 = null;
			while ((line = testR.readLine()) != null) {
				String[] cols = line.split(",");
				if (cols.length == 1) {
					name = cols[0];
				}
				else {
					if (scores1 == null) {
						scores1 = new double[cols.length];
						for (int i = 0; i < cols.length; ++i) {
							scores1[i] = Double.parseDouble(cols[i]);
						}
					}
					else {
						scores2 = new double[cols.length];
						for (int i = 0; i < cols.length; ++i) {
							scores2[i] = Double.parseDouble(cols[i]);
						}
						System.out.printf("%s\t%.4f\t%.4f\t%.4f\n", name, CommonMethod.calKSTest(scores1, scores2),
								kst.kolmogorovSmirnovTest(scores1, scores2), cm.getEsdcScore(scores1, scores2));
						scores1 = null;
						scores2 = null;
					}
				}
//				readL.add(new double[] {Double.parseDouble(cols[0]), Double.parseDouble(cols[1])});
			}
//			scores1 = new double[readL.size()];
//			scores2 = new double[readL.size()];
//			for (int i = 0; i < readL.size(); ++i) {
//				scores1[i] = readL.get(i)[0];
//				scores2[i] = readL.get(i)[1];
//			}
//			System.out.println(kst.kolmogorovSmirnovStatistic(scores1, scores2));
//			System.out.println(kst.kolmogorovSmirnovTest(scores1, scores2));
//			System.out.println(Method.calKSTest(scores1, scores2));
//			System.out.println(Method.calKSTest(scores1, scores2));
			
//			List<String> output = new ArrayList<>();
//			
//			cm.bases = new HashSet<>();
//			cm.bases.add('A');
//			cm.bases.add('C');
//			cm.bases.add('T');
//			cm.bases.add('G');
////			int start = 109078;
////			int end = 110946;
////			int start = 113348;
////			int end = 118417;
//			int start = 108879;
//			int end = 118616;
//			
////			double[][] counts = new double[4][end - start + 1];
//			int [][] countsF = new int[4][end - start + 1];
//			boolean[] truth = new boolean[end - start + 1];
//			int index = 0;
//			while ((line = reader.readLine()) != null) {
//				String[] cols = line.split("\t");
//				
//				for (int i = 0; i < 4; ++i) {
//					countsF[i][index] = Integer.parseInt(cols[i + 3]);
//				}
//				truth[index] = "0".equals(cols[2]);
//				++index;
//			}
//			
//			double[] maxAuc = {0.0, 0.0};
//			double[] minAuc = {1.0, 1.0};
//			int[] maxAucWin = {0, 0, 0, 0};
//			int[] minAucWin = {0, 0, 0, 0};
//			double[] maxSubF = {0.00, 0.00};
//			double[] minSubF = {0.00, 0.00};
//			double[][][] allAuc18 = new double[39][10][21];
//			double[][][] allAuc28 = new double[39][10][21];
//			for (cm.window = 10; cm.window <= 200; cm.window += 5) {
//				for (cm.winStep = 1; cm.winStep <= 10; ++cm.winStep) {
//					for (cm.sub_fac = 0.00; cm.sub_fac < 1.005; cm.sub_fac += 0.05) {
//						double[] scoresF = cm.calScores(countsF, start, end);
//						double auc18 = Method.calAUC(Method.calROC(Arrays.copyOfRange(scoresF, 109078 - start, 110947 - start),
//								Arrays.copyOfRange(truth, 109078 - start, 110947 - start)));
//						double auc28 = Method.calAUC(Method.calROC(Arrays.copyOfRange(scoresF, 113348 - start, 118418 - start),
//								Arrays.copyOfRange(truth, 113348 - start, 118418 - start)));
//						allAuc18[cm.window / 5 - 2][cm.winStep - 1][(int) Math.round(cm.sub_fac * 20.0)] = auc18;
//						allAuc28[cm.window / 5 - 2][cm.winStep - 1][(int) Math.round(cm.sub_fac * 20.0)] = auc28;
//						if (maxAuc[0] < auc18) {
//							maxAuc[0] = auc18;
//							maxAucWin[0] = cm.window;
//							maxAucWin[1] = cm.winStep;
//							maxSubF[0] = cm.sub_fac;
//						}
//						if (minAuc[0] > auc18) {
//							minAuc[0] = auc18;
//							minAucWin[0] = cm.window;
//							minAucWin[1] = cm.winStep;
//							minSubF[0] = cm.sub_fac;
//						}
//						if (maxAuc[1] < auc28) {
//							maxAuc[1] = auc28;
//							maxAucWin[2] = cm.window;
//							maxAucWin[3] = cm.winStep;
//							maxSubF[1] = cm.sub_fac;
//						}
//						if (minAuc[1] > auc28) {
//							minAuc[1] = auc28;
//							minAucWin[2] = cm.window;
//							minAucWin[3] = cm.winStep;
//							minSubF[1] = cm.sub_fac;
//						}
//					}
//				}
//			}
//
//			System.out.println("18s:\tMax: " + maxAuc[0] + "\tWindow:\t" + maxAucWin[0] + "\tStep:\t" + maxAucWin[1] + "\tSub:\t" + maxSubF[0]);
//			System.out.println("18s:\tMin: " + minAuc[0] + "\tWindow:\t" + minAucWin[0] + "\tStep:\t" + minAucWin[1] + "\tSub:\t" + minSubF[0]);
//			System.out.println("28s:\tMax: " + maxAuc[1] + "\tWindow:\t" + maxAucWin[2] + "\tStep:\t" + maxAucWin[3] + "\tSub:\t" + maxSubF[1]);
//			System.out.println("28s:\tMin: " + minAuc[1] + "\tWindow:\t" + minAucWin[2] + "\tStep:\t" + minAucWin[3] + "\tSub:\t" + minSubF[1]);
//			for (int i = 0; i < 21; ++i) {
//				System.out.println(i);
//				for (int j = 0; j < 10; ++j) {
//					StringBuilder sb18 = new StringBuilder();
//					StringBuilder sb28 = new StringBuilder();
//					for (int k = 0; k < 39; ++k) {
//						sb18.append(String.format("%.3f", allAuc18[k][j][i]));
//						sb18.append(',');
//						sb28.append(String.format("%.3f", allAuc28[k][j][i]));
//						sb28.append(',');
//					}
//					sb18.setLength(sb18.length() - 1);
//					sb28.setLength(sb28.length() - 1);
//					sb18.append('\t');
//					sb18.append(sb28);
//					System.out.println(sb18);
//				}
//			}
//			cm.window = 50;
//			cm.winStep = 3;
//			cm.sub_fac = 0.25;
//			double[] scoresF = cm.calScores(countsF, start, end);
//			double auc28 = Method.calAUC(Method.calROC(Arrays.copyOfRange(scoresF, 113348 - start, 118418 - start),
//					Arrays.copyOfRange(truth, 113348 - start, 118418 - start)));
//			
//			
//			if (scoresF != null) {
//				for (int i = 0; i < scoresF.length; ++i) {
//					StringBuilder sb = new StringBuilder();
//					sb.append(i + start);
//					sb.append('\t');
//					sb.append(countsF[0][i]);
//					sb.append('\t');
//					sb.append(countsF[1][i]);
//					sb.append('\t');
//					sb.append(countsF[2][i]);
//					sb.append('\t');
//					sb.append(countsF[3][i]);
//					sb.append('\t');
//					sb.append(String.format("%.3f", scoresF[i]));
////						sb.append('\t');
////						sb.append(counts[0][i]);
////						sb.append('\t');
////						sb.append(counts[0][i]);
////						sb.append('\t');
////						sb.append(counts[0][i]);
//					output.add(sb.toString());
//				}
//			}
//		Method.writeFile(args[4], output, null);
//		output = new ArrayList<>();
//		start = Integer.parseInt(args[5]);
//		end = Integer.parseInt(args[6]);
//		Genome genome = Method.getInstance().loadGenomeInfo(null, args[3]);
////			for (int i = start; i <= end; ++i) {
////				if (mutBases.contains(Character.toUpperCase(genome.getChr("NT_167214.1").getSeq().charAt(i - 1)))) {
////					System.out.println(i);
////				}
////			}
//			
////			cm.sub_fac = 0.25;
////			for (int i = 0; i < counts.length; ++i) {
////				cm.normalize(counts[i]);
////			}			
////			double[] tscores = cm.getRTstopScroes(counts[1], counts[3], counts[0], counts[2]);
////			cm.normalize2(tscores);
////			for (double d : tscores) {
////				System.out.println(d);
////			}
//				if (!Method.isCoordinate(readerT.getFileHeader()) || !Method.isCoordinate(readerC.getFileHeader())) {
//					return;
//				}
//				cm.chrs = new ArrayList<>();
//				cm.chrs.add("NT_167214.1");
////				List<Chromosome> chrs = Method.getChrs(readerT.getFileHeader());
//				List<Iterator<SAMRecord>> itsT = new ArrayList<>(1);
//				List<Iterator<SAMRecord>> itsC = new ArrayList<>(1);
//				itsT.add(readerT.iterator());
//				itsC.add(readerC.iterator());
//				List<SAMRecord> recordsT = new ArrayList<>(1);
//				List<SAMRecord> recordsC = new ArrayList<>(1);
//				recordsT.add(itsT.get(0).hasNext() ? itsT.get(0).next() : null);
//				recordsC.add(itsC.get(0).hasNext() ? itsC.get(0).next() : null);
//				for (int chrIndex = 0; chrIndex < cm.chrs.size(); ++chrIndex) {
//					String chr = cm.chrs.get(chrIndex);
//					IntervalTree<List<String>> treeT = new IntervalTree<>();
//					IntervalTree<List<String>> treeC = new IntervalTree<>();
//					cm.ip_flag = true;
//					cm.addRegion(treeT, itsT, recordsT, chr, start, end);
//					if (args.length > 7) {
//						cm.ip_flag = false;
//						cm.addRegion(treeC, itsC, recordsC, chr, start, end);
//					}
//					int[][][][] counts = cm.getRTstopCounts(null, treeT, treeC, null, genome.getChr(chr).getSeq(), start, end, true);
//					cm.removeRegion(treeT, treeC, end + 1);
//					
//					IntervalTree<double[]> scoreTree = new IntervalTree<>();
//					int winStart = start;
//					double[] scores = new double[end - start + 1];
//					int scoreStart = 0;
//					for (int winEnd = start + cm.window - 1; winEnd < end; winEnd += cm.winStep) {
//						scoreTree.put(winStart, winEnd, cm.normalize2(cm.getRTstopScores(counts[0][0], winStart - start, winEnd - start + 1)));
//						int nextStart = Math.min(winEnd + cm.winStep, end) - cm.window + 1;
//						double[] winScore = cm.calScores(scoreTree, winStart, nextStart);
//						for (int i = 0; i < winScore.length; ++i, ++scoreStart) {
//							scores[scoreStart] = winScore[i];
//						}
//						winStart = nextStart;
//					}
//					scoreTree.put(winStart, end, cm.normalize2(cm.getRTstopScores(counts[0][0], winStart - start, end - start + 1)));
//					double[] winScore = cm.calScores(scoreTree, winStart, end + 1);
//					for (int i = 0; i < winScore.length; ++i, ++scoreStart) {
//						scores[scoreStart] = winScore[i];
//					}
////					double max_auc = 0.0;
////					double min_auc = 1.0;
////					double[] indexs = new double[2];
////					for (cm.sub_fac = 0.0; cm.sub_fac < 1.0; cm.sub_fac += 0.01) {
////						double[] scores = cm.getRTmutScores(treeT, treeC, mutBases, genome.getChr(chr).getSeq(), start, end, true);
////						cm.normalize2(scores);
////						double auc = Method.calAUC(Method.calROC(scores, truth));
////						if (auc < min_auc) {
////							indexs[0] = cm.sub_fac;
////							min_auc = auc;
////						}
////						if (auc > max_auc) {
////							indexs[1] = cm.sub_fac;
////							max_auc = auc;
////						}
////					}
////					cm.sub_fac = indexs[0];
////					cm.sub_fac = 0.25;
////					double[] scores = cm.getRTmutScores(treeT, mutBases, genome.getChr(chr).getSeq(), start, end, true);
////					double[] scores = cm.getRTstopScores(treeT, treeC, genome.getChr(chr).getSeq(), start, end, true);
////					cm.normalize2(scores);
////					for (double d : scores) {
////						System.out.println(Double.isFinite(d) ? d : -1);
////					}
//					
////					int[][] counts = cm.getRTstopCounts(treeT, treeC, genome.getChr(chr).getSeq(), start, end, true);
//					
////					int[][] counts = cm.getRTstopCounts(treeT, treeC, null, 1, chrs.get(chrIndex).getLength(), true);
////					double[] scores = cm.getTestRTstopScores(counts);
////					double[] normal = cm.normalize2(scores);
////					double[] baseScores = new double[4];
////					if (normal != null) {
//						for (int i = 0; i < scores.length; ++i) {
//							StringBuilder sb = new StringBuilder();
//							sb.append(chr);
//							sb.append('\t');
//							sb.append(i + start);
//							sb.append('\t');
//							sb.append(counts[0][i]);
//							sb.append('\t');
//							sb.append(counts[1][i]);
//							if (args.length > 7) {
//								sb.append('\t');
//								sb.append(counts[2][i]);
//								sb.append('\t');
//								sb.append(counts[3][i]);
//							}
//							sb.append('\t');
//							sb.append(String.format("%.3f", scores[i]));
////							sb.append('\t');
////							sb.append(String.format("%.3f", normal[i]));
//	//						sb.append('\t');
//	//						sb.append(counts[0][i]);
//	//						sb.append('\t');
//	//						sb.append(counts[0][i]);
//	//						sb.append('\t');
//	//						sb.append(counts[0][i]);
//							output.add(sb.toString());
//						}
////					}
//				}
//				Method.writeFile(args[4], output, null);
			} catch (IOException e) {
				e.printStackTrace();
			}
	}
	
//	public static void countSNP() {
//		try {
//			Genome genome = Method.getInstance().loadGenomeInfo(InParam.getParams().getExonFile(), InParam.getParams().getGenomeFile());
//			HashMap<String, IntervalTree<Mutation>> mut_table = Method.readMutation(InParam.getParams().getMutFile(), null);
//			mut_table.forEach((chr, mutTree) -> {
//				if (genome.getChr(chr) == null) {
//					return;
//				}
//				Method.setStrictScript(mutTree, genome.getChr(chr));
//			});
//			cm.getSNPinSortBam(InParam.getParams().getTreatFile(), InParam.getParams().getControlFile(), genome, mut_table);
////			Method.writeFile(InParam.getParams().getOutPrefix() + "_mut.txt", Method.getInstance().mutOutput(genome, mut_table), Mutation.getSimpleHeader());
//			
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//	}
	
	private void addSeg(IntervalTree<List<String>> align_tree, IntervalTree<List<String>> repTree, IntervalTree<Mutation> mutTree, Alignment align, int spliceStat) {
		for (align = align.getNext(); align != null; align = align.getNext()) {
			if (align.getSeg() == null) {
				continue;
			}
			List<String> old = align_tree.put(align.getSeg().getStart() + 1, align.getSeg().getEnd(), null);
			old = old == null ? new ArrayList<>(1) : old;
			String seq = align.getSeg().getSeq();
			seq = spliceStat == 1 ? "S" + seq : spliceStat == 2 ? seq + "S" : spliceStat == 3 ? "S" + seq + "S" : seq;
			old.add(seq);
			align_tree.put(align.getSeg().getStart() + 1, align.getSeg().getEnd(), old);
			if (repTree != null) {
				old = repTree.put(align.getSeg().getStart() + 1, align.getSeg().getEnd(), null);
				old = old == null ? new ArrayList<>(1) : old;
				old.add(seq);
				repTree.put(align.getSeg().getStart() + 1, align.getSeg().getEnd(), old);
			}
			for (Iterator<Node<Mutation>> mutNodes = mutTree.overlappers(align.getSeg().getStart() + 1, align.getSeg().getEnd());
				mutNodes.hasNext();) {
				Node<Mutation> mutNode = mutNodes.next();
				char key = align.getSeg().getSeq().charAt(mutNode.getStart() - align.getSeg().getStart() - 1);
				mutNode.getValue().incCount(key, ip_flag);
			}
		}
	}
	
	private void addRegion(IntervalTree<List<String>> alignTree, List<Iterator<SAMRecord>> its, 
			List<SAMRecord> records, String chr, int start, int end) {
		int mid = (start + end) / 2;
		for (int i = 0; i < its.size(); ++i) {
			while (records.get(i) != null && chrs.indexOf(records.get(i).getReferenceName()) < chrs.indexOf(chr)) {
				if (alignLen < records.get(i).getReadLength()) {
					alignLen = records.get(i).getReadLength();
					start = mid - alignLen + 1;
					end = mid + alignLen - 1;
				}
				records.set(i, its.get(i).hasNext() ? its.get(i).next() : null);
			}
			while (records.get(i) != null && sameChr(records.get(i), chr) && records.get(i).getAlignmentStart() <= end) {
				if (alignLen < records.get(i).getReadLength()) {
					alignLen = records.get(i).getReadLength();
					start = mid - alignLen + 1;
					end = mid + alignLen - 1;
				}
				if (records.get(i).getAlignmentEnd() >= start) {
					Alignment align = new Alignment(null);
					align.addAlignment(records.get(i));
					addSeg(alignTree, ip_flag && replicate && i == 0 ? treeTR1 : ip_flag && replicate && i == 1 ? treeTR2 : null, mutTable.getOrDefault(chr, emptyMutTree),
							align, args.isRtMut() ? 0 : getClipStat(records.get(i)));
				}
				records.set(i, its.get(i).hasNext() ? its.get(i).next() : null);
			}
		}
	}
	
	private void removeRegion(IntervalTree<List<String>> alignTree, int nextStart) {
		for (Iterator<Node<List<String>>> alignNodes = alignTree.overlappers(0, nextStart - 1); alignNodes.hasNext();) {
			Node<List<String>> alignNode = alignNodes.next();
			if (alignNode.getEnd() < nextStart) {
				alignNode.getValue().clear();
				alignNode.setValue(null);
				alignTree.remove(alignNode.getStart(), alignNode.getEnd());
			}
		}
	}

	private boolean removeSNP(IntervalTree<Mutation> snpTree, Node<Mutation> snpNode) {
		if (!snpNode.getValue().isValidSNP(args.getControlFiles().isEmpty(), 2)) {
			snpTree.remove(snpNode.getStart(), snpNode.getEnd());
			return true;
		}
		for (Iterator<Node<Mutation>> snpWinNodes = snpTree.overlappers(snpNode.getStart() - window + 1, snpNode.getEnd() + window - 1);
				snpWinNodes.hasNext();) {
			Node<Mutation> snpWinNode = snpWinNodes.next();
			if (snpWinNode == snpNode) {
				continue;
			}
			if (snpWinNode.getValue().isValidSNP(args.getControlFiles().isEmpty(), 2)) {
				snpWinNode.getValue().setNeedRemove();
				snpNode.getValue().setNeedRemove();
			}
		}
		if (snpNode.getValue().needRemove()) {
			snpTree.remove(snpNode.getStart(), snpNode.getEnd());
			return true;
		}
		return false;
	}

//	private double[] calScores(int[][] counts, int start, int end) {
//		IntervalTree<double[]> scoreTree = new IntervalTree<>();
//		int winStart = start;
//		double[] scores = new double[end - start + 1];
//		int scoreStart = 0;
//		for (int winEnd = start + window - 1; winEnd < end; winEnd += winStep) {
//			scoreTree.put(winStart, winEnd, normalize2(getRTstopScores(counts, winStart - start, winEnd - start + 1)));
//			int nextStart = Math.min(winEnd + winStep, end) - window + 1;
//			double[] winScore = calScores(scoreTree, winStart, nextStart);
//			for (int i = 0; i < winScore.length; ++i, ++scoreStart) {
//				scores[scoreStart] = winScore[i];
//			}
//			winStart = nextStart;
//		}
//		scoreTree.put(winStart, end, normalize2(getRTstopScores(counts, winStart - start, end - start + 1)));
//		double[] winScore = calScores(scoreTree, winStart, end + 1);
//		for (int i = 0; i < winScore.length; ++i, ++scoreStart) {
//			scores[scoreStart] = winScore[i];
//		}
//		return scores;
//	}
	
	private double[] calScores(IntervalTree<double[]> scoreTree, int start, int nextStart) {
		double[] scores = new double[nextStart - start];
		int[] counts = new int[nextStart - start];
		for (Iterator<Node<double[]>> scoreNodes = scoreTree.overlappers(0, nextStart - 1); scoreNodes.hasNext();) {
			Node<double[]> scoreNode =  scoreNodes.next();
			if (scoreNode.getEnd() >= start) {
				int end = Math.min(scoreNode.getEnd() + 1, nextStart);
				for (int i = start; i < end; ++i) {
					double score = scoreNode.getValue()[i - scoreNode.getStart()];
					if (score < 0.0) {
						continue;
					}
					scores[i - start] += score;
					++counts[i - start];
				}
			}
			if (scoreNode.getEnd() < nextStart) {
				scoreNode.setValue(null);
				scoreTree.remove(scoreNode.getStart(), scoreNode.getEnd());
			}
		}
		for (int i = 0; i < scores.length; ++i) {
			scores[i] = counts[i] == 0 ? -1.0 : scores[i] / counts[i];
		}
		return scores;
	}
	
	public double[] getWinShapeScores(int[][] counts, int start, int end, String seq) {
		double[] scores = new double[end - start + 1];
		IntervalTree<double[]> scoreTree = new IntervalTree<>();
		int winStart = start;
			
		int scoreStart = 0;
		for (int winEnd = start + window - 1; winEnd < end; winEnd += winStep) {
			double[] winScore = normalize2(getShapeScores(counts, winStart - start, winEnd - start + 1,
					seq == null ? null : seq.substring(winStart - start, winEnd - start + 1)));
			if (winScore != null) {
				winScore = filtScoreFill(winScore, seq == null ? null : seq.substring(winStart - start, winEnd - start + 1));
				scoreTree.put(winStart, winEnd, winScore);
			}
			int nextStart = Math.min(winEnd + winStep, end) - window + 1;
			winScore = calScores(scoreTree, winStart, nextStart);
			for (int i = 0; i < winScore.length; ++i, ++scoreStart) {
				scores[scoreStart] = winScore[i];
			}
			winStart = nextStart;
		}
		double[] winScore = normalize2(getShapeScores(counts, winStart - start, end - start + 1, 
				seq == null ? null : seq.substring(winStart - start, end - start + 1)));
		if (winScore != null) {
			winScore = filtScoreFill(winScore, seq == null ? null : seq.substring(winStart - start, end - start + 1));
			scoreTree.put(winStart, end, winScore);
		}
		winScore = calScores(scoreTree, winStart, end + 1);
		for (int i = 0; i < winScore.length; ++i, ++scoreStart) {
			scores[scoreStart] = winScore[i];
		}
		return scores;
	}
	
//	private void getSNPinSortBam(String treat_file, Genome genome,
//			HashMap<String, IntervalTree<Mutation>> mut_table) throws IOException {
//		try(SamReader readerT = SamReaderFactory.makeDefault().open(new File(treat_file))) {
//			if (!Method.isCoordinate(readerT.getFileHeader())) {
//				return;
//			}
//			List<Chromosome> chrs = Method.getChrs(readerT.getFileHeader());
//			HashMap<String, Integer> chrMap = new HashMap<>();
//			Method.appendFile(InParam.getParams().getOutPrefix() + "_mut.txt", chrMap.entrySet(), Mutation.getSimpleHeader());
//			for (int i = 0; i < chrs.size(); ++i) {
//				chrMap.put(chrs.get(i).getChr(), i);
//			}
//			Iterator<SAMRecord> itT = readerT.iterator();
//			SAMRecord recordT = itT.hasNext() ? itT.next() : null;
//			for (int chrIndex = 0; chrIndex < chrs.size(); ++chrIndex) {
//				String chr = chrs.get(chrIndex).getChr();
//				Method.printNow("Dealing with " + chr);
//				if (!mut_table.containsKey(chr)) {
//					continue;
//				}
//				IntervalTree<List<String>> treeT = new IntervalTree<>();
//				
//				while (itT.hasNext() && chrIndex > chrMap.getOrDefault(recordT.getReferenceName(), -1)){
//					recordT = itT.next();
//				}
//				ip_flag = true;
//				while (itT.hasNext() && sameChr(recordT, chr)) {
//					Alignment align = new Alignment(null);
//					align.addMut(mut_table.get(chr), recordT);
//					addSeg(treeT, mut_table.get(chr), align);
//					window = Math.max(window, recordT.getReadLength());
//					recordT = itT.next();
//				}
//				for (Node<Mutation> mutNode : mut_table.get(chr)) {
//					Mutation mutation = mutNode.getValue();
//					int windowSize = InParam.getParams().getWindowSize() > 0 ? InParam.getParams().getWindowSize() : this.window;
//					mutation.remakeMostAlt();
//					if (!matchSNPCounts(mutNode, InParam.getParams().getExam())) {
//						mut_table.get(chr).remove(mutNode.getStart(), mutNode.getEnd());
//						continue;
//					}
//					boolean up = mutation.getScript().getGene().getStrand() == '+';
//					mutation.setScores(getRTstopScores(treeT, null, null, genome.getChr(chr).getSeq(), mutNode.getStart() - windowSize + 1, mutNode.getEnd() + windowSize - 1, up));
////					mutation.setAltScores(getRTstopSNPScores(treeT, mutNode, windowSize, true, up));
////					mutation.setRefScores(getRTstopSNPScores(treeT, mutNode, windowSize, false, up));
////					mutation.setAltTestScores(testRTstopSNPScores(treeT, mutNode, windowSize, 1000, mutation.getAltScores(), up));
////					mutation.setRefTestScores(testRTstopSNPScores(treeT, mutNode, windowSize, 1000, mutation.getRefScores(), up));
////					mutation.setAltNoneScores(getRTstopNoSNPScores(treeT, mutNode, windowSize, true, up));
////					mutation.setRefNoneScores(getRTstopNoSNPScores(treeT, mutNode, windowSize, false, up));
//				}
//				HashMap<String, IntervalTree<Mutation>> out_table = new HashMap<>();
//				out_table.put(chr, mut_table.get(chr));
//				Method.appendFile(InParam.getParams().getOutPrefix() + "_mut.txt", Method.getInstance().mutOutput(genome, out_table), null);
//				mut_table.remove(chr);
//			}
//		}
//	}
//	
//	
//	
//	private void getSNPinSortBam(String treat_file, String control_file, Genome genome,
//			HashMap<String, IntervalTree<Mutation>> mut_table) throws IOException {		
//		try(SamReader readerT = SamReaderFactory.makeDefault().open(new File(treat_file));
//			SamReader readerC = SamReaderFactory.makeDefault().open(new File(control_file))) {
//			if (!Method.isCoordinate(readerT.getFileHeader()) || !Method.isCoordinate(readerC.getFileHeader())) {
//				return;
//			}
//			List<Chromosome> chrs = Method.getChrs(readerT.getFileHeader());
//			HashMap<String, Integer> chrMap = new HashMap<>();
//			Method.appendFile(InParam.getParams().getOutPrefix() + "_mut.txt", chrMap.entrySet(), Mutation.getSimpleHeader());
//			for (int i = 0; i < chrs.size(); ++i) {
//				chrMap.put(chrs.get(i).getChr(), i);
//			}
//			Iterator<SAMRecord> itT = readerT.iterator();
//			Iterator<SAMRecord> itC = readerC.iterator();
//			SAMRecord recordT = itT.hasNext() ? itT.next() : null;
//			SAMRecord recordC = itC.hasNext() ? itC.next() : null;
//			for (int chrIndex = 0; chrIndex < chrs.size(); ++chrIndex) {
//				String chr = chrs.get(chrIndex).getChr();
//				Method.printNow("Dealing with " + chr);
//				if (!mut_table.containsKey(chr)) {
//					continue;
//				}
//				IntervalTree<List<String>> treeT = new IntervalTree<>();
//				IntervalTree<List<String>> treeC = new IntervalTree<>();
//				
//				while (itT.hasNext() && chrIndex > chrMap.getOrDefault(recordT.getReferenceName(), -1)){
//					recordT = itT.next();
//				}
//				while (itC.hasNext() && chrIndex > chrMap.getOrDefault(recordC.getReferenceName(), -1)) {
//					recordC = itC.next();
//				}
//				ip_flag = true;
//				while (itT.hasNext() && sameChr(recordT, chr)) {
//					Alignment align = new Alignment(null);
//					align.addMut(mut_table.get(chr), recordT);
//					addSeg(treeT, mut_table.get(chr), align);
//					window = Math.max(window, recordT.getReadLength());
//					recordT = itT.next();
//				}
//				ip_flag = false;
//				while (itC.hasNext() && sameChr(recordC, chr)) {
//					Alignment align = new Alignment(null);
//					align.addMut(mut_table.get(chr), recordC);
//					addSeg(treeC, mut_table.get(chr), align);
//					window = Math.max(window, recordC.getReadLength());
//					recordC = itC.next();
//				}
//				for (Node<Mutation> mutNode : mut_table.get(chr)) {
//					Mutation mutation = mutNode.getValue();
//					int windowSize = InParam.getParams().getWindowSize() > 0 ? InParam.getParams().getWindowSize() : this.window;
//					mutation.remakeMostAlt();
//					if (!matchSNPCounts(mutNode, InParam.getParams().getExam())) {
//						mut_table.get(chr).remove(mutNode.getStart(), mutNode.getEnd());
//						continue;
//					}
//					boolean up = mutation.getScript().getGene().getStrand() == '+';
//					mutation.setScores(getRTstopScores(treeT, treeC, null, genome.getChr(chr).getSeq(), mutNode.getStart() - windowSize + 1, mutNode.getEnd() + windowSize - 1, up));
//					mutation.setAltScores(getRTstopSNPScores(treeT, treeC, mutNode, windowSize, true, up));
//					mutation.setRefScores(getRTstopSNPScores(treeT, treeC, mutNode, windowSize, false, up));
//					mutation.setAltTestScores(testRTstopSNPScores(treeT, treeC, mutNode, windowSize, 1000, mutation.getAltScores(), up));
//					mutation.setRefTestScores(testRTstopSNPScores(treeT, treeC, mutNode, windowSize, 1000, mutation.getRefScores(), up));
//					mutation.setAltNoneScores(getRTstopNoSNPScores(treeT, treeC, mutNode, windowSize, true, up));
//					mutation.setRefNoneScores(getRTstopNoSNPScores(treeT, treeC, mutNode, windowSize, false, up));
//				}
//				HashMap<String, IntervalTree<Mutation>> out_table = new HashMap<>();
//				out_table.put(chr, mut_table.get(chr));
//				Method.appendFile(InParam.getParams().getOutPrefix() + "_mut.txt", Method.getInstance().mutOutput(genome, out_table), null);
//				mut_table.remove(chr);
//			}
//		}
//	}
//	
//	private void getSNPinSortBam(List<String> treatFiles, Genome genome,
//			HashMap<String, IntervalTree<Mutation>> mut_table) throws IOException {
//		List<Iterator<SAMRecord>> itsT = new ArrayList<>();
//		List<Chromosome> chrs = null;
//		try {
//			for (int i = 0; i < treatFiles.size(); ++i) {
//				SamReader readerT = SamReaderFactory.makeDefault().open(new File(treatFiles.get(i)));
//				if (!Method.isCoordinate(readerT.getFileHeader())) {
//					return;
//				}
//				if (i == 0) {
//					chrs = Method.getChrs(readerT.getFileHeader());
//				}
//				itsT.add(readerT.iterator());
//			}
//			if (chrs == null) {
//				return;
//			}
//			HashMap<String, Integer> chrMap = new HashMap<>();
//			Method.appendFile(InParam.getParams().getOutPrefix() + "_mut.txt", chrMap.entrySet(), Mutation.getSimpleHeader());
//			for (int i = 0; i < chrs.size(); ++i) {
//				chrMap.put(chrs.get(i).getChr(), i);
//			}
//			
//			List<SAMRecord> recordsT = new ArrayList<>(itsT.size());
//			for (int i = 0; i < itsT.size(); ++i) {
//				SAMRecord recordT = itsT.get(i).hasNext() ? itsT.get(i).next() : null;
//				recordsT.add(recordT);
//			}
//			for (int chrIndex = 0; chrIndex < chrs.size(); ++chrIndex) {
//				String chr = chrs.get(chrIndex).getChr();
//				Method.printNow("Dealing with " + chr);
//				if (!mut_table.containsKey(chr)) {
//					continue;
//				}
//				IntervalTree<List<String>> treeT = new IntervalTree<>();
//				
//				ip_flag = true;
//				for (int i = 0; i < itsT.size(); ++i) {
//					while (itsT.get(i).hasNext() && chrIndex > chrMap.getOrDefault(recordsT.get(i).getReferenceName(), -1)){
//						recordsT.set(i, itsT.get(i).next());
//					}
//					while (itsT.get(i).hasNext() && sameChr(recordsT.get(i), chr)) {
//						Alignment align = new Alignment(null);
//						align.addMut(mut_table.get(chr), recordsT.get(i));
//						addSeg(treeT, mut_table.get(chr), align);
//						window = Math.max(window, recordsT.get(i).getReadLength());
//						recordsT.set(i, itsT.get(i).next());
//					}
//				}
//				
//				for (Node<Mutation> mutNode : mut_table.get(chr)) {
//					Mutation mutation = mutNode.getValue();
//					int windowSize = InParam.getParams().getWindowSize() > 0 ? InParam.getParams().getWindowSize() : this.window;
//					mutation.remakeMostAlt();
//					if (!matchSNPCounts(mutNode, InParam.getParams().getExam())) {
//						mut_table.get(chr).remove(mutNode.getStart(), mutNode.getEnd());
//						continue;
//					}
//					boolean up = mutation.getScript().getGene().getStrand() == '+';
//					mutation.setScores(getRTstopScores(treeT, null, null, genome.getChr(chr).getSeq(), mutNode.getStart() - windowSize + 1, mutNode.getEnd() + windowSize - 1, up));
//					mutation.setAltScores(getRTstopSNPScores(treeT, null, mutNode, windowSize, true, up));
//					mutation.setRefScores(getRTstopSNPScores(treeT, null, mutNode, windowSize, false, up));
//					mutation.setAltTestScores(testRTstopSNPScores(treeT, null, mutNode, windowSize, 1000, mutation.getAltScores(), up));
//					mutation.setRefTestScores(testRTstopSNPScores(treeT, null, mutNode, windowSize, 1000, mutation.getRefScores(), up));
//					mutation.setAltNoneScores(getRTstopNoSNPScores(treeT, null, mutNode, windowSize, true, up));
//					mutation.setRefNoneScores(getRTstopNoSNPScores(treeT, null, mutNode, windowSize, false, up));
//				}
//				HashMap<String, IntervalTree<Mutation>> out_table = new HashMap<>();
//				out_table.put(chr, mut_table.get(chr));
//				Method.appendFile(InParam.getParams().getOutPrefix() + "_mut.txt", Method.getInstance().mutOutput(genome, out_table), null);
//				mut_table.remove(chr);
//			}
//		} finally {
//		}
//	}
//	
	private void getSamplesinSortBam(Genome genome) {
		List<SamReader> readersT = new ArrayList<>(args.getTreatFiles().size());
		List<SamReader> readersC = new ArrayList<>(args.getControlFiles().size());
		List<SamReader> readersT2 = new ArrayList<>(args.getTreatFiles2().size());
		List<SamReader> readersC2 = new ArrayList<>(args.getControlFiles2().size());
		try {
			ip_flag = true;
			readersT = Method.openBamFiles(args.getTreatFiles());
			readersT2 = Method.openBamFiles(args.getTreatFiles2());
			ip_flag = false;
			readersC = Method.openBamFiles(args.getControlFiles());
			readersC2 = Method.openBamFiles(args.getControlFiles2());
			if (readersT.size() < args.getTreatFiles().size() || readersC.size() < args.getControlFiles().size() ||
					readersT2.size() < args.getTreatFiles2().size() || readersC2.size() < args.getControlFiles2().size()) {
				return;
			}
			chrs = Method.getChrs(readersT.get(0).getFileHeader());
			itsT = new ArrayList<>(args.getTreatFiles().size());
			itsC = new ArrayList<>(args.getControlFiles().size());
			List<Iterator<SAMRecord>> itsT2 = new ArrayList<>(args.getTreatFiles2().size());
			List<Iterator<SAMRecord>> itsC2 = new ArrayList<>(args.getControlFiles2().size());
			for (SamReader reader : readersT) {
				itsT.add(reader.iterator());
			}
			for (SamReader reader : readersC) {
				itsC.add(reader.iterator());
			}
			for (SamReader reader : readersT2) {
				itsT2.add(reader.iterator());
			}
			for (SamReader reader : readersC2) {
				itsC2.add(reader.iterator());
			}
			
			String outFile = args.getPeakFile() == null ? args.getOutPrefix() + "_gene.txt" : args.getOutPrefix() + "_region.txt";
			Method.writeFile(outFile, new ArrayList<>(), Mutation.getSimpleHeader());
			Map<String, IntervalTree<Peak>> regionTable = args.getPeakFile() == null ? null : 
				Method.getInstance().readPeak(args.getPeakFile(), null);
			
			recordsT = new ArrayList<>(itsT.size());
			recordsC = new ArrayList<>(itsC.size());
			List<SAMRecord> recordsT2 = new ArrayList<>(itsT2.size());
			List<SAMRecord> recordsC2 = new ArrayList<>(itsC2.size());
			for (int i = 0; i < itsT.size(); ++i) {
				SAMRecord recordT = itsT.get(i).hasNext() ? itsT.get(i).next() : null;
				recordsT.add(recordT);
			}
			for (int i = 0; i < itsC.size(); ++i) {
				SAMRecord recordC = itsC.get(i).hasNext() ? itsC.get(i).next() : null;
				recordsC.add(recordC);
			}
			for (int i = 0; i < itsT2.size(); ++i) {
				SAMRecord recordT2 = itsT2.get(i).hasNext() ? itsT2.get(i).next() : null;
				recordsT2.add(recordT2);
			}
			for (int i = 0; i < itsC2.size(); ++i) {
				SAMRecord recordC2 = itsC2.get(i).hasNext() ? itsC2.get(i).next() : null;
				recordsC2.add(recordC2);
			}
			
			for (int chrIndex = 0; chrIndex < chrs.size(); ++chrIndex) {
				String chr = chrs.get(chrIndex);
				if (genome.getChr(chr) == null) {
					continue;
				}
				treeT = new IntervalTree<>();
				treeC = new IntervalTree<>();
				treeTH = new IntervalTree<>();
				treeCH = new IntervalTree<>();
				List<String> outList = new ArrayList<>();
				IntervalTree<Gene> geneTree = genome.getChr(chr).getGeneTree();
				if (regionTable != null) {
					geneTree = new IntervalTree<>();
					for (Node<Peak> regionNode : regionTable.getOrDefault(chr, new IntervalTree<>())) {
						Gene gene = new Gene(null, null, null, null, '.', regionNode.getStart() - 1, regionNode.getEnd());
						Node<Gene> geneNode = genome.getChr(chr).getGeneTree() != null ? null : 
							genome.getChr(chr).getGeneTree().minOverlapper(regionNode.getStart(), regionNode.getEnd());
						gene.setStrand(geneNode == null ? '.' : geneNode.getValue().getStrand());
						geneTree.put(regionNode.getStart(), regionNode.getEnd(), gene);
					}
				}
				Method.printNow("Search " + chr + " with " + geneTree.size() + " regions");
				for (Iterator<Node<Gene>> geneIt = geneTree.overlappers(0, Integer.MAX_VALUE);
						geneIt.hasNext();) {
					Node<Gene> geneNode = geneIt.next();
					int start = geneNode.getStart();
					int end = geneNode.getEnd();
					boolean positive = geneNode.getValue().getStrand() != '-';
					removeRegion(treeT, start);
					removeRegion(treeC, start);
					removeRegion(treeTH, start);
					removeRegion(treeCH, start);
					
					ip_flag = true;
					addRegion(treeT, itsT, recordsT, chr, start, end);
					addRegion(treeTH, itsT2, recordsT2, chr, start, end);
					ip_flag = false;
					addRegion(treeC, itsC, recordsC, chr, start, end);
					addRegion(treeCH, itsC2, recordsC2, chr, start, end);
					
					int[][][] counts = getCounts(null, treeT, treeC, null, getChrSeq(genome, chr), start, end, positive);
					int[][][] counts2 = getCounts(null, treeTH, treeCH, null, getChrSeq(genome, chr), start, end, positive);
					
 					outList.add(getOutGeneScores(counts, counts2, chr, genome, start, end, positive));
				}
				Method.appendFile(outFile, outList, null);
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				for (SamReader reader : readersT) {
					reader.close();
				}
				for (SamReader reader : readersC) {
					reader.close();
				}
			} catch (IOException e) {
			}
		}
	}
	
	private void getSampleSnpinSortBam(Genome genome) {
		List<SamReader> readersT = new ArrayList<>(args.getTreatFiles().size());
		List<SamReader> readersC = new ArrayList<>(args.getControlFiles().size());
		List<SamReader> readersT2 = new ArrayList<>(args.getTreatFiles2().size());
		List<SamReader> readersC2 = new ArrayList<>(args.getControlFiles2().size());
		try {
			ip_flag = true;
			readersT = Method.openBamFiles(args.getTreatFiles());
			readersT2 = Method.openBamFiles(args.getTreatFiles2());
			ip_flag = false;
			readersC = Method.openBamFiles(args.getControlFiles());
			readersC2 = Method.openBamFiles(args.getControlFiles2());
			if (readersT.size() < args.getTreatFiles().size() || readersC.size() < args.getControlFiles().size() ||
					readersT2.size() < args.getTreatFiles2().size() || readersC2.size() < args.getControlFiles2().size()) {
				return;
			}
			chrs = Method.getChrs(readersT.get(0).getFileHeader());
			replicate = args.isReplicate() && (readersT.size() > 1 || readersT2.size() > 1);
			itsT = new ArrayList<>(args.getTreatFiles().size());
			itsC = new ArrayList<>(args.getControlFiles().size());
			List<Iterator<SAMRecord>> itsT2 = new ArrayList<>(args.getTreatFiles2().size());
			List<Iterator<SAMRecord>> itsC2 = new ArrayList<>(args.getControlFiles2().size());
			for (SamReader reader : readersT) {
				itsT.add(reader.iterator());
			}
			for (SamReader reader : readersC) {
				itsC.add(reader.iterator());
			}
			for (SamReader reader : readersT2) {
				itsT2.add(reader.iterator());
			}
			for (SamReader reader : readersC2) {
				itsC2.add(reader.iterator());
			}
			
			String outFile = args.getOutPrefix() + "_riboSNitch.txt";
			String detailFile = args.getOutPrefix() + "_riboSNitchDetail.txt";
			Method.writeFile(outFile, new ArrayList<>(), 
					"#Chr\tPosition\tHeter_Ref\tHeter_Alt\tHomo\tName\tScore\tStrand\tGene\tTranscript\tAnnotation\tRegion_in_Genome|Transcriptome");
			Method.writeFile(detailFile, new ArrayList<>(), 
					"#Name\teSDC\teSDC_pValue\tHeter_diffBase_ratio\tHeter_Replicate_eSDC\tHeter_Counts(A|T|C|G)\tHomo_Counts(A|T|C|G)\tRef_Seq\tHeter_Score|Homo_Score\tStructDiff_pValue");
			
			recordsT = new ArrayList<>(itsT.size());
			recordsC = new ArrayList<>(itsC.size());
			List<SAMRecord> recordsT2 = new ArrayList<>(itsT2.size());
			List<SAMRecord> recordsC2 = new ArrayList<>(itsC2.size());
			
			for (int i = 0; i < itsT.size(); ++i) {
				SAMRecord recordT = itsT.get(i).hasNext() ? itsT.get(i).next() : null;
				recordsT.add(recordT);
			}
			for (int i = 0; i < itsC.size(); ++i) {
				SAMRecord recordC = itsC.get(i).hasNext() ? itsC.get(i).next() : null;
				recordsC.add(recordC);
			}
			for (int i = 0; i < itsT2.size(); ++i) {
				SAMRecord recordT2 = itsT2.get(i).hasNext() ? itsT2.get(i).next() : null;
				recordsT2.add(recordT2);
			}
			for (int i = 0; i < itsC2.size(); ++i) {
				SAMRecord recordC2 = itsC2.get(i).hasNext() ? itsC2.get(i).next() : null;
				recordsC2.add(recordC2);
			}
			
			String snpLine = "";
			BufferedReader snpReader = new BufferedReader(new FileReader(args.getMutFile()));
			
			for (int chrIndex = 0; chrIndex < chrs.size(); ++chrIndex) {
				String chr = chrs.get(chrIndex);
				snpLine = Method.readMutation(snpReader, mutTable, args.getMutFile().endsWith("vcf"), chr, snpLine);
				Method.setStrictScript(mutTable, genome, chr);
				treeT = new IntervalTree<>();
				treeC = new IntervalTree<>();
				treeTH = new IntervalTree<>();
				treeCH = new IntervalTree<>();
				treeTR1 = replicate ? new IntervalTree<>() : null;
				treeTR2 = replicate ? new IntervalTree<>() : null;
				List<String> outList = new ArrayList<>();
				List<String> outDetailList = new ArrayList<>();
				IntervalTree<Mutation> snpTree = mutTable.getOrDefault(chr, emptyMutTree);
				Method.printNow(String.format("Search %s with %d SNP sites", chr, snpTree.size()));
				int validSNP = 0;
				for (Iterator<Node<Mutation>> snpIt = snpTree.overlappers(0, Integer.MAX_VALUE);
						snpIt.hasNext();) {
					Node<Mutation> snpNode = snpIt.next();
					int start = snpNode.getStart() - alignLen + 1;
					int end = snpNode.getEnd() + alignLen - 1 ;
					boolean positive = snpNode.getValue().getScript().getGene().getStrand() != '-';
					removeRegion(treeT, start);
					removeRegion(treeC, start);
					removeRegion(treeTH, start);
					removeRegion(treeCH, start);
					
					ip_flag = true;
					addRegion(treeT, itsT, recordsT, chr, start, end);
					if (!itsC.isEmpty()) {
						ip_flag = false;
						addRegion(treeC, itsC, recordsC, chr, start, end);
					}
					if (removeSNP(snpTree, snpNode)) {
//						StringBuilder sb = new StringBuilder();
//						sb.append(chr);
//						sb.append('\t');
//						sb.append(snpNode.getStart() - 1);
//						sb.append('\t');
//						sb.append(snpNode.getEnd());
//						sb.append('\t');
//						sb.append(snpNode.getValue().toSimpleString(-3.0, -3.0));
//						outList.add(sb.toString());
						continue;
					}
					snpNode.getValue().remakeMostAlt(args.getControlFiles().isEmpty());
					if (snpNode.getValue().getSNP().length != 1) {
//						StringBuilder sb = new StringBuilder();
//						sb.append(chr);
//						sb.append('\t');
//						sb.append(snpNode.getStart() - 1);
//						sb.append('\t');
//						sb.append(snpNode.getEnd());
//						sb.append('\t');
//						sb.append(snpNode.getValue().toSimpleString(-4.0, -4.0));
//						outList.add(sb.toString());
						continue;
					}
					
					ip_flag = true;
					addRegion(treeTH, itsT2, recordsT2, chr, start, end);
					if (!itsC2.isEmpty()) {
						ip_flag = false;
						addRegion(treeCH, itsC2, recordsC2, chr, start, end);
					}
					
					char snp = snpNode.getValue().getSNP()[0];
					int[] region = changeSampleSnpRegion(snpNode.getStart(), snpNode.getEnd(), snpNode, snp);
					if (region[0] <= region[1]) {
						start = region[0];
						end = region[1];
					}
					
					int[][][] counts = getCounts(null, treeT, treeC, snpNode, getChrSeq(genome, chr), start, end, positive);
					int[][][] counts2 = getCounts(null, treeTH, treeCH, snpNode, getChrSeq(genome, chr), start, end, positive);
					
					int[] baseCounts = new int[counts2.length];
					int[] baseCounts2 = new int[counts.length];
					int baseSum = 0;
					int maxIndex = 0;
					for (int i = 0; i < counts2.length; ++i) {
						baseCounts[i] = counts2[i].length < 4 ? counts2[i][1][snpNode.getStart() - start] : counts2[i][3][snpNode.getStart() - start];
						baseCounts2[i] = counts[i].length < 4 ? counts[i][1][snpNode.getStart() - start] : counts[i][3][snpNode.getStart() - start] ;
						baseSum += baseCounts[i];
						maxIndex = baseCounts[maxIndex] < baseCounts[i] ? i : maxIndex;
					}
					
//					int[] debugCounts = baseCounts.clone();
					Arrays.sort(baseCounts);
					if ((maxIndex != getBaseIndex(snp) && maxIndex != getBaseIndex(snpNode.getValue().getRef()))
							|| baseCounts[baseCounts.length - 1] * 2 <= baseSum || baseCounts[baseCounts.length - 2] * 10 >= 3 *baseSum) {
//						StringBuilder sb = new StringBuilder();
//						sb.append(chr);
//						sb.append('\t');
//						sb.append(snpNode.getStart() - 1);
//						sb.append('\t');
//						sb.append(snpNode.getEnd());
//						sb.append('\t');
//						sb.append(snpNode.getValue().toSimpleString(-5.0, -5.0));
//						sb.append('\t');
//						sb.append(Arrays.toString(baseCounts2));
//						sb.append('\t');
//						sb.append(Arrays.toString(debugCounts));
//						outList.add(sb.toString());
						continue;
					}
					
					String[] outLines = getOutSnpScores(counts, counts2, chr, getChrSeq(genome, chr), snpNode, snp,
							getChar(maxIndex), start, end, positive);
					outList.add(outLines[0]);
					outDetailList.add(outLines[1]);
					++validSNP;
//					outList.add(getOutGeneScores(counts, counts2, geneNode.getValue(), getChrSeq(genome, chr), start, end, positive));
//				}
//				Method.appendFile(InParam.getParams().getOutPrefix() + "_mut.txt", outList, null);
				}
				Method.appendFile(outFile, outList, null);
				Method.appendFile(detailFile, outDetailList, null);
				Method.printNow(String.format("Found %d existed SNPs in %s", validSNP, chr));
				genome.getChrMap().remove(chr);
				mutTable.remove(chr);
			}
			snpReader.close();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				for (SamReader reader : readersT) {
					reader.close();
				}
				for (SamReader reader : readersC) {
					reader.close();
				}
			} catch (IOException e) {
			}
		}
	}
	
	private void getSNPinSortBam(List<String> treatFiles, List<String> controlFiles, Genome genome) {
		List<SamReader> readersT = new ArrayList<>(treatFiles.size());
		List<SamReader> readersC = new ArrayList<>(controlFiles.size());
		try {
			ip_flag = true;
			readersT = Method.openBamFiles(treatFiles);
			ip_flag = false;
			readersC = Method.openBamFiles(controlFiles);
			if (readersT.size() < treatFiles.size() || readersC.size() < controlFiles.size()) {
				return;
			}
			chrs = Method.getChrs(readersT.get(0).getFileHeader());
			itsT = new ArrayList<>(treatFiles.size());
			itsC = new ArrayList<>(controlFiles.size());
			for (SamReader reader : readersT) {
				itsT.add(reader.iterator());
			}
			for (SamReader reader : readersC) {
				itsC.add(reader.iterator());
			}
			
			String outFile = args.getOutPrefix() + "_riboSNitch.txt";
			String detailFile = args.getOutPrefix() + "_riboSNitchDetail.txt";
			Method.writeFile(outFile, new ArrayList<>(), 
					"#Chr\tPosition\tRef\tAlt\tName\tScore\tStrand\tGene\tTranscript\tAnnotation\tRegion_in_Genome|Transcriptome");
			Method.writeFile(detailFile, new ArrayList<>(), 
					"#Name\teSDC\teSDC_pValue\tSNP_Reads_Ratio\tReplicate_eSDC\tTreat_Counts(A|T|C|G)\tControl_Counts(A|T|C|G)\tRef_Seq\tRef_Score|Alt_Score\tStructDiff_pValue");
			
			
			recordsT = new ArrayList<>(itsT.size());
			recordsC = new ArrayList<>(itsC.size());
			for (int i = 0; i < itsT.size(); ++i) {
				SAMRecord recordT = itsT.get(i).hasNext() ? itsT.get(i).next() : null;
				recordsT.add(recordT);
			}
			for (int i = 0; i < itsC.size(); ++i) {
				SAMRecord recordC = itsC.get(i).hasNext() ? itsC.get(i).next() : null;
				recordsC.add(recordC);
			}
			replicate = args.isReplicate() && treatFiles.size() > 1;
			String snpLine = "";
//			String exonLine = "";
//			String genomeLine = null;
			BufferedReader snpReader = new BufferedReader(new FileReader(args.getMutFile()));
//			BufferedReader exonReader = new BufferedReader(new FileReader(args.getExonFile()));
//			BufferedReader genomeReader = new BufferedReader(new FileReader(args.getGenomeFile()));
			for (int chrIndex = 0; chrIndex < chrs.size(); ++chrIndex) {
				String chr = chrs.get(chrIndex);
				snpLine = Method.readMutation(snpReader, mutTable, args.getMutFile().endsWith("vcf"), chr, snpLine);
//				exonLine = genome.readAnnoteFile(exonReader, chr, exonLine);
//				genomeLine = genome.readSeqFile(genomeReader, chr, genomeLine);
				Method.setStrictScript(mutTable, genome, chr);
				treeT = new IntervalTree<>();
				treeC = new IntervalTree<>();
				treeTR1 = replicate ? new IntervalTree<>() : null;
				treeTR2 = replicate ? new IntervalTree<>() : null;
				List<String> outList = new ArrayList<>();
				List<String> outDetailList = new ArrayList<>();
				List<String> outOddList = args.isSharedPeak() ? new ArrayList<>() : null;
				IntervalTree<Mutation> snpTree = mutTable.getOrDefault(chr, emptyMutTree);
				Method.printNow(String.format("Search %s with %d SNP sites", chr, snpTree.size()));
				int validSNP = 0;
				for (Iterator<Node<Mutation>> snpIt = snpTree.overlappers(0, Integer.MAX_VALUE);
						snpIt.hasNext();) {
					Node<Mutation> snpNode = snpIt.next();
					int start = snpNode.getStart() - alignLen + 1;
					int end = snpNode.getEnd() + alignLen - 1 ;
					boolean positive = snpNode.getValue().getScript().getGene().getStrand() != '-';
					removeRegion(treeT, start);
					removeRegion(treeC, start);
					if (replicate) {
						removeRegion(treeTR1, start);
						removeRegion(treeTR2, start);
					}
					
					ip_flag = true;
					addRegion(treeT, itsT, recordsT, chr, start, end);
					if (!itsC.isEmpty()) {
						ip_flag = false;
						addRegion(treeC, itsC, recordsC, chr, start, end);
					}
					
					if (outOddList != null) {
						outOddList.add(getOddLine(chr, snpNode));
					}
					if (removeSNP(snpTree, snpNode)) {
//						StringBuilder sb = new StringBuilder();
//						sb.append(chr);
//						sb.append('\t');
//						sb.append(snpNode.getStart() - 1);
//						sb.append('\t');
//						sb.append(snpNode.getEnd());
//						sb.append('\t');
//						sb.append(snpNode.getValue().toSimpleString(-3.0, -3.0));
//						outList.add(sb.toString());
						continue;
					}
					snpNode.getValue().remakeMostAlt(args.getControlFiles().isEmpty());
					if (snpNode.getValue().getSNP().length != 1) {
//						StringBuilder sb = new StringBuilder();
//						sb.append(chr);
//						sb.append('\t');
//						sb.append(snpNode.getStart() - 1);
//						sb.append('\t');
//						sb.append(snpNode.getEnd());
//						sb.append('\t');
//						sb.append(snpNode.getValue().toSimpleString(-4.0, -4.0));
//						outList.add(sb.toString());
						continue;
					}
					
					char snp = snpNode.getValue().getSNP()[0];
					int[] region = changeSnpRegion(snpNode.getStart(), snpNode.getEnd(), snpNode, snp);
					if (region[0] <= region[1]) {
						start = region[0];
						end = region[1];
					}
					if (!args.isRtMut()) {
						if (positive) {
							end = snpNode.getEnd();
						}
						else {
							start = snpNode.getStart();
						}
					}
					
					int[][][] counts = getCounts(null, treeT, treeC, snpNode, getChrSeq(genome, chr), start, end, positive);
					
					String[] outLines = getOutSnpScores(counts, chr, getChrSeq(genome, chr), snpNode, snp, start, end, positive);
					outList.add(outLines[0]);
					outDetailList.add(outLines[1]);
					++validSNP;
					
//					for (int i = 0; i < snpNode.getValue().getSNP().length; ++i) {
//						char snp = snpNode.getValue().getSNP()[i];
//						int[] region = changeSnpRegion(snpNode.getStart(), snpNode.getEnd(), snpNode, snp);
//						if (region[0] <= region[1]) {
//							start = region[0];
//							end = region[1];
//						}
//						if (!args.isAnnote()) {
//							if (positive) {
//								end = snpNode.getEnd();
//							}
//							else {
//								start = snpNode.getStart();
//							}
//						}
// 						outList.add(getOutSnpScores(counts, chr, getChrSeq(genome, chr), snpNode, snp, start, end, positive));
//					}
					

					
//					int[][][][] counts = args.isAnnote() ? getRTmutCounts(null, treeT, treeC, snpNode, getChrSeq(genome, chr), start, end)
//							: getRTstopCounts(null, treeT, treeC, snpNode, getChrSeq(genome, chr), start, end, positive);
//					double[][] repScores = replicate ? new double[2 * counts.length][] : null;
//					if (replicate) {
//						int [][][][] repCounts = args.isAnnote() ? getRTmutCounts(null, treeT1, treeC, snpNode, getChrSeq(genome, chr), start, end)
//								: getRTstopCounts(null, treeT1, treeC, snpNode, getChrSeq(genome, chr), start, end, positive);
//						for (int i = 0; i < repCounts.length; ++i) {
//							repScores[2 * i] = getWinShapeScores(Method.addArray(repCounts[i]), start, end);
//						}
//						repCounts = args.isAnnote() ? getRTmutCounts(null, treeT2, treeC, snpNode, getChrSeq(genome, chr), start, end)
//								: getRTstopCounts(null, treeT2, treeC, snpNode, getChrSeq(genome, chr), start, end, positive);
//						for (int i = 0; i < repCounts.length; ++i) {
//							repScores[2 * i + 1] = getWinShapeScores(Method.addArray(repCounts[i]), start, end);
//						}
//					}
//					
//					for (int i = 0; i < counts.length; ++i) {
//						if (counts[i].length > 2) {
// 							double snpCount = snpNode.getValue().getCount(snpNode.getValue().getRef(), true)
//									+ snpNode.getValue().getCount(snpNode.getValue().getRef(), false) 
//									+ snpNode.getValue().getCount(snpNode.getValue().getSNP()[i], true)
//									+ snpNode.getValue().getCount(snpNode.getValue().getSNP()[i], false);
//							
// 							if (repScores != null) {
// 								outList.add(getOutCalScores(counts[i], chr, getChrSeq(genome, chr), snpNode, start, end, positive,
//									snpCount, repScores[i * 2], repScores[i * 2 + 1]));
// 							}
// 							else {
// 								outList.add(getOutCalScores(counts[i], chr, getChrSeq(genome, chr), snpNode, start, end, positive,
// 									snpCount, null, null));
// 							}
//						}
//					}
					
				}
				if (outOddList != null) {
					Method.appendFile(args.getOutPrefix() + "_odd.txt", outOddList, null);
				}
				Method.appendFile(outFile, outList, null);
				Method.appendFile(detailFile, outDetailList, null);
				Method.printNow(String.format("Found %d existed SNPs in %s", validSNP, chr));
				genome.getChrMap().remove(chr);
				mutTable.remove(chr);
			}
//			for (Node<Mutation> mutNode : mut_table.get(chr)) {
//				Mutation mutation = mutNode.getValue();
//				int windowSize = InParam.getParams().getWindowSize() > 0 ? InParam.getParams().getWindowSize() : this.window;
//				mutation.remakeMostAlt();
//				if (!matchSNPCounts(mutNode, InParam.getParams().getExam())) {
//					mut_table.get(chr).remove(mutNode.getStart(), mutNode.getEnd());
//					continue;
//				}
//				boolean up = mutation.getScript().getGene().getStrand() == '+';
//				mutation.setScores(getRTstopScores(treeT, treeC, null, genome.getChr(chr).getSeq(), mutNode.getStart() - windowSize + 1, mutNode.getEnd() + windowSize - 1, up));
//				mutation.setAltScores(getRTstopSNPScores(treeT, treeC, mutNode, windowSize, true, up));
//				mutation.setRefScores(getRTstopSNPScores(treeT, treeC, mutNode, windowSize, false, up));
//				mutation.setAltTestScores(testRTstopSNPScores(treeT, treeC, mutNode, windowSize, 1000, mutation.getAltScores(), up));
//				mutation.setRefTestScores(testRTstopSNPScores(treeT, treeC, mutNode, windowSize, 1000, mutation.getRefScores(), up));
//				mutation.setAltNoneScores(getRTstopNoSNPScores(treeT, treeC, mutNode, windowSize, true, up));
//				mutation.setRefNoneScores(getRTstopNoSNPScores(treeT, treeC, mutNode, windowSize, false, up));
//			}
//			HashMap<String, IntervalTree<Mutation>> out_table = new HashMap<>();
//			out_table.put(chr, mut_table.get(chr));
//			Method.appendFile(InParam.getParams().getOutPrefix() + "_mut.txt", Method.getInstance().mutOutput(genome, out_table), null);
//			mut_table.remove(chr);
			snpReader.close();
//			exonReader.close();
//			genomeReader.close();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				for (SamReader reader : readersT) {
					reader.close();
				}
				for (SamReader reader : readersC) {
					reader.close();
				}
			} catch (IOException e) {
			}
		}
	}
	
	private boolean sameChr(SAMRecord record, String chr) {
		return record != null && chr.equals(record.getReferenceName());
	}
	
	private int[] changeSampleSnpRegion(int start, int end, Node<Mutation> snpNode, char snp) {
		int[] region = new int[] {Integer.MAX_VALUE, 0};
		if (start > end) {
			return region;
		}
		IntervalTree<List<String>> tree = args.getControlFiles2().size() > 0 ? treeCH : treeTH;
		for (Iterator<Node<List<String>>> segNodes = tree.overlappers(start, end); segNodes.hasNext();) {
			Node<List<String>> segNode = segNodes.next();
			region[0] = Math.min(region[0], segNode.getStart());
			region[1] = Math.max(region[1], segNode.getEnd());
		}
		int[] snpRegion = new int[] {Integer.MAX_VALUE, 0};
		tree = args.getControlFiles().size() > 0 ? treeC : treeT;
		for (Iterator<Node<List<String>>> segNodes = tree.overlappers(start, end); segNodes.hasNext();) {
			Node<List<String>> segNode = segNodes.next();
			snpRegion[0] = Math.min(snpRegion[0], segNode.getStart());
			snpRegion[1] = Math.max(snpRegion[1], segNode.getEnd());
		}
		region[0] = Math.max(region[0], snpRegion[0]);
		region[1] = Math.min(region[1], snpRegion[1]);
		return region;
	}
	
	private int[] changeSnpRegion(int start, int end, Node<Mutation> snpNode, char snp) {
		int[] region = new int[] {Integer.MAX_VALUE, 0};
		if (start > end) {
			return region;
		}
		int[] snpRegion = new int[] {Integer.MAX_VALUE, 0};
		IntervalTree<List<String>> tree = args.getControlFiles().size() > 0 ? treeC : treeT;
		for (Iterator<Node<List<String>>> segNodes = tree.overlappers(start, end); segNodes.hasNext();) {
			Node<List<String>> segNode = segNodes.next();
			for (String seg : segNode.getValue()) {
				int startAlign = seg.charAt(0) == 'S' ? 1 : 0;
				char snpAlign = snpNode == null || snpNode.getStart() > segNode.getEnd() || snpNode.getEnd() < segNode.getStart() ?
						0 : Character.toUpperCase(seg.charAt(snpNode.getStart() - segNode.getStart() + startAlign));
				int snpIndex = snpNode == null ? 0 : getSnpStat(snpAlign, snpNode.getValue().getRef(), snp);
				if (snpIndex == 1) {
					region[0] = Math.min(region[0], segNode.getStart());
					region[1] = Math.max(region[1], segNode.getEnd());
				}
				else if (snpIndex == 2) {
					snpRegion[0] = Math.min(snpRegion[0], segNode.getStart());
					snpRegion[1] = Math.max(snpRegion[1], segNode.getEnd());
				}
			}
		}
		region[0] = Math.max(region[0], snpRegion[0]);
		region[1] = Math.min(region[1], snpRegion[1]);
		return region;
	}
	
//	private int[][][][] getRTstopCounts(int[][][][] counts,IntervalTree<List<String>> treeT, IntervalTree<List<String>> treeC, 
//			Node<Mutation> snpNode, String seq, int start, int end, boolean upstream) {
//		int snpAlt = snpNode == null ? 1 : snpNode.getValue().getSNP().length;
//		int snpNum = snpNode == null ? 1 : 5;
//		int tc = InParam.getParams().getControlFiles().isEmpty() ? 2 : 4;
//		counts = counts != null ? counts : new int[snpAlt][snpNum][tc][end - start + 1];
//		for (int i = 0; i < counts.length; ++i) {
//			getRTstopCounts(counts[i], 0, treeT, snpNode, seq, start, end, upstream);
//			treatReads = controlReads;
//			if (!InParam.getParams().getControlFiles().isEmpty()) {
//				getRTstopCounts(counts[i], 2, treeC, snpNode, seq, start, end, upstream);
//			}
//		}
//		return counts;
//	}
	
	private int[][][] getCounts(int[][][] counts, IntervalTree<List<String>> treeT, IntervalTree<List<String>> treeC, 
			Node<Mutation> snpNode, String seq, int start, int end, boolean upstream) {
		counts = counts == null ? new int[snpNode == null ? 1 : 5][args.getControlFiles().isEmpty() ? 2 : 4][end - start + 1] : counts;
		if (args.isRtMut()) {
			getRTmutCounts(counts, 0, treeT, snpNode, seq, start, end);
			treatReads = controlReads;
			if (!args.getControlFiles().isEmpty()) {
				getRTmutCounts(counts, 2, treeC, snpNode, seq, start, end);
			}
		}
		else {
			getRTstopCounts(counts, 0, treeT, snpNode, seq, start, end, upstream);
			treatReads = controlReads;
			if (!args.getControlFiles().isEmpty()) {
				getRTstopCounts(counts, 2, treeC, snpNode, seq, start, end, upstream);
			}
		}
		return counts;
	}
	
	private int[][][] getRTstopCounts(int[][][] counts, int countIndex, IntervalTree<List<String>> tree, Node<Mutation> snpNode,
			String seq, int start, int end, boolean upstream) {
		boolean base_flag = true;
		controlReads = 0;
		for (Iterator<Node<List<String>>> nodes = tree.overlappers(start, end); nodes.hasNext();) {
			Node<List<String>> node = nodes.next();
			base_flag = args.isVaildBase(seq == null ? 0 : 
				Character.toUpperCase(upstream ? seq.charAt(node.getStart() - 2) : seq.charAt(node.getEnd())));
			char ref = seq == null ? 0 : 
				Character.toUpperCase(upstream ? seq.charAt(node.getStart() - 1) : seq.charAt(node.getEnd() - 1));
			int stop_index = upstream ? (node.getStart() > start ? node.getStart() - start - 1 : -1)
					: (node.getEnd() < end ? node.getEnd() - start + 1 : -1);
			for (String seg : node.getValue()) {
				if ((upstream && seg.startsWith("S")) || (!upstream && seg.endsWith("S")) 
						|| (ref != 0 && ref != Character.toUpperCase(upstream ? seg.charAt(0) : seg.charAt(seg.length() - 1)))) {
					continue;
				}
				++controlReads;
				int baseIndex = CommonMethod.randInt(2) == 0 ? 0 : 1;
				if (!permutate) {
					int startAlign = seg.startsWith("S") ? 1 : 0;
					char snpAlign = snpNode == null || snpNode.getStart() > node.getEnd() || snpNode.getEnd() < node.getStart() ?
							0 : Character.toUpperCase(seg.charAt(snpNode.getStart() - node.getStart() + startAlign));
					baseIndex = getBaseIndex(snpAlign);
				}
				if (stop_index >= 0 && base_flag) {
					++counts[baseIndex][countIndex][stop_index];
				}
				for (int i = Math.max(node.getStart(), start); i <= Math.min(node.getEnd(), end); ++i) {
					++counts[baseIndex][countIndex + 1][i - start];
				}
			}
		}
		return counts;
	}
	
//	private int[][][] getRTstopCountsPre(int[][][] counts, int countIndex, IntervalTree<List<String>> tree, Node<Mutation> snpNode,
//			String seq, int start, int end, boolean upstream) {
//		boolean base_flag = true;
//		counts = counts != null ? counts : new int[snpNode == null ? 1 : 5][args.getControlFiles().isEmpty() ? 2 : 4][end - start + 1];
//		controlReads = 0;
//		for (Iterator<Node<List<String>>> nodes = tree.overlappers(start, end); nodes.hasNext();) {
//			Node<List<String>> node = nodes.next();
//			base_flag = isVaildBase(seq == null ? 0 : 
//				Character.toUpperCase(upstream ? seq.charAt(node.getStart() - 2) : seq.charAt(node.getEnd())));
//			char ref = seq == null ? 0 : 
//				Character.toUpperCase(upstream ? seq.charAt(node.getStart() - 1) : seq.charAt(node.getEnd() - 1));
//			int[] valid = getValidStopAlign(node, snpNode, ref, upstream);
//			for (int i = 0; i < valid.length; ++i) {
//				controlReads += valid[i];
//			}
//			int stop_index = upstream ? (node.getStart() > start ? node.getStart() - start - 1 : -1)
//					: (node.getEnd() < end ? node.getEnd() - start + 1 : -1);
//			if (stop_index >= 0 && base_flag) {
//				for (int i = 0; i < counts.length; ++i) {
//					counts[i][countIndex][stop_index] += valid[i];
//				}
//			}
//			for (int i = Math.max(node.getStart(), start); i <= Math.min(node.getEnd(), end); ++i) {
//				for (int j = 0; j < counts.length; ++j) {
//					counts[j][countIndex + 1][i - start] += valid[j];
//				}
//			}
//		}
//		return counts;
//	}
	
//	private int[] getValidStopAlign(Node<List<String>> alignNode, Node<Mutation> snpNode, char ref, boolean upstream) {
//		int[] valid = {0, 0, 0, 0, 0};
//		for (String seg : alignNode.getValue()) {
//			if (permutate) {
//				int baseIndex = Math.random() < 0.5 ? 0 : 1;
//				valid[baseIndex] += ref == 0 || ref == Character.toUpperCase(upstream ? seg.charAt(0) : seg.charAt(seg.length() - 1)) ? 1 : 0;
//				continue;
//			}
//			if ((upstream && seg.startsWith("S")) || (!upstream && seg.endsWith("S"))) {
//				continue;
//			}
//			int startAlign = seg.charAt(0) == 'S' ? 1 : 0;
//			char snpAlign = snpNode == null || snpNode.getStart() > alignNode.getEnd() || snpNode.getEnd() < alignNode.getStart() ?
//					0 : Character.toUpperCase(seg.charAt(snpNode.getStart() - alignNode.getStart() - startAlign));
//			int baseIndex = getBaseIndex(snpAlign);
//			valid[baseIndex] += ref == 0 || ref == Character.toUpperCase(upstream ? seg.charAt(0) : seg.charAt(seg.length() - 1)) ? 1 : 0;
//		}
//		return valid;
//	}
	
//	private int[][][][] getRTmutCounts(int[][][][] counts,IntervalTree<List<String>> treeT, IntervalTree<List<String>> treeC, 
//			Node<Mutation> snpNode, String seq, int start, int end) {
//		int snpAlt = snpNode == null ? 1 : snpNode.getValue().getSNP().length;
//		int snpNum = snpNode == null ? 1 : 5;
//		int tc = args.getControlFiles().isEmpty() ? 2 : 4;
//		counts = counts != null ? counts : new int[snpAlt][snpNum][tc][end - start + 1];
//		for (int i = 0; i < counts.length; ++i) {
//			getRTmutCounts(counts[i], 0, treeT, snpNode, seq, start, end);
//			treatReads = controlReads;
//			if (!InParam.getParams().getControlFiles().isEmpty()) {
//				getRTmutCounts(counts[i], 2, treeC, snpNode, seq, start, end);
//			}
//		}
//		return counts;
//	}
	
	private int[][][] getRTmutCounts(int[][][] counts, int countIndex, IntervalTree<List<String>> tree, Node<Mutation> snpNode,
			String seq, int start, int end) {
		controlReads = 0;
		for (Iterator<Node<List<String>>> nodes = tree.overlappers(start, end); nodes.hasNext();) {
			Node<List<String>> node = nodes.next();
			int overStart = Math.max(node.getStart(), start);
			int overEnd = Math.min(node.getEnd(), end);
			for (String seg : node.getValue()) {
				++controlReads;
				int startAlign = seg.charAt(0) == 'S' ? 1 : 0;
				char snpAlign = snpNode == null || snpNode.getStart() > node.getEnd() || snpNode.getEnd() < node.getStart() ?
						0 : Character.toUpperCase(seg.charAt(snpNode.getStart() - node.getStart() + startAlign));
				int baseIndex = CommonMethod.randInt(2) == 0 ? 0 : 1;
				baseIndex = permutate ? baseIndex : getBaseIndex(snpAlign);
				for (int i = overStart; i <= overEnd; ++i) {
					char c = seq == null ? 0 : Character.toUpperCase(seq.charAt(i - 1));
					counts[baseIndex][countIndex][i - start] += args.isVaildBase(c) && 
							c != Character.toUpperCase(seg.charAt(i - node.getStart() - startAlign)) ? 1 : 0;
					++counts[baseIndex][countIndex + 1][i - start];
				}
			}
		}
		return counts;
	}
	
//	private double[] testRTstopSNPScores(IntervalTree<List<String>> treat_tree, IntervalTree<List<String>> control_tree,
//			Node<Mutation> mutNode, int window, int time, double[] socres, boolean upstream) {
//		int start = mutNode.getStart() - window + 1;
//		int end = mutNode.getEnd() + window - 1;
//		double[] permutation = new double[socres.length];
//		for (int i = 0; i < permutation.length; ++i) {
//			permutation[i] /= time;
//		}
//		return permutation;
//	}
//	
//	private void compareScores(double[] permutation, double[] socres, double[] stopScroes) {
//		for (int i = 0; i < permutation.length; ++i) {
//			if (stopScroes[i] >= socres[i]) {
//				++permutation[i];
//			}
//		}
//	}
//
//	private int[][] getRTstopSNPCounts(IntervalTree<List<String>> treat_tree, IntervalTree<List<String>> control_tree, String seq,
//			Node<Mutation> mutNode, int window, boolean alt, boolean upstream) {
//		int[][] countsT = getRTstopSNPCounts(treat_tree, seq, mutNode, window, alt, upstream);
//		if (control_tree != null) {
//			int[][] countsC = getRTstopSNPCounts(control_tree, seq, mutNode, window, alt, upstream);
//			int[][] counts = new int[4][];
//			counts[0] = countsT[0];
//			counts[1] = countsT[1];
//			counts[2] = countsC[0];
//			counts[3] = countsC[1];
//			return counts;
//		}
//		return countsT;
//	}
//	
//	private int[][] getRTstopSNPCounts(IntervalTree<List<String>> tree, String seq, Node<Mutation> mutNode, int window,
//			boolean alt, boolean upstream) {
//		int start = mutNode.getStart() - window + 1;
//		int end = mutNode.getEnd() + window - 1;
//		int[][] counts = new int[2][end - start + 1];
//		int[] count = counts[1];
//		int[] stop = counts[0];
//		boolean base_flag = true;
//		for (Iterator<Node<List<String>>> nodes = tree.overlappers(mutNode.getStart(), mutNode.getEnd()); nodes.hasNext();) {
//			Node<List<String>> node = nodes.next();
//			char ref = seq == null ? 0 : 
//				Character.toUpperCase(upstream ? seq.charAt(node.getStart() - 2) : seq.charAt(node.getEnd()));
//			base_flag = isVaildBase(ref);
//			int valid = getSNPStopAlign(node.getValue(), ref, mutNode.getValue(), 
//					mutNode.getStart() - node.getStart(), alt, upstream);
//			int stop_index = upstream ? (node.getStart() > start ? node.getStart() - start - 1 : -1)
//					: (node.getEnd() < end ? node.getEnd() - start + 1 : -1);
//			if (stop_index > 0 && base_flag) {
//				stop[stop_index] += valid;
//			}
//			for (int i = Math.max(node.getStart(), start); i <= Math.min(node.getEnd(), end); ++i) {
//				count[i - start] += valid;
//			}
//		}
//		return counts;
//	}
//	
//	private int getSNPStopAlign(List<String> alignments, char ref, Mutation mut, int index, boolean alt, boolean upstream) {
//		int valid = 0;
//		for (String seg : alignments) {
//			char segBase = Character.toUpperCase(upstream ? seg.charAt(0) : seg.charAt(seg.length() - 1));
//			char mutBase = seg.charAt(index);
//			valid += (ref == 0 || segBase == ref) && (alt ? mut.isMut(mutBase) : mut.getRef() == mutBase) ? 1 : 0;
//		}
//		return valid;
//	}
	
	private double[] normalizeFac(int[][] counts) {
		double[] facts = new double[counts.length];
		if (counts[0].length == 0) {
			return facts;
		}
		for (int i = 0; i < counts.length; ++i) {
			int[] copy = counts[i].clone();
			Arrays.sort(copy);
			int start = copy.length - 1 - copy.length * 2 / 10;
			int end = copy.length - 1 - copy.length / 10;
			facts[i] = (start + end) % 2 == 0 ? (double) copy[(start + end) / 2] 
					: (double) (copy[(start + end) / 2] + copy[(start + end + 1) / 2]) / 2.0;
		}
		return facts;
	}
	
	private int[][] baseCountFilter(int[][] counts, String seq) {
		if (seq == null || counts == null || counts.length == 0 || seq.length() < counts[0].length) {
			return counts;
		}
		List<int[]> filtList = new ArrayList<>(); 
		for (int i = 0; i < counts[0].length; ++i) {
			if (args.isVaildBase(seq.charAt(i))) {
				int[] baseCount = new int[counts.length];
				for (int j = 0; j < counts.length; ++j) {
					baseCount[j] = counts[j][i];
				}
				filtList.add(baseCount);
			}
		}
		int[][] filtCount = new int[counts.length][filtList.size()];
		for (int i = 0; i < filtList.size(); ++i) {
			for (int j = 0; j < counts.length; ++j) {
				filtCount[j][i] = filtList.get(i)[j];
			}
		}
		return filtCount;
	}
	
	private double[] scoreFilter(double[] scores, String seq) {
		if (seq == null || scores == null) {
			return scores;
		}
		List<Double> filtList = new ArrayList<>();
		for (int i = 0; i < scores.length; ++i) {
			if (args.isVaildBase(seq.charAt(i))) {
				filtList.add(scores[i]);
			}
		}
		double[] fillScore = new double[filtList.size()];
		for (int i = 0; i < filtList.size(); ++i) {
			fillScore[i] = filtList.get(i);
		}
		return fillScore;
	}
	
	private double[] filtScoreFill(double[] scores, String seq) {
		if (seq == null || scores == null) {
			return scores;
		}
		List<Double> filtList = new ArrayList<>();
		int j = -1;
		for (int i = 0; i < seq.length(); ++i) {
			if (args.isVaildBase(seq.charAt(i))) {
				filtList.add(scores[++j]);
			}
			else {
				filtList.add(-1.0);
			}
		}
		double[] fillScore = new double[filtList.size()];
		for (int i = 0; i < filtList.size(); ++i) {
			fillScore[i] = filtList.get(i);
		}
		return fillScore;
	}

	
	private double[] normalize2(double[] nums) {
		double[] copy = nums.clone();
		Arrays.sort(copy);
		int start = 0;
		for (int i = 0; i < copy.length; ++i) {
			if (Double.isFinite(copy[i]) || copy[i] >= 0.0) {
				start = i;
				break;
			}
		}
		int end = 0;
		for (int i = copy.length - 1; i >= 0; --i) {
			if (Double.isFinite(copy[i])) {
				end = i;
				break;
			}
		}
		if (end - start + 1 < 20) {
			return null;
		}
		start = (int) (win_fac * (end - start + 1));
		double max = copy[copy.length - start - 1];
		double min = Math.max(0.0, copy[start]);
		if (max <= min) {
			return null;
		}
		for (int i = 0; i < nums.length; ++i) {
			copy[i] = nums[i] < 0.0 ? -1.0 : nums[i] > max ? 1.0 : (nums[i] < min ? 0.0 : (nums[i] - min) / (max - min));
		}
		return copy;
	}
	
	private String getOutGeneScores(int[][][] counts, int[][][] counts2, String chr, Genome genome, 
			int start, int end, boolean positive) {
		StringBuilder sb = new StringBuilder();
		StringBuilder sbd = new StringBuilder();
		sb.append(chr);
		sb.append('\t');
		sb.append(start - 1);
		sb.append('\t');
		sb.append(end);
		String subSeq = getChrSeq(genome, chr) == null ? null : getChrSeq(genome, chr).substring(start - 1, end);
		double[] scores1 = getWinShapeScores(counts[0], start, end, subSeq);
		double[] scores2 = getWinShapeScores(counts2[0], start, end, subSeq);

		
//		double[] scoresShare1 = getWinShapeScores(counts[getBaseIndex(snpNode.getValue().getRef())], start, end, subSeq);
//		double[] scoresShare2 = getWinShapeScores(counts[getBaseIndex(snp)], start, end, subSeq);
		double esdcShareScore = getEsdcScore(scores1, scores2);
		sbd.append('\t');
		sbd.append(String.format("%.4f", esdcShareScore));
		int times = args.getPermutate();
		double esdcSharePer = 0.0;
		double[] diffShareScores = CommonMethod.subNoNegAbsArray(scores1, scores2);
		double[] diffSharePer = new double[diffShareScores.length];
		double[] diffShareTotal = new double[diffShareScores.length];
		int esdcTime = 0;
		permutate = true;
		for (int i = 0; i < times; ++i) {
			double[][] testScores = getTestSnpScores(null, getChrSeq(genome, chr), start, end, positive);
			double testScore = getEsdcScore(testScores[0], testScores[1]);
			if (testScore >= 0.0) {
				esdcSharePer += testScore >= esdcShareScore ? 1.0 : 0.0;
				++esdcTime;
			}
			double[] diffScoresPer = CommonMethod.subNoNegAbsArray(testScores[0], testScores[1]);
	 		for (int j = 0; j < diffScoresPer.length; ++j) {
	 			diffShareTotal[j] += diffScoresPer[j] >= 0.0 ? 1.0 : 0.0;
	 			diffSharePer[j] += diffScoresPer[j] >= diffShareScores[j] ? 1.0 : 0.0;
	 		}
		}
		permutate = false;
		
		sbd.append('\t');
		sbd.append(String.format("%.4f", esdcTime == 0 || esdcShareScore < 0.0 ? -1.0 : esdcSharePer / esdcTime));
		double esdcRepScore = -2.0;
		if (replicate) {
			double[][] repScores = null;
			int[][][] repCounts = getCounts(null, treeTR1, treeC, null, getChrSeq(genome, chr), start, end, positive);
			repScores = new double[2][];
			repScores[0] = getWinShapeScores(CommonMethod.addArray(repCounts), start, end, subSeq);
			repCounts = getCounts(null, treeTR2, treeC, null, getChrSeq(genome, chr), start, end, positive);
			repScores[1] = getWinShapeScores(CommonMethod.addArray(repCounts), start, end, subSeq);

//			double kstRepScore = getKstScore(repScores[0], repScores[1]);
//			sb.append('\t');
//			sb.append(String.format("%.4f", kstRepScore));
			esdcRepScore = getEsdcScore(repScores[0], repScores[1]);
		}
		sbd.append('\t');
		sbd.append(esdcRepScore < -1.5 ? "-" : String.format("%.4f", esdcRepScore));
		sbd.append('\t');
		sbd.append(subSeq);
		sbd.append('\t');
		sbd.append(Method.toString(scores1, "%.3f", ","));
		sbd.append('|');
		sbd.append(Method.toString(scores2, "%.3f", ","));
		for (int i = 0; i < diffSharePer.length; ++i) {
			diffSharePer[i] = diffShareTotal[i] < 0.5 || diffShareScores[i] < 0.0 ? -1.0 : diffSharePer[i] / diffShareTotal[i];
		}
		sbd.append('\t');
		sbd.append(Method.toString(diffSharePer, "%.3f", ","));
		
		
		esdcSharePer = esdcSharePer < 0.5 ? Math.log10(times) + 1 : -Math.log10(esdcSharePer / esdcTime);
		boolean negShareFlag = CommonMethod.retainNoNeg(scores1).length < 2 || esdcShareScore < 0.0 ||
			CommonMethod.retainNoNeg(scores2).length < 2 || Math.abs(esdcRepScore + 1.0) < 0.5 || esdcTime == 0;
		negShareFlag |= (Math.abs(esdcShareScore) > 0.000001 && esdcShareScore <= esdcRepScore);

		double esdcShareAdjustScore = negShareFlag ? -1.0 :	Math.abs(esdcShareScore * esdcSharePer);
		sb.append('\t');
		sb.append(esdcShareAdjustScore);
		
		Transcript script = genome.getChr(chr).getScrpitTree().minOverlapper(start, end) != null ? 
				genome.getChr(chr).getScrpitTree().minOverlapper(start, end).getValue() : null;
		sb.append('\t');
		sb.append(script != null ? script.getGene().getStrand() : '.');
		sb.append('\t');
		sb.append(script != null ? script.getGene().getGene_id() : "None");
		sb.append('\t');
		sb.append(script != null ? script.getScript_id() : "None");
		sb.append('\t');
		sb.append(script != null ? script.getScript_type() : "None");
		sb.append('|');
		sb.append(script != null ? script.getRegionFeature((start + end) / 2) : "None");
		sb.append('\t');
		sb.append(start);
		sb.append(',');
		sb.append(end);
		sb.append('|');
		sb.append(script != null ? start - script.getStart() : 0);
		sb.append(',');
		sb.append(script != null ? end - script.getStart() : 0);
//		if (genomeSeq != null) {
//			sb.append('\t');
//			sb.append(genomeSeq.substring(start - 1, end));
//		}
//		sb.append('\t');
//		sb.append(Method.toString(scores1, "%.3f", ","));
//		sb.append('\t');
//		sb.append(Method.toString(scores2, "%.3f", ","));
		
//		double snpScore = getKstScore(scores1, scores2);
//		sb.append('\t');
//		sb.append(String.format("%.4f", snpScore));
//		sb.append('\t');
//		sb.append(String.format("%.4f", getTestKstScores(snpScore, args.getPermutate(), null, genomeSeq, start, end, positive)));
//		
//		snpScore = getEsdcScore(scores1, scores2);
//		sb.append('\t');
//		sb.append(String.format("%.4f", snpScore));
//		sb.append('\t');
//		sb.append(String.format("%.4f", getTestEsdcScores(snpScore, args.getPermutate(), null, genomeSeq, start, end, positive)));
//		
//		double[] scoresDiff = CommonMethod.subAbsArray(scores1, scores2);
//		sb.append('\t');
//		sb.append(String.format("%.4f", CommonMethod.equalArray(scores1, scores2) ? 1.0 : 
//			CommonMethod.combinePvalue(getTestShapeScores(scoresDiff, args.getPermutate(), null, genomeSeq, start, end, positive))));
		
		return(sb.toString());
	}
	
	private String[] getOutSnpScores(int[][][] counts, int[][][] counts2, String chr, String genomeSeq, Node<Mutation> snpNode,
			char snp, char homo, int start, int end, boolean positive) {
		StringBuilder sb = new StringBuilder();
		StringBuilder sbd = new StringBuilder();
		sb.append(chr);
//		sb.append('\t');
//		sb.append(snpNode.getStart() - 1);
		sb.append('\t');
		sb.append(snpNode.getEnd());
//		sb.append('\t');
//		sb.append(snpNode.getValue().toSimpleString(1.0, 1.0));
		sb.append('\t');
		sb.append(snpNode.getValue().getRef());
		sb.append('\t');
		sb.append(snpNode.getValue().getSNPString());
		sb.append('\t');
		sb.append(homo);
		sb.append('\t');
		sb.append(snpNode.getValue().getDescription());
		sbd.append(snpNode.getValue().getDescription());
		
		
		String subSeq = genomeSeq != null ? genomeSeq.substring(start - 1, end) : null;
		double[] scoresShare1 = getWinShapeScores(CommonMethod.addArray(counts), start, end, subSeq);
		double[] scoresShare2 = getWinShapeScores(CommonMethod.addArray(counts2), start, end, subSeq);
//		double[] scoresShare1 = getWinShapeScores(counts[getBaseIndex(snpNode.getValue().getRef())], start, end, subSeq);
//		double[] scoresShare2 = getWinShapeScores(counts[getBaseIndex(snp)], start, end, subSeq);
		double esdcShareScore = getEsdcScore(scoresShare1, scoresShare2);
		sbd.append('\t');
		sbd.append(String.format("%.4f", esdcShareScore));
		int times = args.getPermutate();
		double esdcSharePer = 0.0;
		double[] diffShareScores = CommonMethod.subNoNegAbsArray(scoresShare1, scoresShare2);
		double[] diffSharePer = new double[diffShareScores.length];
		double[] diffShareTotal = new double[diffShareScores.length];
		int esdcTime = 0;
		permutate = true;
		for (int i = 0; i < times; ++i) {
			double[][] testScores = getTestSnpScores(snpNode, genomeSeq, start, end, positive);
			double testScore = getEsdcScore(testScores[0], testScores[1]);
			if (testScore >= 0.0) {
				esdcSharePer += testScore >= esdcShareScore ? 1.0 : 0.0;
				++esdcTime;
			}
			double[] diffScoresPer = CommonMethod.subNoNegAbsArray(testScores[0], testScores[1]);
	 		for (int j = 0; j < diffScoresPer.length; ++j) {
	 			diffShareTotal[j] += diffScoresPer[j] >= 0.0 ? 1.0 : 0.0;
	 			diffSharePer[j] += diffScoresPer[j] >= diffShareScores[j] ? 1.0 : 0.0;
	 		}
		}
		permutate = false;
		
		sbd.append('\t');
		sbd.append(String.format("%.4f", esdcTime == 0 || esdcShareScore < 0.0 ? -1.0 : esdcSharePer / esdcTime));
		double esdcRepScore = -2.0;
		if (replicate) {
			double[][] repScores = null;
			int[][][] repCounts = getCounts(null, treeTR1, treeC, null, genomeSeq, start, end, positive);
			repScores = new double[2][];
			repScores[0] = getWinShapeScores(CommonMethod.addArray(repCounts), start, end, subSeq);
			repCounts = getCounts(null, treeTR2, treeC, null, genomeSeq, start, end, positive);
			repScores[1] = getWinShapeScores(CommonMethod.addArray(repCounts), start, end, subSeq);

			esdcRepScore = getEsdcScore(repScores[0], repScores[1]);
		}
		
		int homoCount = 0;
		int totalCount = 0;
		for (int i = 0; i < counts[0].length; ++i) {
			homoCount += counts[getBaseIndex(homo)][i][snpNode.getStart() - start];
			totalCount += counts[getBaseIndex(snp)][i][snpNode.getStart() - start];
			totalCount += counts[getBaseIndex(snpNode.getValue().getRef())][i][snpNode.getStart() - start];
		}
		sbd.append('\t');
		sbd.append(String.format("%.4f", (double) homoCount / (double) totalCount));
		sbd.append('\t');
		sbd.append(esdcRepScore < -1.5 ? "-" : String.format("%.4f", esdcRepScore));
		sbd.append('\t');
		sbd.append(counts[getBaseIndex('A')][1][snpNode.getStart() - start]);
		sbd.append('|');
		sbd.append(counts[getBaseIndex('T')][1][snpNode.getStart() - start]);
		sbd.append('|');
		sbd.append(counts[getBaseIndex('C')][1][snpNode.getStart() - start]);
		sbd.append('|');
		sbd.append(counts[getBaseIndex('G')][1][snpNode.getStart() - start]);
		sbd.append('\t');
		sbd.append(counts2[getBaseIndex('A')][1][snpNode.getStart() - start]);
		sbd.append('|');
		sbd.append(counts2[getBaseIndex('T')][1][snpNode.getStart() - start]);
		sbd.append('|');
		sbd.append(counts2[getBaseIndex('C')][1][snpNode.getStart() - start]);
		sbd.append('|');
		sbd.append(counts2[getBaseIndex('G')][1][snpNode.getStart() - start]);
		sbd.append('\t');
		sbd.append(subSeq);
		sbd.append('\t');
		sbd.append(Method.toString(scoresShare1, "%.3f", ","));
		sbd.append('|');
		sbd.append(Method.toString(scoresShare2, "%.3f", ","));
		for (int i = 0; i < diffSharePer.length; ++i) {
			diffSharePer[i] = diffShareTotal[i] < 0.5 || diffShareScores[i] < 0.0 ? -1.0 : diffSharePer[i] / diffShareTotal[i];
		}
		sbd.append('\t');
		sbd.append(Method.toString(diffSharePer, "%.3f", ","));
		sbd.append('\t');
//		sbd.append(genomeSeq.substring(Math.max(0, snpNode.getEnd() - 201), Math.min(genomeSeq.length(), snpNode.getEnd() + 200)));
		
		
		esdcSharePer = esdcSharePer < 0.5 ? Math.log10(times) + 1 : -Math.log10(esdcSharePer / esdcTime);
		boolean negShareFlag = CommonMethod.retainNoNeg(scoresShare1).length < 2 || esdcShareScore < 0.0 ||
			CommonMethod.retainNoNeg(scoresShare2).length < 2 || Math.abs(esdcRepScore + 1.0) < 0.5 || esdcTime == 0;
		negShareFlag |= (Math.abs(esdcShareScore) > 0.000001 && esdcShareScore <= esdcRepScore);

		double esdcShareAdjustScore = negShareFlag ? -1.0 : Math.abs(esdcShareScore * esdcSharePer)
				* (1.0 - Math.log((double) homoCount / (double) totalCount));
		sb.append('\t');
		sb.append(esdcShareAdjustScore);
		
		sb.append('\t');
		sb.append(snpNode.getValue().getScript().getGene().getStrand());
		sb.append('\t');
		sb.append(snpNode.getValue().getScript().getGene().getGene_id());
		sb.append('\t');
		sb.append(snpNode.getValue().getScript().getScript_id());
		sb.append('\t');
		sb.append(snpNode.getValue().getScript().getScript_type());
		sb.append('|');
		sb.append(snpNode.getValue().getScript().getRegionFeature(snpNode.getEnd()));
		sb.append('\t');
		sb.append(start);
		sb.append(',');
		sb.append(end);
		sb.append('|');
		sb.append(start - snpNode.getValue().getScript().getStart());
		sb.append(',');
		sb.append(end - snpNode.getValue().getScript().getStart());
		return(new String[] {sb.toString(), sbd.toString()});
	}
	
	private String[] getOutSnpScores(int[][][] counts, String chr, String genomeSeq, Node<Mutation> snpNode, char snp,
			int start, int end, boolean positive) {
		StringBuilder sb = new StringBuilder();
		StringBuilder sbd = new StringBuilder();
		sb.append(chr);
//		sb.append('\t');
//		sb.append(snpNode.getStart() - 1);
		sb.append('\t');
		sb.append(snpNode.getEnd());
//		sb.append('\t');
//		sb.append(snpNode.getValue().toSimpleString(1.0, 1.0));
		sb.append('\t');
		sb.append(snpNode.getValue().getRef());
		sb.append('\t');
		sb.append(snpNode.getValue().getSNPString());
		sb.append('\t');
		sb.append(snpNode.getValue().getDescription());
		sbd.append(snpNode.getValue().getDescription());
		
		
		String subSeq = genomeSeq != null ? genomeSeq.substring(start - 1, end) : null;
		double snpCount = snpNode.getValue().getCount(snpNode.getValue().getRef(), true)
				+ snpNode.getValue().getCount(snpNode.getValue().getRef(), false) 
				+ snpNode.getValue().getCount(snp, true)
				+ snpNode.getValue().getCount(snp, false);
		double[] scoresShare1 = getWinShapeScores(CommonMethod.addArrayWithout(counts, getBaseIndex(snp)), start, end, subSeq);
		double[] scoresShare2 = getWinShapeScores(CommonMethod.addArrayWithout(
				counts, getBaseIndex(snpNode.getValue().getRef())), start, end, subSeq);
//		double[] scoresShare1 = getWinShapeScores(counts[getBaseIndex(snpNode.getValue().getRef())], start, end, subSeq);
//		double[] scoresShare2 = getWinShapeScores(counts[getBaseIndex(snp)], start, end, subSeq);
		double esdcShareScore = getEsdcScore(scoresShare1, scoresShare2);
		sbd.append('\t');
		sbd.append(String.format("%.4f", esdcShareScore));
		int times = args.getPermutate();
		double esdcSharePer = 0.0;
		double[] diffShareScores = CommonMethod.subNoNegAbsArray(scoresShare1, scoresShare2);
		double[] diffSharePer = new double[diffShareScores.length];
		double[] diffShareTotal = new double[diffShareScores.length];
		int esdcTime = 0;
		permutate = true;
		for (int i = 0; i < times; ++i) {
			double[][] testScores = getTestSnpScores(snpNode, genomeSeq, start, end, positive);
			double testScore = getEsdcScore(testScores[0], testScores[1]);
			if (testScore >= 0.0) {
				esdcSharePer += testScore >= esdcShareScore ? 1.0 : 0.0;
				++esdcTime;
			}
			double[] diffScoresPer = CommonMethod.subNoNegAbsArray(testScores[0], testScores[1]);
	 		for (int j = 0; j < diffScoresPer.length; ++j) {
	 			diffShareTotal[j] += diffScoresPer[j] >= 0.0 ? 1.0 : 0.0;
	 			diffSharePer[j] += diffScoresPer[j] >= diffShareScores[j] ? 1.0 : 0.0;
	 		}
		}
		permutate = false;
		
		sbd.append('\t');
		sbd.append(String.format("%.4f", esdcTime == 0 || esdcShareScore < 0.0 ? -1.0 : esdcSharePer / esdcTime));
		sbd.append('\t');
		sbd.append(String.format("%.4f", snpCount / (treatReads + controlReads)));
		double esdcRepScore = -2.0;
		if (replicate) {
			double[][] repScores = null;
			int[][][] repCounts = getCounts(null, treeTR1, treeC, null, genomeSeq, start, end, positive);
			repScores = new double[2][];
			repScores[0] = getWinShapeScores(CommonMethod.addArray(repCounts), start, end, subSeq);
			repCounts = getCounts(null, treeTR2, treeC, null, genomeSeq, start, end, positive);
			repScores[1] = getWinShapeScores(CommonMethod.addArray(repCounts), start, end, subSeq);

//			double kstRepScore = getKstScore(repScores[0], repScores[1]);
//			sb.append('\t');
//			sb.append(String.format("%.4f", kstRepScore));
			esdcRepScore = getEsdcScore(repScores[0], repScores[1]);
		}
		sbd.append('\t');
		sbd.append(esdcRepScore < -1.5 ? "-" : String.format("%.4f", esdcRepScore));
		sbd.append('\t');
		sbd.append(snpNode.getValue().getIpCount('A'));
		sbd.append('|');
		sbd.append(snpNode.getValue().getIpCount('T'));
		sbd.append('|');
		sbd.append(snpNode.getValue().getIpCount('C'));
		sbd.append('|');
		sbd.append(snpNode.getValue().getIpCount('G'));
		sbd.append('\t');
		sbd.append(snpNode.getValue().getInputCount('A'));
		sbd.append('|');
		sbd.append(snpNode.getValue().getInputCount('T'));
		sbd.append('|');
		sbd.append(snpNode.getValue().getInputCount('C'));
		sbd.append('|');
		sbd.append(snpNode.getValue().getInputCount('G'));
		sbd.append('\t');
		sbd.append(subSeq);
		sbd.append('\t');
		sbd.append(Method.toString(scoresShare1, "%.3f", ","));
		sbd.append('|');
		sbd.append(Method.toString(scoresShare2, "%.3f", ","));
		for (int i = 0; i < diffSharePer.length; ++i) {
			diffSharePer[i] = diffShareTotal[i] < 0.5 || diffShareScores[i] < 0.0 ? -1.0 : diffSharePer[i] / diffShareTotal[i];
		}
		sbd.append('\t');
		sbd.append(Method.toString(diffSharePer, "%.3f", ","));
		
		
		esdcSharePer = esdcSharePer < 0.5 ? Math.log10(times) + 1 : -Math.log10(esdcSharePer / esdcTime);
		boolean negShareFlag = CommonMethod.retainNoNeg(scoresShare1).length < 2 || esdcShareScore < 0.0 ||
			CommonMethod.retainNoNeg(scoresShare2).length < 2 || Math.abs(esdcRepScore + 1.0) < 0.5 || esdcTime == 0;
		negShareFlag |= (Math.abs(esdcShareScore) > 0.000001 && esdcShareScore <= esdcRepScore);

		double esdcShareAdjustScore = negShareFlag ? -1.0 : 
			Math.abs(esdcShareScore * esdcSharePer * (1 - Math.log10(snpCount / (treatReads + controlReads))));
		sb.append('\t');
		sb.append(esdcShareAdjustScore);
		
		sb.append('\t');
		sb.append(snpNode.getValue().getScript().getGene().getStrand());
		sb.append('\t');
		sb.append(snpNode.getValue().getScript().getGene().getGene_id());
		sb.append('\t');
		sb.append(snpNode.getValue().getScript().getScript_id());
		sb.append('\t');
		sb.append(snpNode.getValue().getScript().getScript_type());
		sb.append('|');
		sb.append(snpNode.getValue().getScript().getRegionFeature(snpNode.getEnd()));
		sb.append('\t');
		sb.append(start);
		sb.append(',');
		sb.append(end);
		sb.append('|');
		sb.append(start - snpNode.getValue().getScript().getStart());
		sb.append(',');
		sb.append(end - snpNode.getValue().getScript().getStart());
//		
//		
//		
//		
//		
//		
//		
//		
//		String feature = snpNode.getValue().getScript().getRegionFeature(snpNode.getEnd());
//		sb.append('\t');
//		sb.append(feature);
//		
//		
//		
//		sb.append('\t');
//		sb.append(start);
//		sb.append('\t');
//		sb.append(end);
//		for (int k = 0; k < counts.length; ++k) {
//			sb.append('\t');
//			sb.append(k);
//			for (int i = 0; i < counts[k].length; ++i) {
//				sb.append('\t');
//				sb.append(Arrays.toString(counts[k][i]));
//			}
//		}
//		
////		String subSeq = null;
//		if (genomeSeq != null) {
//			sb.append('\t');
//			subSeq = genomeSeq.substring(start - 1, end);
//			sb.append(subSeq);
//			sb.append('\t');
//			sb.append(snpNode.getValue().getScript().getExonGCContent(feature, snpNode.getEnd()));
//			sb.append('|');
//			sb.append(BamMethod.getGCContent(genomeSeq, snpNode.getEnd(), 1000));
//		}
//		
////		double[] scores1 = getWinShapeScores(counts[getBaseIndex(snpNode.getValue().getRef())], start, end, subSeq);
////		double[] scores2 = getWinShapeScores(counts[getBaseIndex(snp)], start, end, subSeq);
//		
////		double kstScore = getKstScore(scores1, scores2);
////		sb.append('\t');
////		sb.append(String.format("%.4f", kstScore));
////		double esdcScore = getEsdcScore(scores1, scores2);
////		sb.append('\t');
////		sb.append(String.format("%.4f", esdcScore));
//		
////		double[] scoresShare1 = getWinShapeScores(CommonMethod.addArrayWithout(counts, getBaseIndex(snp)), start, end, subSeq);
////		double[] scoresShare2 = getWinShapeScores(CommonMethod.addArrayWithout(
////				counts, getBaseIndex(snpNode.getValue().getRef())), start, end, subSeq);
//		
//		double kstShareScore = getKstScore(scoresShare1, scoresShare2);
//		sb.append('\t');
//		sb.append(String.format("%.4f", kstShareScore));	
//
//		double snpCount = snpNode.getValue().getCount(snpNode.getValue().getRef(), true)
//			+ snpNode.getValue().getCount(snpNode.getValue().getRef(), false) 
//			+ snpNode.getValue().getCount(snp, true)
//			+ snpNode.getValue().getCount(snp, false);
//		sb.append('\t');
//		sb.append(String.format("%.4f", snpCount / (treatReads + controlReads)));
//		sb.append('\t');
//		sb.append(String.format("%.4f", kstShareScore * snpCount / (treatReads + controlReads)));
//		
//		double esdcShareScore = getEsdcScore(scoresShare1, scoresShare2);
//		sb.append('\t');
//		sb.append(String.format("%.4f", esdcShareScore));
//		
//		double esdcRepScore = -1.0;
//		if (replicate) {
//			double[][] repScores = null;
//			int[][][] repCounts = getCounts(null, treeT1, treeC, null, genomeSeq, start, end, positive);
//			repScores = new double[2][];
//			repScores[0] = getWinShapeScores(CommonMethod.addArray(repCounts), start, end, subSeq);
//			repCounts = getCounts(null, treeT2, treeC, null, genomeSeq, start, end, positive);
//			repScores[1] = getWinShapeScores(CommonMethod.addArray(repCounts), start, end, subSeq);
//
////			double kstRepScore = getKstScore(repScores[0], repScores[1]);
////			sb.append('\t');
////			sb.append(String.format("%.4f", kstRepScore));
//			esdcRepScore = getEsdcScore(repScores[0], repScores[1]);
//			sb.append('\t');
//			sb.append(String.format("%.4f", esdcRepScore));
//		}
//		
////		if (args.isPermutate()) {
////			int times = 1000;
////			double kstPer = 0.0;
////			double esdcPer = 0.0;
////			double kstSharePer = 0.0;
////			double esdcSharePer = 0.0;
////			double[] diffScores = filtScoreFill(Method.smoothNoNegArray(scoreFilter(
////					Method.subNoNegAbsArray(scores1, scores2), subSeq), 2), subSeq);
////			double[] diffShareScores = filtScoreFill(CommonMethod.smoothNoNegArray(scoreFilter(
////					CommonMethod.subNoNegAbsArray(scoresShare1, scoresShare2), subSeq), 2), subSeq);
////			double[] diffPer = new double[diffScores.length];
////			double[] diffSmoothPer = new double[diffScores.length];
////			double[] diffSharePer = new double[diffShareScores.length];
////			double[] diffShareSmoothPer = new double[diffShareScores.length];
////			boolean negFlag = Method.retainNoNeg(scores1).length < 2 || Method.retainNoNeg(scores2).length < 2;
////			boolean negShareFlag = CommonMethod.retainNoNeg(scoresShare1).length < 2 || 
////					CommonMethod.retainNoNeg(scoresShare2).length < 2;
////			
////			int kstTime = 0;
////			int esdcTime = 0;
////			permutate = true;
////			for (int i = 0; i < times; ++i) {
////				double[][] testScores = getTestSnpScores(snpNode, genomeSeq, start, end, positive);
////				double testScore = getKstScore(testScores[0], testScores[1]);
////				if (testScore >= 0.0) {
////					kstPer += testScore <= kstScore ? 1.0 : 0.0;
////					kstSharePer += testScore <= kstShareScore ? 1.0 : 0.0;
////					++kstTime;
////				}
//				testScore = getEsdcScore(testScores[0], testScores[1]);
//				if (testScore >= 0.0) {
////					esdcPer += testScore >= esdcScore ? 1.0 : 0.0;
//					esdcSharePer += testScore >= esdcShareScore ? 1.0 : 0.0;
//					++esdcTime;
//				}
//				double[] diffScoresPer = CommonMethod.subNoNegAbsArray(testScores[0], testScores[1]);
//		 		for (int j = 0; j < diffScoresPer.length; ++j) {
////		 			diffPer[j] += diffScoresPer[j] >= diffScores[j] ? 1.0 : 0.0;
//		 			diffSharePer[j] += diffScoresPer[j] >= diffShareScores[j] ? 1.0 : 0.0;
//		 		}
//		 		diffScoresPer = filtScoreFill(CommonMethod.smoothNoNegArray(scoreFilter(diffScoresPer, subSeq), 2), subSeq);
//		 		for (int j = 0; j < diffScoresPer.length; ++j) {
////		 			diffSmoothPer[j] += diffScoresPer[j] >= diffScores[j] ? 1.0 : 0.0;
//		 			diffShareSmoothPer[j] += diffScoresPer[j] >= diffShareScores[j] ? 1.0 : 0.0;
//		 		}
////			}
////			permutate = false;
////			sb.append('\t');
////			sb.append(String.format("%.4f", negFlag || kstTime == 0 ? -1.0 : kstPer / kstTime));
////			sb.append('\t');
////			sb.append(String.format("%.4f", negFlag || esdcTime == 0 ? -1.0 : esdcPer / esdcTime));
//			sb.append('\t');
//			sb.append(String.format("%.4f", negShareFlag || kstTime == 0 ? -1.0 : kstSharePer / kstTime));
//			sb.append('\t');
//			sb.append(String.format("%.4f", negShareFlag || esdcTime == 0 ? -1.0 : esdcSharePer / esdcTime));
//			sb.append('\t');
//			esdcSharePer = esdcSharePer < 0.5 ? Math.log10(times) + 1 : -Math.log10(esdcSharePer / esdcTime);
//			esdcShareScore = negShareFlag || esdcTime == 0 ? -1.0 : 
//				esdcShareScore * esdcSharePer * (1 - Math.log10(snpCount / (treatReads + controlReads)));
//			sb.append('\t');
//			sb.append(String.format("%.4f", Math.abs(esdcShareScore)));
////			sb.append('\t');
////			sb.append(Method.toString(diffPer, "%.0f", ","));
//			sb.append('\t');
//			sb.append(Method.toString(diffSharePer, "%.0f", ","));
////			sb.append('\t');
////			sb.append(Method.toString(diffSmoothPer, "%.0f", ","));
////			sb.append('\t');
////			sb.append(Method.toString(diffShareSmoothPer, "%.0f", ","));
////			sb.append('\t');
////			sb.append(Method.toString(Method.subNoNegAbsArray(scores1, scores2), "%f", ","));
//			sb.append('\t');
//			sb.append(Method.toString(CommonMethod.subNoNegAbsArray(scoresShare1, scoresShare2), "%f", ","));
////			sb.append('\t');
////			sb.append(Method.toString(diffScores, "%f", ","));
//			sb.append('\t');
//			sb.append(Method.toString(diffShareScores, "%f", ","));
//		}
//		
//		if (replicate) {
//			sb.append('\t');
//			sb.append(Method.toString(repScores[0], "%.3f", ","));
//			sb.append('\t');
//			sb.append(Method.toString(repScores[1], "%.3f", ","));
//		}
//		
////		sb.append('\t');
////		sb.append(Method.toString(scores1, "%.3f", ","));
////		sb.append('\t');
////		sb.append(Method.toString(scores2, "%.3f", ","));
//		sb.append('\t');
//		sb.append(Method.toString(scoresShare1, "%.3f", ","));
//		sb.append('\t');
//		sb.append(Method.toString(scoresShare2, "%.3f", ","));
//		sb.append('\t');
//		sb.append(CommonMethod.retainNoNegCol(scoresShare1, scoresShare2).length);
////		for (int i = 0; i < counts[getBaseIndex(snpNode.getValue().getRef())].length; ++i) {
////			sb.append('\t');
////			sb.append(Method.toString(counts[getBaseIndex(snpNode.getValue().getRef())][i], "%d", ","));
////			sb.append('\t');
////			sb.append(Method.toString(counts[getBaseIndex(snp)][i], "%d", ","));
////		}
		return(new String[] {sb.toString(), sbd.toString()});
	}
	
	private String getOddLine(String chr, Node<Mutation> snpnode) {
		StringBuilder sb = new StringBuilder();
		sb.append(chr);
		sb.append('\t');
		sb.append(snpnode.getEnd());
		sb.append('\t');
		sb.append(snpnode.getValue().getDescription());
		sb.append('\t');
		sb.append(snpnode.getValue().getRef());
		sb.append('\t');
		sb.append(snpnode.getValue().getSNPString());
		sb.append('\t');
		sb.append(snpnode.getValue().getCount(snpnode.getValue().getSNP()[0], true));
		sb.append('|');
		sb.append(snpnode.getValue().getCount(snpnode.getValue().getRef(), true));
		sb.append('\t');
		sb.append(snpnode.getValue().getCount(snpnode.getValue().getSNP()[0], false));
		sb.append('|');
		sb.append(snpnode.getValue().getCount(snpnode.getValue().getRef(), false));
		double oddValue = CommonMethod.calORvalue(snpnode.getValue().getCount(snpnode.getValue().getSNP()[0], false), 
				snpnode.getValue().getCount(snpnode.getValue().getRef(), false), 
				snpnode.getValue().getCount(snpnode.getValue().getSNP()[0], true), 
				snpnode.getValue().getCount(snpnode.getValue().getRef(), true));
		oddValue = snpnode.getValue().getDescription().length() >= 1 
				&& Character.toUpperCase(snpnode.getValue().getDescription().charAt(0)) == 'S' ? 1.0 / oddValue : oddValue;
		sb.append('\t');
		sb.append(String.format("%.4f", oddValue));
		return sb.toString();
	}
	
	private double getTestKstScores(double score, int times, Node<Mutation> snpNode, String seq, int start, int end, boolean upstream) {
		if (times <= 0) {
			return -1.0;
		}
		permutate = true;
		double perScore = 0.0;
		double perTime = 1.0 / times;
		for (int i = 0; i < times; ++i) {
			double[][] testScores = getTestSnpScores(snpNode, seq, start, end, upstream);
			double testScore = CommonMethod.calKSTestNoNeg(testScores[0], testScores[1]);
	 		perScore += testScore >= score ? perTime : 0.0;
		}
 		permutate = false;
		return perScore;
	}
	
	private double getTestEsdcScores(double score, int times, Node<Mutation> snpNode, String seq, int start, int end, boolean upstream) {
		if (times <= 0) {
			return -1.0;
		}
		permutate = true;
		double perScore = 0.0;
		double perTime = 1.0 / times;
		for (int i = 0; i < times; ++i) {
			double[][] testScores = getTestSnpScores(snpNode, seq, start, end, upstream);
			double testScore = getEsdcScore(testScores[0], testScores[1]);
	 		perScore += testScore >= score ? perTime : 0.0;
		}
 		permutate = false;
		return perScore;
	}
	
	private double[] getTestShapeScores(double[] scores, int times, Node<Mutation> snpNode, String seq, int start, int end, boolean upstream) {
		double[] perScore = new double[scores.length];
		if (times <= 0) {
			Arrays.fill(perScore, -1.0);
			return perScore;
		}
		double perTime = 1.0 / times;
		permutate = true;
		for (int i = 0; i < times; ++i) {
			double[][] testScores = getTestSnpScores(snpNode, seq, start, end, upstream);
			testScores[0] = CommonMethod.subAbsArray(testScores[0], testScores[1]);
	 		for (int j = 0; j < scores.length; ++j) {
	 			perScore[j] += testScores[0][j] >= scores[j] ? perTime : 0.0;
	 		}
		}
 		permutate = false;
		return perScore;
	}
	
	private double[][] getTestSnpScores(Node<Mutation> snpNode, String seq, int start, int end, boolean upstream) {
		double[][] scores= new double[2][];
		int[][][] counts = new int[2][args.getControlFiles().isEmpty() ? 2 : 4][end - start + 1];
		String subSeq = seq == null ? null : seq.substring(start - 1, end);
		if (!args.getTreatFiles2().isEmpty() && snpNode == null) {
			IntervalTree<List<String>> treeGT = new IntervalTree<>();
			IntervalTree<List<String>> treeGC = new IntervalTree<>();
			for (Node<List<String>> node : treeT) {
				treeGT.put(node.getStart(), node.getEnd(), new ArrayList<>(node.getValue()));
			}
			for (Node<List<String>> node : treeC) {
				treeGC.put(node.getStart(), node.getEnd(), new ArrayList<>(node.getValue()));
			}
			for (Node<List<String>> node : treeTH) {
				List<String> old = treeGT.put(node.getStart(), node.getEnd(), null);
				old = old == null ? new ArrayList<>() : old;
				old.addAll(node.getValue());
				treeGT.put(node.getStart(), node.getEnd(), old);
			}
			for (Node<List<String>> node : treeCH) {
				List<String> old = treeGT.put(node.getStart(), node.getEnd(), null);
				old = old == null ? new ArrayList<>() : old;
				old.addAll(node.getValue());
				treeGT.put(node.getStart(), node.getEnd(), old);
			}
			getCounts(counts, treeGT, treeGC, snpNode, seq, start, end, upstream);
			scores[0] = getWinShapeScores(counts[0], start, end, subSeq);
			scores[1] = getWinShapeScores(counts[1], start, end, subSeq);
			return scores;
		}
		getCounts(counts, treeT, treeC, snpNode, seq, start, end, upstream);
		scores[0] = getWinShapeScores(counts[0], start, end, subSeq);
		scores[1] = getWinShapeScores(counts[1], start, end, subSeq);
 		return scores;
	}
	
//	private double[] getTryRTstopScores(int[][] counts) {
//		double[] score = new double[counts[0].length];
//		double[] control = new double[counts[0].length];
//		double treatSum = 0.0;
//		double controlSum = 0.0;
//		if (counts.length >= 4) {
//			for (int i = 0; i < counts[0].length; ++i) {
//				score[i] = Math.log(counts[0][i] + 1);
//				treatSum += score[i];
//				control[i] =  Math.log(counts[2][i] + 1);
//				controlSum += control[i];
//			}
//			for (int i = 0; i < score.length; ++i) {
//				score[i] = Math.max(0, score[i] / treatSum - control[i] / controlSum);
//			}
//		}
//		return score;
//	}
	
//	private double[] getTryRTstopScores2(int[][] counts) {
//		double[] score = new double[counts[0].length];
//		if (counts.length >= 4) {
//			for (int i = 0; i < counts[0].length; ++i) { 
//				score[i] = counts[0][i] == 0 ? 0 : (double) counts[0][i] / counts[1][i] - (double) counts[2][i] / counts[3][i];
//			}
//		}
//		return score;
//	}
	
//	private double[] getTryRTstopScores3(int[][] counts) {
//		double[] score = new double[counts[0].length];
//		if (counts.length >= 4) {
//			for (int i = 0; i < counts[0].length; ++i) { 
//				score[i] = counts[0][i] == 0 ? 0 : (double) counts[0][i] / counts[1][i];
//			}
//		}
//		return score;
//	}
	
//	private double[] getTryRTstopScores4(int[][] counts) {
//		double[] score = new double[counts[0].length];
//		int treatMax = 0;
//		int controlMax = 0;
//		if (counts.length >= 4) {
//			for (int i = 0; i < counts[0].length; ++i) { 
//				treatMax = Math.max(counts[0][i], treatMax);
//				controlMax = Math.max(counts[2][i], controlMax);
//			}
//			for (int i = 0; i < counts[0].length; ++i) { 
//				score[i] = counts[0][i] == 0 ? 0.0 : (double) counts[0][i] * controlMax / treatMax / counts[2][i];
//			}
//		}
//		return score;
//	}
	
//	private double[] getTryRTstopScores5(int[][] counts) {
//		double[] score = new double[counts[0].length];
//		if (counts.length >= 4) {
//			for (int i = 0; i < counts[0].length; ++i) { 
//				score[i] = (double) (counts[0][i] - counts[2][i]) / counts[3][i];
//			}
//		}
//		return score;
//	}
	
//	private double[] getRTstopScores(int[][] counts) {
//		double[] facts = normalizeFac(counts);
//		if (counts.length >= 4) {
//			return getRTstopScores(counts[0], counts[1], counts[2], counts[3], facts);
//		}
//		return getRTstopScores(counts[0], counts[1], facts);
//	}
	
	private double getKstScore(double[] scores1, double[] scores2) {
		return scores1 == null || scores2 == null ? -1.0 : CommonMethod.calKSTestNoNeg(scores1, scores2);
	}
	
	private double getEsdcScore(double[] scores1, double[] scores2) {
		if (scores1 == null || scores2 == null) {
			return -1.0;
		}
		double[][] testScores = CommonMethod.retainNoNegCol(scores1, scores2);
		if (testScores[0].length < 2) {
			return -1.0;
		}
		double score = CommonMethod.calPearsonCor(testScores[0], testScores[1]);
		score = Double.isFinite(score) ? (1.0 - score) * Math.sqrt(testScores[0].length) : -1.0;
		return score;
	}
	
	private double[] getShapeScores(int[][] counts, int start, int end, String seq) {
		int[][] countsCopy = new int[counts.length][];
		for (int i = 0; i < counts.length; ++i) {
			countsCopy[i] = Arrays.copyOfRange(counts[i], start, end);
		}
		countsCopy = baseCountFilter(countsCopy, seq);
		double[] facts = normalizeFac(countsCopy);
		if (counts.length >= 4) {
			return getShapeScores(countsCopy[0], countsCopy[1], countsCopy[2], countsCopy[3], facts);
		}
		return getShapeScores(countsCopy[0], countsCopy[1], facts);
	}
	
	private double[] getShapeScores(int[] treatStop, int[] treatCount, int[] controlStop, int[] controlCount, double[] facts) {
		double[] scores = new double[treatStop.length];
		for (int i = 0; i < facts.length; ++i) {
			if (facts[i] <= 0.0) {
				Arrays.fill(scores, -1);
				return scores;
			}
		}
		for (int i = 0; i < treatStop.length; ++i) {
			scores[i] = controlCount[i] <= 0 ? -1.0 : div_fac * (treatStop[i] / facts[0] - sub_fac * controlStop[i] / facts[2]) * facts[3] / controlCount[i] ;
		}
		return scores;
	}
	
	private double[] getShapeScores(int[] stop, int[] count, double[] facts) {
		double[] scores = new double[count.length];
		for (int i = 0; i < count.length; ++i) {
			scores[i] = count[i] <= 0 ? -1.0 : div_fac * stop[i] / facts[0] * facts [1] / count[i] ;
		}
		return scores;
	}
	
//	private void countSNP(HashMap<String, IntervalTree<SiteCount>> mut_table) {
//		final int len = 100;
//		final List<Integer> mut_rate = new ArrayList<>(len + 1);
//		for (int i = 0; i <= len; ++i) {
//			mut_rate.add(0);
//		}
//		mut_table.entrySet().forEach(entry -> {
//			for (Node<SiteCount> site_node: entry.getValue()) {
//				SiteCount count = site_node.getValue();
//				if (count.getAllCount() >= 10) {
//					int index = count.getCount(count.getMostChar()) * len / count.getAllCount();
//					mut_rate.set(index, mut_rate.get(index) + 1);
//				}
//			}
//		});
//		try {
//			Method.writeFile(InParam.getParams().getOutPrefix() + "_mut_rate.txt", mut_rate, null);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//	}
	
	private int getBaseIndex(char base) {
		switch (Character.toUpperCase(base)) {
		case 'A':
			return 1;
			
		case 'T':
			return 2;
			
		case 'G':
			return 3;
			
		case 'C':
			return 4;

		default:
			return 0;
		}
	}
	
	private char getChar(int index) {
		switch (index) {
		case 1:
			return 'A';

		case 2:
			return 'T';
			
		case 3:
			return 'G';
			
		case 4:
			return 'C';
			
		case 0:
			return 'N';
			
		default:
			return 'N';
		}
	}
	
	private int getClipStat(SAMRecord record) {
		if (record == null || record.getCigarLength() <= 1) {
			return 0;
		}
		int stat = record.getCigar().isLeftClipped() ? 1 : 0;
		stat += record.getCigar().isRightClipped() ? 2 : 0;
		return stat;
	}
	
	private int getSnpStat(char segChar, char ref, char alt) {
		return segChar == ref ? 1 : segChar == alt ? 2 : 0;
	}
	
	private String getChrSeq(Genome genome, String chr) {
		return genome == null || genome.getChr(chr) == null ? null : genome.getChr(chr).getSeq();
	}
}
