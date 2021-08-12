package main;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

import smile.stat.hypothesis.KSTest;

public class CommonMethod {

	static Random rand = new Random();
	private static PearsonsCorrelation pearCor = new PearsonsCorrelation();
	
	private CommonMethod(){
	}
	
	public static int randInt(int n) {
		return rand.nextInt(n);
	}
	
	public static double randDouble() {
		return rand.nextDouble();
	}
	
	public static <T> T randElement(Collection<T> c) {
		if (c == null || c.isEmpty()) {
			return null;
		}
		int index = randInt(c.size()) + 1;
		for (T t : c) {
			if (--index == 0) {
				return t;
			}
		}
		return null;
	}
	
	public static <K, V> K randKeyWithValue(Map<K, V> map, V value) {
		final List<K> keys = new ArrayList<>();
		map.forEach((key, v) -> {
			if (v == value || (value != null && value.equals(v))) {
				keys.add(key);
			}
		});
		return randElement(keys);
	}
	
	public static Set<Integer> randInSize(int l, int k) {
		Set<Integer> out = new HashSet<>();
		if (l <= k) {
			for (int i = 0; i < l; i++) {
				out.add(i);
			}
			return out;
		}
		int[] indexs = new int[l];
		for (int i = 1; i < l; ++i) {
			indexs[i] = i;
		}
		for (int i = 0; i < k; ++i) {
			int index = randInt(l - i) + i;
			out.add(indexs[index]);
			indexs[i] ^= indexs[index];
			indexs[index] ^= indexs[i];
			indexs[i] ^= indexs[index];
		}
		return out;
	}
	
	public static int randWithWeight(Collection<Double> weights) {
		double[] sums = new double[weights.size()];
		int i = 0;
		for (Double d : weights) {
			double pre = i > 0 ? sums[i - 1] : 0.0;
			sums[i] = pre + d;
		}
		int index = Arrays.binarySearch(sums, CommonMethod.randDouble() * sums[weights.size() - 1]);
		index = index < 0 ? - index - 1 : index;
		return index;
	}
	
	public static Set<Integer> randNWithWeight(Collection<Double> weights, int k) {
		Set<Integer> out = new HashSet<>();
		if (k >= weights.size()) {
			for (int i = 0; i < weights.size(); ++i) {
				out.add(i);
			}
			return out;
		}
		double[] allWeight = new double[weights.size()];
		double sum = 0.0;
		int index = -1;
		for (Double d : weights) {
			allWeight[++index] = d;
			sum += d;
		}
		for (int i = 0; i < k; ++i) {
			double rand = sum * randDouble();
			for (index = 0; out.contains(index) || rand < allWeight[index]; ++index) {
				rand -= out.contains(index) ? 0.0 : allWeight[index];
			}
			out.add(index);
			sum -= allWeight[index];
		}
		return out;
	}
	
	public static <T> Collection<T> randNElement(Collection<T> c, int k) {
		if (c == null || c.size() <= k) {
			return c;
		}
		Set<Integer> indexs = randInSize(c.size(), k);
		Collection<T> out = new ArrayList<>(k);
		int index = 0;
		for (T t : c) {
			if (indexs.contains(index)) {
				out.add(t);
			}
			++index;
		}
		return out;
	}

	public static int indexOfKth(String source, String needle, int k) {
		for (int index = source.indexOf(needle); index != -1; index = source.indexOf(needle, index + 1)) {
			if (--k == 0) {
				return index;
			}
		}
		return -1;
	}
	
	public static int getHammingDistance(String s1, String s2) {
		if (s1 == null) {
			return s2 == null ? 0 : s2.length();
		}
		if (s2 == null) {
			return s1.length();
		}
		if (s2.length() < s1.length()) {
			String s = s1;
			s1 = s2;
			s2 = s;
		}
		int dist = s2.length() - s1.length();
		for (int i = 0; i < s1.length(); ++i) {
			dist += s1.charAt(i) == s2.charAt(i) ? 0 : 1;
		}
		return dist;
	}
	
	public static double calORvalue(int in_a1, int in_a2, int ip_a1, int ip_a2) {
		return (double) (ip_a1 * in_a2) / (double) (ip_a2 * in_a1);
	}
	
	public static double calPvalue(int in_a1, int in_a2, int ip_a1, int ip_a2) {
		if (ip_a1 + ip_a2 + in_a1 + in_a2 == 0) {
			return 1.0;
		}
		HypergeometricDistribution hd = new HypergeometricDistribution(ip_a1 + ip_a2 + in_a1 + in_a2, ip_a1 + in_a1, ip_a1 + ip_a2);
		return hd.cumulativeProbability(ip_a1 - 1, ip_a1 + in_a1);
//		return (double) ip_a1 / (double) ip_a2 > (double) in_a1 / (double) in_a2 ? hd.cumulativeProbability(ip_a1 - 1, ip_a1 + in_a1)
//		: hd.cumulativeProbability(-1, ip_a1);
	}
	
	public static double calPvalue(double p, int ip_a1, int ip_a2) {
		BinomialDistribution bd = new BinomialDistribution(ip_a1 + ip_a2, p);
		boolean gain = (double) ip_a1 / (double) (ip_a1 + ip_a2) > p;
		if (gain) {
			return 1.0 - bd.cumulativeProbability(ip_a1 - 1);
		}
		else {
			return bd.cumulativeProbability(ip_a1);
		}
	}
	
	public static void adjustPValue(List<Double> p_value, String method){
		if ("bon".equals(method)) {
			double n = (double) p_value.size();
			for (int i = 0; i < p_value.size() ; ++i) {
				p_value.set(i, p_value.get(i) * n / (i + 1));
			}
		}
		else if ("bh".equals(method)) {
			List<double[]> sort = new ArrayList<>(p_value.size());
			for (int i = 0; i < p_value.size(); ++i) {
				sort.add(new double[] {p_value.get(i), i});
			}
			sort.sort((p1, p2) -> {
				Double d = p1[0];
				return d.compareTo(p2[0]);
			});
			double value = Double.MAX_VALUE;
			for (int i = sort.size(); i > 0 ; --i) {
				value = Math.min(value, sort.get(i - 1)[0] * (double) sort.size() / (double) i);
				p_value.set((int) sort.get(i - 1)[1], value);
			}
		}
		else {
			Method.logGlobal("Warning: unkown method of adjusting p-value, adjusting is disabled");
		}
	}
	
	public static double combinePvalue(double[] p_values) {
		if (p_values == null || p_values.length == 0) {
			return Double.NaN;
		}
		double p_value = 0.0;
		switch (InParam.getParams().getCombineP()) {
		case "ave":
			for (double d : p_values) {
				p_value += d;
			}
			p_value = p_value / p_values.length;
			break;
			
		case "fm":
			for (double d : p_values) {
				p_value += d <= 0.0 ? Math.log(Double.MIN_VALUE) : Math.log(d);
			}
			p_value *= -2.0;
			ChiSquaredDistribution csd = new ChiSquaredDistribution(2.0 * p_values.length);
			p_value = 1.0 - csd.cumulativeProbability(p_value);
			break;
			
		default:
			Method.logGlobal("Warning: Unknown method of combining p-value, using average of -log");
		case "log":
			for (double d : p_values) {
				p_value += d <= 0.0 ? -Math.log(Double.MIN_VALUE) : -Math.log(d);
			}
			p_value = p_value / p_values.length;
			break;
		}
		return p_value;
	}
	
	public static double calAUC(List<double[]> roc) {
		double auc = 0.0;
		double[] last = roc.get(0);
		for (int i = 1; i < roc.size(); ++i) {
			auc += (last[1] + roc.get(i)[1]) / 2 * (roc.get(i)[0] - last[0]);
			last = roc.get(i);
		}
		return auc;
	}
	
	public static List<double[]> calROC(double[] scores, boolean[] labels) {
		int len = Math.min(scores.length, labels.length);
		double[][] data = new double[len][2];
		int positive = 0;
		for (int i = 0; i < len; ++i) {
			data[i][0] = scores[i];
			data[i][1] = labels[i] ? 1.0 : 0.0;
			positive += labels[i] ? 1 : 0;
		}
		return calROC(data, positive, 0.5);
	}
	
	public static List<double[]> calROC(double[] scores, double[] labels, double thrL) {
		int len = Math.min(scores.length, labels.length);
		double[][] data = new double[len][2];
		int positive = 0;
		for (int i = 0; i < len; ++i) {
			data[i][0] = scores[i];
			data[i][1] = labels[i];
			positive += labels[i] >= thrL ? 1 : 0;
		}
		return calROC(data, positive, thrL);
	}
	
	public static List<double[]> calROC(double[][] data, double thrL) {
		int len = data.length;
		int positive = 0;
		for (int i = 0; i < len; ++i) {
			positive += data[i][1] >= thrL ? 1 : 0;
		}
		return calROC(data, positive, thrL);
	}
	
	public static List<double[]> calROC(double[][] data, int positive, double thrL) {
		int len = data.length;
		int usedPositive = 0;
		List<double[]> curve = new ArrayList<>();
		Arrays.sort(data, (o1, o2) -> o1[0] - o2[0] < 0.0 ? 1 : o1[0] - o2[0] > 0.0 ? -1 : 0);
		curve.add(new double[] {0.0, 0.0, data[0][0]});
		int lastDirection = 0;
		for (int i = 0; i < len;) {
			double thr = data[i][0];
			int direction = 3;
			while (i < len && data[i][0] >= thr) {
				direction &= data[i][1] >= thrL ? 1 : 2;
				usedPositive += data[i][1] >= thrL ? 1 : 0;
				++i;
			}
			double tpr = (double) (usedPositive) / (double) positive;
			double fpr = (double) (i - usedPositive) / (double) (len - positive);
			switch (lastDirection & direction) {
			case 0:
				curve.add(new double[] {fpr, tpr, thr});
				break;

			case 1:
				curve.get(curve.size() - 1)[1] = tpr;
				break;
				
			case 2:
				curve.get(curve.size() - 1)[0] = fpr;
				break;
				
			default:
				break;
			}
			lastDirection = direction;
		}
		return curve;
	}
	
	public static List<Integer> lineToList(String line, String sep) {
		List<Integer> out = new ArrayList<>();
		String[] nums = line.split(sep);
		for (String num : nums) {
			out.add(Integer.parseInt(num));
		}
		return out;
	}
	
	public static double pow2(double x) {
		return x * x;
	}
	
	public static double calKSTest(double[] scores1, double[] scores2) {
		if (scores1.length < 2 || scores2.length < 2) {
			return -1.0;
		}
		
		try {
			return KSTest.test(scores1, scores2).pvalue;
		} catch (Exception e) {
			return -2.0;
		}
	}
	
	public static double calKSTestNoNeg(double[] scores1, double[] scores2) {
		double[] testScores1 = retainNoNeg(scores1);
		double[] testScores2 = retainNoNeg(scores2);
		return calKSTest(testScores1, testScores2);
	}
	
	public static double[][] retainNoNegCol(double[] ... arrays) {
		List<List<Double>> tmp = new ArrayList<>();
		for (int i = 0; i < arrays.length; ++i) {
			tmp.add(new ArrayList<>());
		}
		for (int i = 0; i < arrays[0].length; ++i) {
			int j = 0;
			for (; j < arrays.length; ++j) {
				if (arrays[j][i] < 0.0) {
					break;
				}
			}
			if (j == arrays.length) {
				for (j = 0; j < arrays.length; ++j) {
					tmp.get(j).add(arrays[j][i]);
				}
			}
		}
		double[][] noNeg = new double[tmp.size()][tmp.get(0).size()];
		for (int i = 0; i < noNeg.length; ++i) {
			for (int j = 0; j < noNeg[0].length; ++j) {
				noNeg[i][j] = tmp.get(i).get(j);
			}
		}
		return noNeg;
	}
	
	public static double[][] retainNoNeg(double[][] array) {
		List<double[]> tmp = new ArrayList<>();
		for (int i = 0; i < array.length; ++i) {
			int j = 0;
			double[] tmpLine = new double[array[i].length];
			for (; j < array[i].length; ++j) {
				if (array[i][j] < 0.0) {
					break;
				}
				tmpLine[j] = array[i][j];
			}
			if (j == array[i].length) {
				tmp.add(tmpLine);
			}
		}
		double[][] noNeg = new double[tmp.size()][];
		for (int i = 0; i < noNeg.length; ++i) {
			noNeg[i] = tmp.get(i);
		}
		return noNeg;
	}
	
	public static double[] retainNoNeg(double[] array) {
		List<Double> tmp = new ArrayList<>();
		for (int i = 0; i < array.length; ++i) {
			if (array[i] >= 0.0) {
				tmp.add(array[i]);
			}
		}
		double[] noNeg = new double[tmp.size()];
		for (int i = 0; i < noNeg.length; ++i) {
			noNeg[i] = tmp.get(i);
		}
		return noNeg;
	}
	
	public static double[][] arrayTrans(double[] ... arrays) {
		double[][] trans = new double[arrays[0].length][arrays.length];
		for (int i = 0; i < arrays.length; ++i) {
			for (int j = 0; j < arrays[0].length; ++j) {
				trans[j][i] = arrays[i][j];
			}
		}
		return trans;
	}
	
//	public static double[] smoothArray(double[] array, int rad) {
//		double[] smooth = new double[array.length];
//		double sum = 0.0;
//		int len = 2 * rad + 1;
//		if (array.length < len) {
//			return array;
//		}
//		for (int i = 0; i < rad; ++i) {
//			sum += array[i];
//		}
//		for (int i = rad; i <= 2 * rad; ++i) {
//			sum += array[i];
//			smooth[i - rad] = sum / (i + 1);
//		}
//		for (int i = len; i < array.length; ++i) {
//			sum += array[i];
//			sum -= array[i - len];
//			smooth[i - rad] = sum / len;
//		}
//		for (int i = array.length - rad; i < array.length; ++i) {
//			sum -= array[i - rad - 1];
//			smooth[i] = sum / (array.length - i + rad);
//		}
//		return smooth;
//	}
	
	public static double[] smoothNoNegArray(double[] array, int rad) {
		double[] smooth = new double[array.length];
		double sum = 0.0;
		int len = 2 * rad + 1;
		if (array.length < len) {
			return array;
		}
		for (int i = 0; i < len; ++i) {
			sum += Math.max(0.0, array[i]);
			smooth[(i + 1) / 2] = sum / (i + 1);
		}
		for (int i = len; i < array.length; ++i) {
			sum += Math.max(0.0, array[i]);
			sum -= Math.max(0.0, array[i - len]);
			smooth[i - rad] = sum / len;
		}
		for (int i = array.length - rad; i < array.length; ++i) {
			int sub = array.length - i;
			sum -= Math.max(0.0, array[i - sub - 1]);
			sum -= Math.max(0.0, array[i - sub]);
			smooth[i] = sum / (2 * sub - 1);
		}
		for (int i = 0; i < array.length; ++i) {
			smooth[i] = array[i] < 0.0 ? -1.0 : smooth[i];
		}
		return smooth;
	}
	
//	public static double[] smoothNoNegRegion(double[] array, int rad) {
//		double[] smooth = array.clone();
//		int start = 0;
//		int len = 2 * rad + 1;
//		for (int i = 0; i < array.length; ++i) {
//			if (array[i] < 0.0) {
//				if (i - start > len) {
//					double[] copy = Arrays.copyOfRange(array, start, i);
//					copy = smoothNoNegArray(copy, rad);
//					for (int j = start; j < i; ++j) {
//						smooth[j] = copy[j - start];
//					}
//				}
//				start = i + 1;
//			}
//		}
//		if (array.length - start > len) {
//			double[] copy = Arrays.copyOfRange(array, start, array.length);
//			copy = smoothNoNegArray(copy, rad);
//			for (int j = start; j < array.length; ++j) {
//				smooth[j] = copy[j - start];
//			}
//		}
//		return smooth;
//	}
	
	public static boolean equalArray(double[] darray1, double[] darray2) {
		if (darray1 == darray2) {
			return true;
		}
		if (darray1 == null || darray2 == null || darray1.length != darray2.length) {
			return false;
		}
		for (int i = 0; i < darray1.length; ++i) {
			if (darray1[i] != darray2[i]) {
				return false;
			}
		}
		return true;
	}

	public static int[][] addArray(int[][] ... arrays) {
		int[][] sum = new int[arrays[0].length][arrays[0][0].length];
		for (int i = 0; i < arrays.length; ++i) {
			for (int j = 0; j < arrays[0].length; ++j) {
				for (int k = 0; k < arrays[0][0].length; ++k) {
					sum[j][k] += arrays[i][j][k];
				}
			}
		}
		return sum;
	}
	
	public static int[][] addArrayWithout(int[][][] arrays, int index) {
		int[][] sum = new int[arrays[0].length][arrays[0][0].length];
		for (int i = 0; i < arrays.length; ++i) {
			if (i == index) {
				continue;
			}
			for (int j = 0; j < arrays[0].length; ++j) {
				for (int k = 0; k < arrays[0][0].length; ++k) {
					sum[j][k] += arrays[i][j][k];
				}
			}
		}
		return sum;
	}
	
	public static double[] subAbsArray(double[] arrays1, double[] arrays2) {
		int len = Math.max(arrays1.length, arrays2.length);
		double[] sum = new double[len];
		for (int i = 0; i < len; ++i) {
			double d1 = i < arrays1.length ? arrays1[i] : 0.0;
			double d2 = i < arrays2.length ? arrays2[i] : 0.0;
			sum[i] = Math.abs(d1 - d2);
		}
		return sum;
	}
	
	public static double[] subNoNegAbsArray(double[] arrays1, double[] arrays2) {
		int len = Math.max(arrays1.length, arrays2.length);
		double[] sum = new double[len];
		for (int i = 0; i < len; ++i) {
			double d1 = i < arrays1.length ? arrays1[i] : 0.0;
			double d2 = i < arrays2.length ? arrays2[i] : 0.0;
			sum[i] = d1 < 0.0 || d2 < 0.0 ? -1.0 : Math.abs(d1 - d2);
		}
		return sum;
	}
	
	public static boolean equalMatrix(int[][][] array1, int[][][] array2) {
		if (array1 == array2) {
			return true;
		}
		if (array1 == null || array2 == null || array1.length != array2.length) {
			return false;
		}
		for (int i = 0; i < array1.length; ++i) {
			if (array1[i].length != array2[i].length) {
				return false;
			}
			for (int j = 0; j < array1[i].length; ++j) {
				if (array1[i][j].length != array2[i][j].length) {
					return false;
				}
				for (int k = 0; k < array1[i][j].length; ++k) {
					if (array1[i][j][k] != array2[i][j][k]) {
						return false;
					}
				}
			}
		}
		return true;
	}
	
	public static double calPearsonCor(double[] scores1, double[] scores2) {
		if (scores1.length < 2 || scores2.length < 2) {
			return -1.0;
		}
		return pearCor.correlation(scores1, scores2);
	}
	
	public static double calPearsonCorNoNeg(double[] scores1, double[] scores2) {
		double[][] testScores1 = retainNoNegCol(scores1, scores2);
		return pearCor.correlation(testScores1[0], testScores1[1]);
	}
}
