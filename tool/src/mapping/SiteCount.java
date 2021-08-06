package mapping;

import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

public class SiteCount{
	
	private int[] counts = null;
	
	public SiteCount() {
		counts = new int[6];
	}
	
	public void incCount(char c) {
		++counts[getIndex(c)];
	}
	
	public void addCount(char c, int i) {
		counts[getIndex(c)] += i;
	}
	
	public void setCount(char c, int i) {
		counts[getIndex(c)] = i;
	}
	
	public boolean setCountByTotal(char c, int total) {
		int sum = 0;
		for (int i = 0; i < 4; ++i) {
			sum += counts[i];
		}
		if (total < sum) {
			return false;
		}
		counts[getIndex(c)] += total - sum;
		return true;
	}
	
//	public boolean setCountByTotalWithBlank(char c, int total) {
//		int sum = getAllCount();
//		if (total < sum) {
//			return false;
//		}
//		counts[getIndex(c)] += total - sum;
//		return true;
//	}
	
	public int[] getBaseCounts() {
		return Arrays.copyOf(counts, 4);
	}
	
	public int getCount(char c) {
		return counts[getIndex(c)];
	}
	
	public char getMostBase() {
		int max_index = 0;
		for (int i = 1; i < 4; ++i) {
			max_index = counts[i] > counts[max_index] ? i : max_index;
		}
		return getChar(max_index);
	}
	
	public char getMostChar() {
		int max_index = 0;
		for (int i = 1; i < counts.length; ++i) {
			max_index = counts[i] > counts[max_index] ? i : max_index;
		}
		return getChar(max_index);
	}
	
	public int getAllCount() {
		int sum = 0;
		for (int i = 0; i < counts.length; ++i) {
			sum += counts[i];
		}
		return sum;
	}
	
	private int getIndex(char c) {
		switch (c) {
		case 'A':
		case 'a':
			return 0;
			
		case 'T':
		case 't':
			return 1;
			
		case 'G':
		case 'g':
			return 2;
			
		case 'C':
		case 'c':
			return 3;

		case 'D':
		case 'd':
			return 4;
			
		case 'N':
		case 'n':
			return 5;
			
		default:
			Logger.getLogger("debug").log(Level.INFO, "Unknown character: {0}", c);
			return 5;
		}
	}
	
	private char getChar(int index) {
		switch (index) {
		case 0:
			return 'A';

		case 1:
			return 'T';
			
		case 2:
			return 'G';
			
		case 3:
			return 'C';
			
		case 4:
			return 'D';
			
		default:
			return 'N';
		}
	}
}
