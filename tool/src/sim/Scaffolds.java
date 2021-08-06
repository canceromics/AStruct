package sim;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import genome.IntRegion;
import main.CommonMethod;

public class Scaffolds extends IntRegion{
	
	public class Scaffold {

		private String sitePair = null;
		private int pairSiteNum = 0;
		private double weight = 0.0;
		
		private Scaffold(String sitePair, double weight) {
			setSitePair(sitePair);
			this.weight = weight;
		}
		
		public List<Integer> getPairSites() {
			List<Integer> out = new ArrayList<>();
			for (int i = 0; i < sitePair.length(); ++i) {
				if (sitePair.charAt(i) != '.') {
					out.add(i + getStart() + 1);
				}
			}
			return out;
		}
		
		public List<Integer> getUnpairSites() {
			List<Integer> out = new ArrayList<>();
			for (int i = 0; i < sitePair.length(); ++i) {
				if (sitePair.charAt(i) == '.') {
					out.add(i + getStart() + 1);
				}
			}
			return out;
		}
		
		public boolean isSitePair(int site) {
			return sitePair.charAt(site - getStart() - 1) != '.';
		}
		
		private void setSitePair(String sitePair) {
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < getLength(); ++i) {
				char c = sitePair != null && sitePair.length() > i ? sitePair.charAt(i) : '.';
				sb.append(c);
				pairSiteNum += c == '.' ? 0 : 1;
			}
			this.sitePair = sb.toString();
		}
		
		public String getSitePair() {
			return sitePair;
		}
		
		public int getPairSiteNum() {
			return pairSiteNum;
		}
		
		public double getWeight() {
			return weight;
		}

	}
	
	private Map<String, Scaffold> scarMap = null;
	private String seq = null;
	private double totalWeight = 0.0;
	
	public Scaffolds(int start, int end) {
		super(start, end);
		this.scarMap = new HashMap<>();
	}
	
	public void putScaf(String pairSites) {
		putScaf(pairSites, 1.0);
	}
	
	public void putScaf(String pairSites, double d) {
		Scaffold scaf = new Scaffold(pairSites, d);
		Scaffold old = scarMap.put(pairSites, scaf);
		totalWeight += d - (old == null ? 0.0 : old.getWeight());
	}
	
	public int size() {
		return scarMap.size();
	}
	
	public Scaffold getRandScaf() {
		return CommonMethod.randElement(scarMap.values());
	}
	
	public Scaffold[] getRandScafPair() {
		Scaffold[] scafPair = new Scaffold[2];
		int i = -1;
		for (Iterator<Scaffold> scafs = CommonMethod.randNElement(scarMap.values(), 2).iterator();
				scafs.hasNext();) {
			Scaffold scaf = scafs.next();
			scafPair[++i] = scaf;
		}		
		if (scafPair[1] == null) {
			scafPair[1] = scafPair[0];
		}
		return scafPair;
	}
	
	public Scaffold[] getDiffScafPair() {
		Scaffold[] scafPair = new Scaffold[2];
		int diff = -1;
		List<String> scafs = new ArrayList<>(scarMap.keySet());
		for (int i = 0; i < scafs.size(); ++i) {
			for (int j = i; j < scafs.size(); ++j) {
				int dist = CommonMethod.getHammingDistance(scafs.get(i), scafs.get(j));
				if (diff < dist) {
					scafPair[0] = scarMap.get(scafs.get(i));
					scafPair[1] = scarMap.get(scafs.get(j));
					diff = dist;
				}
			}
		}
		return scafPair;
	}
	
	public Collection<Scaffold> getRandScafs(int n) {
		return CommonMethod.randNElement(scarMap.values(), n);
	}
	
	public List<Scaffold> getAllScaf() {
		List<Scaffold> scafs = new ArrayList<>();
		scarMap.forEach((pair, scaf) -> {
			scafs.add(scaf);
		});
		return scafs;
	}
	
	public double getTotalWeight() {
		return totalWeight;
	}

	public void setSeq(String seq) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < getLength(); ++i) {
			sb.append(seq != null && seq.length() > i ? Character.toUpperCase(seq.charAt(i)) : 'N');
		}
		this.seq = sb.toString();
	}
	
	public String getSeq() {
		return seq;
	}
	
}
