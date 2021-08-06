package bed;

import java.util.ArrayList;
import java.util.List;

import genome.Transcript;
import main.CommonMethod;
import main.InParam;
import mapping.SiteCount;

public class Mutation {

	private static final char[] BLANK = {'N'};
//	private static final char[] BASES = {'A', 'T', 'G', 'C'};
	private char[] mut = null;
	private char ref = 'N';
	private String description = null;
	private SiteCount counts = null;
	private SiteCount counts_ip = null;
	private Boolean allele_mut = null;
	private Transcript script = null;
	private double[] scores = null;
	private double[] altScores = null;
	private double[] refScores = null;
	private double[] altTestScores = null;
	private double[] refTestScores = null;
	private double[] altNoneScores = null;
	private double[] refNoneScores = null;
	private double ref_maf = 0.0;
	private boolean needRemove = false;
	
	public Mutation(char ref, char[] mut, String des) {
		super();
		this.ref = ref;
		this.mut = mut;
		this.description = des;
	}

	public char getRef() {
		return ref;
	}
	
	public char[] getSNP() {
		return mut == null || mut.length == 0 ? BLANK : mut;
	}
	
	public boolean isMut(char c) {
		if (mut == null) {
			return false;
		}
		for (char a : mut) {
			if (a == c) {
				return true;
			}
		}
		return false;
	}
	
	public String getSNPString() {
		StringBuilder sb = new StringBuilder();
		for (char c : mut) {
			sb.append(c);
			sb.append(',');
		}
		if (sb.length() > 0) {
			sb.setLength(sb.length() - 1);
		}
		return sb.toString();
	}
	
	public int getInputCount(char c) {
		return counts == null ? 0 : counts.getCount(c);
	}
	
	public int getIpCount(char c) {
		return counts_ip == null ? 0 : counts_ip.getCount(c);
	}
	
	public int getRefCount(boolean ip) {
		return ip ? counts_ip.getCount(ref) : counts.getCount(ref);
	}
	
	public int getMutCount(boolean ip) {
		int out = 0;
		for (char c : mut) {
			out += ip ? getIpCount(c) : getInputCount(c);
		}
		return out;
	}
	
	public int getCount(char c, boolean ip) {
		return ip ? getIpCount(c) : getInputCount(c);
	}
	
	public int[] getBaseCounts(boolean ip) {
		return ip ? counts_ip.getBaseCounts() : counts.getBaseCounts();
	}
	
	public void incCount(char c, boolean ip) {
		if (ip) {
			if (counts_ip == null) {
				counts_ip = new SiteCount();
			}
			counts_ip.incCount(c);
		}
		else {
			if (counts == null) {
				counts = new SiteCount();
			}
			counts.incCount(c);
		}
	}
	
	public String getDescription() {
		return description;
	}

	public void setDescription(String description) {
		this.description = description;
	}
	
	public boolean isMutMajor() {
		return allele_mut == null ? allele_mut = getMutCount(false) >= getInputCount(ref) : allele_mut;
	}
	
	public boolean hasAlt(boolean ip) {
		return mut != null && getMutCount(ip) > 0;
	}
	
	public boolean hasRef(boolean ip) {
		return getCount(ref, ip) > 0;
	}
	
	public boolean isBothAllele(boolean ip) {
		return hasAlt(ip) && hasRef(ip); 
	}
	
	public boolean isValidSNP(boolean ip, int thr) {
		if (getCount(ref, ip) < thr || mut == null || getCount(ref, true) <= 0) {
			return false;
		}
		for (char c : mut) {
			if (getCount(c, ip) >= thr && getCount(c, true) > 0) {
				return isMostRS(ip);
			}
		}
		return false;
	}
	
	public boolean isMostRS(boolean ip) {
		char most = ip ? counts_ip.getMostChar() : counts.getMostChar();
		if (most == ref) {
			return true;
		}
		for (char c : mut) {
			if (most == c) {
				return true;
			}
		}
		return false;
	}
	
	public void setNeedRemove() {
		needRemove = true;
	}
	
	public boolean needRemove() {
		return needRemove;
	}
	
	public boolean remakeMostAlt(boolean ip) {
		if (mut == null || mut.length == 0) {
			return false;
		}
		if (mut.length == 1) {
			return true;
		}
		List<Character> mostSNP = new ArrayList<>();
		char mostChar = mut[0];
		mostSNP.add(mut[0]);
		for (int i = 1; i < mut.length; ++i) {
			if (getCount(mut[i], ip) > getCount(mostChar, ip)) {
				mostChar = mut[i];
				mostSNP.clear();
				mostSNP.add(mut[i]);
			}
			else if (getCount(mut[i], ip) == getCount(mostChar, ip)) {
				if (!ip && getCount(mut[i], true) > getCount(mostChar, true)) {
					mostChar = mut[i];
					mostSNP = new ArrayList<>();
					mostSNP.add(mut[i]);
				}
				else if (ip || getCount(mut[i], true) == getCount(mostChar, true)) {
					mostSNP.add(mut[i]);
				}
			}
		}
		mut = new char[mostSNP.size()];
		for (int i = 0; i < mut.length; ++i) {
			mut[i] = mostSNP.get(i);
		}
		return true;
	}
	
	public boolean remakeRef(boolean ip) {
		SiteCount count = ip ? this.counts_ip : this.counts;
		if (mut == null && count != null) {
			char ref = count.getMostChar();
			if (count.getCount(ref) * 2 >= count.getAllCount()) {
				return remakeRef(ref, count);
			}
		}
		return true;
	}
	
	private boolean remakeRef(char c, SiteCount count) {
		if (c == ref) {
			remakeMostAlt(false);
		}
		char mut_base = ref;
		int ref_index = -1;
		for (int i = 0; i < mut.length; ++i) {
			if (mut[i] == c) {
				ref_index = i;
			}
			else if (count.getCount(mut_base) < count.getCount(mut[i])) {
				mut_base = mut[i];
			}
		}
		if (ref_index < 0) {
			return false;
		}
		ref = c;
		mut = new char[] {mut_base};
		return true;
	}
	
	public int getMajorAlleleCount(boolean ip) {
		return isMutMajor() ? getMutCount(ip) : getCount(ref, ip);
	}
	
	public int getTotalCount(boolean ip) {
		return ip ? (counts_ip == null ? 0 : counts_ip.getAllCount()) : counts == null ? 0 : counts.getAllCount();
	}
	
	public void setMaf(double ref_maf) {
		this.ref_maf = ref_maf;
	}
	
	public double getMaf() {
		return ref_maf;
	}
	
	public double[] getScores() {
		return scores;
	}

	public void setScores(double[] scores) {
		this.scores = scores;
	}

	public double[] getAltScores() {
		return altScores;
	}

	public void setAltScores(double[] altScores) {
		this.altScores = altScores;
	}

	public double[] getRefScores() {
		return refScores;
	}

	public void setRefScores(double[] refScores) {
		this.refScores = refScores;
	}
	
	public double[] getAltNoneScores() {
		return altNoneScores;
	}

	public double[] getAltTestScores() {
		return altTestScores;
	}

	public void setAltTestScores(double[] altTestScores) {
		this.altTestScores = altTestScores;
	}

	public double[] getRefTestScores() {
		return refTestScores;
	}

	public void setRefTestScores(double[] refTestScores) {
		this.refTestScores = refTestScores;
	}

	public void setAltNoneScores(double[] altNoneScores) {
		this.altNoneScores = altNoneScores;
	}

	public double[] getRefNoneScores() {
		return refNoneScores;
	}

	public void setRefNoneScores(double[] refNoneScores) {
		this.refNoneScores = refNoneScores;
	}

	public Boolean getAllele() {
		return allele_mut;
	}

	public void setAllele(Boolean allele_mut) {
		this.allele_mut = allele_mut;
	}

	public Transcript getScript() {
		return script;
	}
	
	public void setScript(Transcript script) {
		this.script = script;
	}
	
	public static String getSimpleHeader() {
		return 
//	"#Chr\tStart\tEnd\tRef\tAlt\tID\tTag\tAllSNP\tControl\tTreat\tGeneID\tTranscriptID\tTranscriptType\tTranscriptOffset\tRegion\tUpStream\tDownStream\tGCContent";
	"#Chr\tStart\tEnd\tRef\tAlt\tID\tpValue\tFDR\tOR\tTag\tScores\tAltScores\tRefScores\tControl\tTreat\tGeneID\tStrand\tTranscriptID\tTranscriptType\tTranscriptOffset\tRegion\tUpStream\tDownStream\tGCContent";
	}
	
	public static String getHeader() {
		return InParam.getParams().getControlFile() == null ?
				"#Chr\tStart\tEnd\tRef\tAlt\tID\tTag\tAllSNP\tA\tT\tG\tC"
				: "#Chr\tStart\tEnd\tRef\tAlt\tID\tTag\tAllSNP\tAin\tTin\tGin\tCin\tA\tT\tG\tC";
//		return "#Chr\tStart\tEnd\tRef\tAlt\tRefMAF\tMAF\tRealMAF\tRefStat\tStat\tAltCount\tRefCount\tIPAltCount\tIPRefCount";
	}
	
	private String getTag(double maf, int ip_a1, int ip_a2) {
		double ip_maf = (double) ip_a1 / (double) (ip_a1 + ip_a2);
		if (Double.isNaN(ip_maf) || Double.isNaN(maf)){
			return "none";
		}
		else if (CommonMethod.calPvalue(maf, ip_a1, ip_a2) < 0.05) {
			if (ip_maf < maf) {
				return "loss";
			}
			else {
				return "gain";
			}
		}
		else {
			return "unbiased";
		}
	}
	
	private String getTag(double maf, int a1, int a2, int ip_a1, int ip_a2) {
		double ip_maf = (double) ip_a1 / (double) (ip_a1 + ip_a2);
		if (Double.isNaN(ip_maf) || Double.isNaN(maf)){
			return "none";
		}
		else if (CommonMethod.calPvalue(a1, a2, ip_a1, ip_a2) < 0.05) {
			if (ip_maf < maf) {
				return "loss";
			}
			else {
				return "gain";
			}
		}
		else {
			return "unbiased";
		}
	}
	
	private String getTag(double or, double p) {
		if (Double.isNaN(or)){
			return "none";
		}
		else if (p < 0.05) {
			if (or < 1) {
				return "loss";
			}
			else {
				return "gain";
			}
		}
		else {
			return "unbiased";
		}
	}
	
//	private void allBasesTag(StringBuffer sb) {
//		for (char c : BASES) {
//			sb.append(String.format("%.4f", Method.calPvalue(0.5, getIpCount(c), getIpCount(ref))));
//			sb.append('|');
//			sb.append(getTag(0.5, getIpCount(c), getIpCount(ref)));
//			sb.append(';');
//		}
//		sb.setLength(sb.length() - 1);
//	}
	
	private void varientTag(StringBuilder sb, char c) {
		sb.append(c);
		sb.append(':');
		if (InParam.getParams().getControlFile() == null) {
			sb.append(String.format("%.4f", CommonMethod.calPvalue(0.5, getIpCount(c), getIpCount(ref))));
			sb.append('|');
			sb.append(String.format("%.4f", CommonMethod.calORvalue(1, 1, getIpCount(c), getIpCount(ref))));
			sb.append('|');
			sb.append(getTag(0.5, getIpCount(c), getIpCount(ref)));
		}
		else {
			sb.append(String.format("%.4f", CommonMethod.calPvalue(getInputCount(c), getInputCount(ref),
					getIpCount(c), getIpCount(ref))));
			sb.append('|');
			sb.append(String.format("%.4f", CommonMethod.calORvalue(getInputCount(c), getInputCount(ref),
					getIpCount(c), getIpCount(ref))));
			sb.append('|');
			double maf = script == null ? 0.5 : script.getGene().getSimMaf();
			sb.append(getTag(isMutMajor() ? maf : 1.0 - maf, getInputCount(c), getInputCount(ref), 
					getIpCount(c), getIpCount(ref)));
		}
	}
	
	private void allVarientTag(StringBuilder sb) {
		for (char c : mut) {
			varientTag(sb, c);
			sb.append(';');
		}
		sb.setLength(sb.length() - 1);
	}
	
	public String toSimpleString() {
		StringBuilder sb = new StringBuilder();
		sb.append(ref);
		sb.append('\t');
		sb.append(getSNPString());
		sb.append('\t');
		if (description != null) {
			sb.append(description);
		}
		else {
			sb.append('.');
		}
		sb.append('\t');
		char most_mut = mut[0];
		for (int i = 1; i < mut.length; ++i) {
			most_mut = getInputCount(most_mut) + getIpCount(most_mut) < getIpCount(mut[i]) + getInputCount(mut[i]) ?
				mut[i] : most_mut;
		}
		varientTag(sb, most_mut);
		sb.append('\t');
		allVarientTag(sb);
		if (InParam.getParams().getControlFile() != null) {
			sb.append('\t');
			sb.append(getInputCount('A'));
			sb.append('|');
			sb.append(getInputCount('T'));
			sb.append('|');
			sb.append(getInputCount('G'));
			sb.append('|');
			sb.append(getInputCount('C'));
		}
		sb.append('\t');
		sb.append(getIpCount('A'));
		sb.append('|');
		sb.append(getIpCount('T'));
		sb.append('|');
		sb.append(getIpCount('G'));
		sb.append('|');
		sb.append(getIpCount('C'));
		if (script == null) {
			sb.append("\tNone\t.\tNone\tNone");
		}
		else {
			sb.append('\t');
			sb.append(script.getGene().getGene_id());
			sb.append('\t');
			sb.append(script.getGene().getStrand());
			sb.append('\t');
			sb.append(script.getScript_id());
			sb.append('\t');
			sb.append(script.getScript_type());
		}
		return sb.toString();
	}
	
	public String toSimpleString(double p_value, double adjust_p) {
		StringBuilder sb = new StringBuilder();
		sb.append(ref);
		sb.append('\t');
		char most_mut = mut[0];
		for (int i = 1; i < mut.length; ++i) {
			most_mut = getInputCount(most_mut) + getIpCount(most_mut) < getIpCount(mut[i]) + getInputCount(mut[i]) ?
				mut[i] : most_mut;
		}
		sb.append(getSNPString());
		sb.append('\t');
		if (description != null) {
			sb.append(description);
		}
		else {
			sb.append('.');
		}
		sb.append('\t');
		sb.append(String.format("%.4f", p_value));
		sb.append('\t');
		sb.append(String.format("%.4f", adjust_p));
		sb.append('\t');
		double or_value = CommonMethod.calORvalue(getInputCount(most_mut), getInputCount(ref),
				getIpCount(most_mut), getIpCount(ref));
		sb.append(String.format("%.4f", or_value));
		sb.append('\t');
		sb.append(getTag(or_value, adjust_p));
		if (scores != null && scores.length > 0) {
			sb.append('\t');
			for (double score : scores) {
				sb.append(String.format("%.3f", score));
				sb.append('|');
			}
			sb.setLength(sb.length() - 1);
		}
		if (altScores != null && altScores.length > 0) {
			sb.append('\t');
			for (double score : altScores) {
				sb.append(String.format("%.3f", score));
				sb.append('|');
			}
			sb.setLength(sb.length() - 1);
		}
		if (altTestScores != null && altTestScores.length > 0) {
			sb.append('\t');
			for (double score : altTestScores) {
				sb.append(String.format("%.3f", score));
				sb.append('|');
			}
			sb.setLength(sb.length() - 1);
		}
		if (refScores != null && refScores.length > 0) {
			sb.append('\t');
			for (double score : refScores) {
				sb.append(String.format("%.3f", score));
				sb.append('|');
			}
			sb.setLength(sb.length() - 1);
		}
		if (refTestScores != null && refTestScores.length > 0) {
			sb.append('\t');
			for (double score : refTestScores) {
				sb.append(String.format("%.3f", score));
				sb.append('|');
			}
			sb.setLength(sb.length() - 1);
		}
		if (altNoneScores != null && altNoneScores.length > 0) {
			sb.append('\t');
			for (double score : altNoneScores) {
				sb.append(String.format("%.3f", score));
				sb.append('|');
			}
			sb.setLength(sb.length() - 1);
		}
		if (refNoneScores != null && refNoneScores.length > 0) {
			sb.append('\t');
			for (double score : refNoneScores) {
				sb.append(String.format("%.3f", score));
				sb.append('|');
			}
			sb.setLength(sb.length() - 1);
		}
		if (InParam.getParams().getControlFile() != null) {
			sb.append('\t');
			sb.append(getInputCount('A'));
			sb.append('|');
			sb.append(getInputCount('T'));
			sb.append('|');
			sb.append(getInputCount('G'));
			sb.append('|');
			sb.append(getInputCount('C'));
		}
		sb.append('\t');
		sb.append(getIpCount('A'));
		sb.append('|');
		sb.append(getIpCount('T'));
		sb.append('|');
		sb.append(getIpCount('G'));
		sb.append('|');
		sb.append(getIpCount('C'));
		if (script == null) {
			sb.append("\t.\tNone\tNone\tNone");
		}
		else {
			sb.append('\t');
			sb.append(script.getGene().getStrand());
			sb.append('\t');
			sb.append(script.getGene().getGene_id());
			sb.append('\t');
			sb.append(script.getScript_id());
			sb.append('\t');
			sb.append(script.getScript_type());
		}
		return sb.toString();
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(ref);
		sb.append('\t');
		sb.append(getSNPString());
		if (description != null) {
			sb.append('\t');
			sb.append(description);
		}
		int ip_mut_count = getMutCount(true), ip_ref_count = getIpCount(ref),
				mut_count = getMutCount(false), ref_count = getInputCount(ref);
		double maf = script == null ? 0.5 : script.getGene().getSimMaf();
		sb.append('\t');
		sb.append(String.format("%.3f", maf));
		sb.append('\t');
		sb.append(String.format("%.3f", ref_maf));
		sb.append('\t');
		double ip_maf = (double) ip_mut_count / (double) (ip_ref_count + ip_mut_count);
		sb.append(String.format("%.3f", ip_maf));
		sb.append('\t');
		sb.append(getTag(maf, ip_mut_count, ip_ref_count));
		sb.append('\t');
		sb.append(getTag(ref_maf, ip_mut_count, ip_ref_count));
		sb.append('\t');
		sb.append(mut_count);
		sb.append('\t');
		sb.append(ref_count);
		sb.append('\t');
		sb.append(ip_mut_count);
		sb.append('\t');
		sb.append(ip_ref_count);	
		sb.append('\t');
		sb.append(script);
		sb.append('\t');
		maf = script == null ? 0.5 : script.getGene().getMAF();
		sb.append(String.format("%.3f", maf));
		sb.append('\t');
		sb.append(String.format("%.4f", CommonMethod.calPvalue(isMutMajor() ? maf : 1.0 - maf, 
				ip_mut_count, ip_ref_count)));
		sb.append('\t');
		sb.append(getTag(isMutMajor() ? maf : 1.0 - maf, 
				ip_mut_count, ip_ref_count));
		sb.append('\t');
		sb.append(String.format("%.4f", CommonMethod.calPvalue(mut_count, ref_count,
				ip_mut_count, ip_ref_count)));
		sb.append('\t');
		sb.append(getTag(isMutMajor() ? maf : 1.0 - maf, mut_count, ref_count,
				ip_mut_count, ip_ref_count));
		sb.append('\t');
		sb.append(String.format("%.4f", CommonMethod.calPvalue((double) ip_mut_count / (double) (ip_mut_count + ip_ref_count),
				(int) Math.round(getTotalCount(false) * (isMutMajor() ? maf : 1.0 - maf)), 
				(int) Math.round(getTotalCount(false) * (isMutMajor() ? 1.0 - maf : maf)))));
		sb.append('\t');
		sb.append(getTag((double) ip_mut_count / (double) (ip_mut_count + ip_ref_count), 
				(int) Math.round(getTotalCount(false) * (isMutMajor() ? maf : 1.0 - maf)), 
				(int) Math.round(getTotalCount(false) * (isMutMajor() ? 1.0 - maf : maf))));
		maf = (double) mut_count / (double) (mut_count + ref_count);
		sb.append('\t');
		sb.append(String.format("%.3f", maf));
		sb.append('\t');
		sb.append(String.format("%.4f", CommonMethod.calPvalue(maf, 
				ip_mut_count, ip_ref_count)));
		sb.append('\t');
		sb.append(getTag(maf, 
				ip_mut_count, ip_ref_count));
		sb.append('\t');
		sb.append(String.format("%.4f", CommonMethod.calPvalue((double) ip_mut_count / (double) (ip_mut_count + ip_ref_count),
				(int) Math.round(getTotalCount(false) * maf ), 
				(int) Math.round(getTotalCount(false) * (1.0 - maf)))));
		sb.append('\t');
		sb.append(getTag((double) ip_mut_count / (double) (ip_mut_count + ip_ref_count), 
				(int) Math.round(getTotalCount(false) * maf), 
				(int) Math.round(getTotalCount(false) * (1.0 - maf))));
		return sb.toString();
	}
}
