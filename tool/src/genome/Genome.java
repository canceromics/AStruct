package genome;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import main.InParam;
import main.Method;

public class Genome {

	private Map<String, Chromosome> chrs = new HashMap<>();
	
	public Map<String, Chromosome> getChrMap(){
		return chrs;
	}
	
	public Chromosome getChr(String key) {
		return chrs.get(key);
	}
	
	public Map<String, Gene> getGeneMap() {
		final HashMap<String, Gene> out = new HashMap<>();
		chrs.forEach((chr_string, chr) -> {
			out.putAll(chr.getGeneMap());
		});
		return out;
	}
	
	public Map<String, Transcript> getScriptMap() {
		final HashMap<String, Transcript> out = new HashMap<>();
		chrs.forEach((chr_string, chr) -> {
			out.putAll(chr.getScriptMap());
		});
		return out;
	}
	
	public String readAnnoteFile(String fileName) throws IOException {
		try (BufferedReader reader = new BufferedReader(new FileReader(fileName))) {
			return readAnnoteFile(reader, null, "", Method.readGeneFile(InParam.getParams().getGeneFile()));
		}
	}
	
	public String readAnnoteFile(BufferedReader reader) throws IOException {
		return readAnnoteFile(reader, null, "", null);
	}
	
	public String readAnnoteFile(BufferedReader reader, String chr) throws IOException {
		return readAnnoteFile(reader, chr, "", null);
	}
	
	public String readAnnoteFile(BufferedReader reader, String chr, String line) throws IOException {
		return readAnnoteFile(reader, chr, line, null);
	}
	
	public String readAnnoteFile(BufferedReader reader,	String chr, String line, Map<String, Gene> gene_set) throws IOException {
		line = readAnnote(reader, chr, line, gene_set);
		if (chrs.get(chr) != null) {
			chrs.get(chr).sort();
			return line;
		}
		if (chr == null) {
			for (Entry<String, Chromosome> entry : chrs.entrySet()) {
				entry.getValue().sort();
			}
		}
		return line;
	}
	
	private String readAnnote(BufferedReader reader, String chr, String line,
			Map<String, Gene> gene_set) throws IOException {
//		Map<String, Gene> gene_set = Method.getInstance().readGeneFile(InParam.getParams().getGeneFile(), true);
		Map<String, Gene> gene_map = new HashMap<>();
		Map<String, Transcript> script_map = new HashMap<>();
		
		while (line != null) {
			line = handleLine(line, script_map, gene_map, gene_set, chr);
			if (line == null) {
				line = reader.readLine();
			}
			else {
				return line;
			}
		}
		return null;
	}
	
	private String handleLine(String line, Map<String, Transcript> script_map, Map<String, Gene> gene_map,
			Map<String, Gene> gene_set, String chr_need) {
		String exon_file = InParam.getParams().getExonFile().toLowerCase();
		if (exon_file.endsWith("gtf")) {
			return handleLineGTF(line, gene_map, script_map, gene_set, chr_need);
		}
		else if (exon_file.endsWith("gff")) {
			return handleLineGFF(line, gene_map, script_map, gene_set, chr_need);
		}
		else {
			return null;
		}
	}

	private String handleLineGTF(String line, Map<String, Gene> gene_map, Map<String, Transcript> script_map, 
			Map<String, Gene> gene_set, String chr_need) {
		if (line == null || line.length() == 0 || line.charAt(0) == '#') {
			return null;
		}
		String[] cols = line.split("\t");
		if (cols.length < 9) {
			System.err.println("Warining: There is less than 9 columns in line: \n" + line);
			return null;
		}
		String chr_string = cols[0].intern();
		if (chr_need != null && chr_need.intern() != chr_string) {
			return line;
		}
		Chromosome chr = getChr(chr_string);
		if (chr == null) {
			chr = new Chromosome(0, chr_string);
			chrs.put(chr_string, chr);
		}
		String gene_id = getIDSym(line, "gene_id", " \"", "\";");
		if (gene_id.length() < 1) {
			return null;
		}
		String gene_symbol = getIDSym(line, "gene_name", " \"", "\";");
		String gene_type = getIDSym(line, "gene_type", " \"", "\";");
		String script_id = getIDSym(line, "transcript_id", " \"", "\";");
		String script_type = getIDSym(line, "transcript_type", " \"", "\";");
		Transcript script = script_map.get(script_id);
		Gene gene = gene_map.get(gene_id);
		int start = Integer.parseInt(cols[3]) - 1;
		int end = Integer.parseInt(cols[4]);
		if ("gene".equals(cols[2])) {
			if (gene != null) {
				System.err.println("Warning: Duplicated gene_id in line: \n" + line);
				return null;
			}
			gene = new Gene(chr, gene_id, gene_symbol, gene_type, cols[6].charAt(0), start, end);
			putGeneInChr(chr, gene, gene_set, gene_map);
			return null;
		}
		if (script_id.length() < 1) {
			return null;
		}
		if ("exon".equals(cols[2]) || "CDS".equals(cols[2])) {
			Exon exon = new Exon(script, start, end);
			if (script == null) {
				script = new Transcript(gene, script_id, script_type, 0, 0);
				exon.setScript(script);
				script_map.put(script_id, script);
				if (gene == null) {
					gene = new Gene(chr, gene_id, gene_symbol, gene_type, cols[6].charAt(0), 0, 0);
					putGeneInChr(chr, gene, gene_set, gene_map);
					script.setGene(gene);
				}
				gene.addScript(script);
			}
			script.addLastExon(cols[2], exon);
		}
		else if (cols[2].contains("transcript") || cols[2].contains("RNA") || cols[2].contains("gene_segment")) {
			script = new Transcript(gene, script_id, script_type, start, end);
			if (gene == null) {
				gene = new Gene(chr, gene_id, gene_symbol, gene_type, cols[6].charAt(0), 0, 0);
				putGeneInChr(chr, gene, gene_set, gene_map);
				script.setGene(gene);
			}
			gene.addScript(script);
			script_map.put(script_id, script);
		}
		return null;
	}
	
	private String handleLineGFF(String line, Map<String, Gene> gene_map, Map<String, Transcript> script_map,
			Map<String, Gene> gene_set, String chr_need) {
		if (line == null || line.length() == 0 || line.charAt(0) == '#') {
			return null;
		}
		String[] cols = line.split("\t");
		if (cols.length < 9) {
			System.err.println("Warining: There is less than 9 columns in line: \n" + line);
			return null;
		}
		String chr_string = cols[0].intern();
		if (chr_need != null && chr_need.intern() != chr_string) {
			return line;
		}
		Chromosome chr = getChr(chr_string);
		if (chr == null) {
			chr = new Chromosome(0, chr_string);
			getChrMap().put(chr_string, chr);
		}
		String id = getIDSym(line, "ID", "=", ";");
		if (id.length() < 1) {
			return null;
		}
		String gene_symbol = getIDSym(line, "gene", "=", ";");
		String gene_type = getIDSym(line, "gene_biotype", "=", ";");
		String parent = getIDSym(line, "Parent", "=", ";");
		Gene gene = gene_map.get(parent);
		Transcript script = script_map.get(parent);
		int start = Integer.parseInt(cols[3]) - 1;
		int end = Integer.parseInt(cols[4]);
		if ("gene".equals(cols[2])) {
			gene = new Gene(chr, id, gene_symbol, gene_type, cols[6].charAt(0), start, end);
			putGeneInChr(chr, gene, gene_set, null);
			gene_map.put(id, gene);
			return null;
		}
		if (parent.length() < 1) {
			return null;
		}
		if (("exon".equals(cols[2]) || "CDS".equals(cols[2]))) {
			if (script != null) {
				Exon exon = new Exon(script, start, end);
				script.addLastExon(cols[2], exon);
			}
		}
		else {
			script = new Transcript(gene, id, cols[2], start, end);
			if (gene == null) {
				gene = new Gene(chr, parent, gene_symbol, gene_type, cols[6].charAt(0), 0, 0);
				putGeneInChr(chr, gene, gene_set, null);
				gene_map.put(parent, gene);
				script.setGene(gene);
			}
			gene.addScript(script);
			script_map.put(id, script);
		}
		return null;
	}
	
	private void putGeneInChr(Chromosome chr, Gene gene, Map<String, Gene> gene_set, Map<String, Gene> gene_map) {
		if (gene_set == null || gene_set.containsKey(gene.getGene_id())) {
			if (gene_map != null) {
				gene_map.put(gene.getGene_id(), gene);
			}
			gene.setAnotherGene(chr.getGeneTree().put(Math.min(gene.getStart() + 1, gene.getEnd()), gene.getEnd(), gene));
			if (gene_set != null && gene_set.get(gene.getGene_id()) != null) {
				gene.setInputReadCount(gene_set.get(gene.getGene_id()).getInputReadCount());
				gene.setIpReadCount(gene_set.get(gene.getGene_id()).getIpReadCount());
				gene.setIp_enrich(gene_set.get(gene.getGene_id()).getIp_enrich());
			}
		}
	}
	
	private String getIDSym(String line, String key, String head, String tail) {
		int index = line.indexOf(key + head);
		if (index != -1) {
			int right_index = line.indexOf(tail, index);
			right_index = right_index == -1 ? line.length() + 1 - head.length() : right_index;
			index += key.length() + head.length();
			return line.substring(index, right_index);
		}
		return "";
	}
	
	public String readSeqFile(String file_name) throws IOException {
		if (file_name == null) {
			return null;
		}
		try (BufferedReader reader = new BufferedReader(new FileReader(file_name))) {
			return readSeqFile(reader, null, null);
		}
	}
	
	public String readSeqFile(BufferedReader reader) throws IOException {
		return readSeqFile(reader, null, null);
	}
	
	public String readSeqFile(BufferedReader reader, String chr_need, String lastLine) throws IOException {
		StringBuilder sb = new StringBuilder();
		String line = null;
		String chr = lastLine != null && lastLine.charAt(0) == '>' ? lastLine.split(" ")[0].substring(1).intern() : null;
		while ((line = reader.readLine()) != null) {
			if (line.length() == 0 || line.charAt(0) == '#') {
				continue;
			}
			if (line.charAt(0) == '>') {
				if (chr_need != null && chr_need.equals(chr) && getChr(chr_need) != null) {
					getChr(chr_need).setSeq(sb.toString());
					return line;
				}
				if (chr != null) {
					if (!getChrMap().containsKey(chr)) {
						getChrMap().put(chr, new Chromosome(0, chr));
					}
					getChr(chr).setSeq(sb.toString());
					sb = new StringBuilder();
				}
				chr = line.split(" ")[0].substring(1).intern();
			}
			else {
				sb.append(line);
			}
		}
		if (chr != null) {
			if (!getChrMap().containsKey(chr)) {
				getChrMap().put(chr, new Chromosome(0, chr));
			}
			getChr(chr).setSeq(sb.toString());
		}
		return line;
	}
}
