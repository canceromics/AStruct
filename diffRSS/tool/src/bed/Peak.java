package bed;

import genome.IntRegion;
import genome.Transcript;

public class Peak extends IntRegion {
	
	private String description = null;
	private Transcript script = null;
	private Boolean near_circ = null;

	public Peak(int pos, String description) {
		super(pos - 1, pos);
		this.description = description;
	}
	
	public boolean isCirc() {
		return near_circ != null;
	}
	
	public boolean isNearCirc() {
		return near_circ != null && near_circ;
	}
	
	public void setCirc(boolean circ) {
		this.near_circ = circ;
	}

	public Transcript getScript() {
		return script;
	}
	
	public void addScript(Transcript script) {
		if (this.script == null) {
			this.script = script;
		}
		else {
			this.script.addGeneScript(script);
		}
	}
	
	public static String getHeader() {
		return "#Chr\tStart\tEnd\tCirc";
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(super.toString());
		if (near_circ != null) {
			sb.append('\t');
			sb.append(near_circ);
		}
		if (description != null) {
			sb.append('\t');
			sb.append(description);
		}
		if (script != null) {
			sb.append('\t');
			sb.append(script.getScript_id());
			for (Transcript s = script.getGeneScript(); s != null; s = s.getGeneScript()) {
				sb.append(',');
				sb.append(s.getScript_id());
			}
		}
		return sb.toString();
	}

}
